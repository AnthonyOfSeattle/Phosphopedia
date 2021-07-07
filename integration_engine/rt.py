import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.isotonic import IsotonicRegression

class AugmentedLinearRegression(LinearRegression):
    def invert(self, X, y):
        return (y-self.intercept_[0])/self.coef_[0][0]
    
    def get_input_grad(self, X, y):
        pred_y = self.predict(X)
        input_grad = -2 * self.coef_[0][0] * (y - pred_y)
        
        return input_grad


class AugmentedIsotonicRegression(IsotonicRegression):
    def invert(self, X, y):
        indices = np.minimum(
            np.searchsorted(self.y_thresholds_, y), 
            len(self.y_thresholds_) - 1
        )
        delta = self.y_thresholds_[indices] - self.X_thresholds_[indices]
        return y - delta
    
    def get_input_grad(self, X, y):
        pred_y = self.predict(X)
        input_grad = -2 * (y - pred_y)
        
        return input_grad


class RTAligner:
    def __init__(self, base_regressor, br_kws={}, method="min", lr=1e-1, max_epochs=10):
        # Store regressor class and any keywords to be initialized later
        self._base_regressor = base_regressor
        self._br_kws = br_kws          
        
        # Create internal algorithm components
        self._method = method
        
        self._lr = lr
        self._max_epochs = max_epochs
        self._history = []
        
        self._sample_names = np.array([])
        self._sample_weights = {}
        self._peptide_df = None
        
        self._target_shape = None
        self._target_dict = {}
        self._select_dict = {}
        self._model_dict = {}
        
    def _initialize_params(self, X, pep_ids, sample_names):
        # initialize peptide data
        print("Initialize peptide data...")
        external_id, count = np.unique(pep_ids, return_counts=True)
        self._peptide_df = pd.DataFrame(dict(external_id = external_id, occurence = count))
        self._peptide_df["internal_id"] = np.arange(self._peptide_df.shape[0])
        self._peptide_df = self._peptide_df.join(
            pd.DataFrame(dict(learned_rt = X.flatten()))\
                        .groupby(pep_ids.flatten())\
                        .mean(),
            on="external_id"
        )
        self._peptide_df["fit_weight"] = 1.
        
        # initialize internals
        print("Initializing internal model targets...")
        self._target_shape = (-1,) if X.ndim == 1 else (-1, 1)
        self._sample_names = np.unique(sample_names)
        
        # temporary dataframe useful for subsetting data
        tmp_df = pd.DataFrame(
            dict(rt = X.flatten(), sample_name=sample_names.flatten()), 
            index=pep_ids.flatten()
            ).join(self._peptide_df.set_index("external_id"))\
             .sort_values("sample_name")
        sample_bounds = np.arange(tmp_df.shape[0])[~tmp_df.sample_name.duplicated()]
        
        # run through temporary data frame and subset data by run
        for ind in range(len(sample_bounds)):
            front = sample_bounds[ind]
            if ind < len(sample_bounds) - 1:
                back = sample_bounds[ind + 1]
            else:
                back = tmp_df.shape[0]
                
            sample = tmp_df.sample_name.values[front]
            self._sample_weights[sample] = 1.
            self._target_dict[sample] = tmp_df.rt.values[front:back]
            self._select_dict[sample] = tmp_df.internal_id.values[front:back]
            self._model_dict[sample] = self._base_regressor(**self._br_kws)
            
    def _scale_rts(self):
        self._peptide_df.learned_rt -= np.mean(self._peptide_df.learned_rt.values)
        self._peptide_df.learned_rt /= np.std(self._peptide_df.learned_rt.values)
            
    def _run_descent_epoch(self):
        score = 0
        self._scale_rts()
        np.random.shuffle(self._sample_names)
        for sample in self._sample_names:
            # grab run specific target
            target = self._target_dict[sample].reshape(*self._target_shape)
                
            # subset to just peptides relevant to fitting
            select = self._select_dict[sample]
            p = self._peptide_df.learned_rt.values[select].reshape(*self._target_shape)
            weight = self._peptide_df.fit_weight.values[select]
             
            # fit model
            model = self._model_dict[sample]
            model.fit(p, target, weight)
            score += model.score(p, target, weight)
                
            # update latent retention times
            grad_p = model.get_input_grad(p, target).flatten()
            self._peptide_df.learned_rt.values[select] -= weight * self._lr * grad_p
            
        return score / len(self._sample_names)
            
    def _run_minimization_epoch(self):
        score = 0
        self._scale_rts()
        retention_sum = np.zeros_like(self._peptide_df.learned_rt.values)
        weight_sum = np.zeros_like(self._peptide_df.learned_rt.values)
        for sample in self._sample_names:
            # grab run specific target
            target = self._target_dict[sample].reshape(*self._target_shape)
            
            # subset to just peptides relevant to fitting
            select = self._select_dict[sample]
            p = self._peptide_df.learned_rt.values[select].reshape(*self._target_shape)
            weight = self._peptide_df.fit_weight.values[select]
            
            # fit model and accumulate predictions
            model = self._model_dict[sample]
            model.fit(p, target, weight)
            score += model.score(p, target, weight)
            retention_sum[select] += model.invert(p, target).flatten()
            weight_sum[select] += 1.
         
        # update latent retention times
        self._peptide_df.learned_rt = retention_sum/weight_sum
        return score / len(self._sample_names)
            
    def _finalize_fit(self):
        self._scale_rts()
        for sample in self._sample_names:
            # grab run specific target
            target = self._target_dict[sample].reshape(*self._target_shape)
            
            # Subset to just peptides relevant to fitting
            select = self._select_dict[sample]
            p = self._peptide_df.learned_rt.values[select].reshape(*self._target_shape)
            weight = self._peptide_df.fit_weight.values[select]
            
            # Fit model
            model = self._model_dict[sample]
            model.fit(p, target, weight)

    def fit(self, X = None, pep_ids = None, sample_names = None, method=None):
        assert all([X is not None, pep_ids is not None, sample_names is not None])
        self._initialize_params(X, pep_ids, sample_names)
            
        method = method if method is not None else self._method
        assert method in ("min", "descent")
        if method == "min":
            print("Training by global minimization...")
            epoch_fn = self._run_minimization_epoch
        else:
            print("Training by gradient descent...")
            epoch_fn = self._run_descent_epoch
        
        for epoch in range(self._max_epochs):
            epoch_score = epoch_fn()
            self._history.append(epoch_score)
            print(
                "Score on epoch {}:".format(len(self._history)), 
                self._history[-1]
            )
            
        self._finalize_fit()
    
    def get_retention(self):
        return self._peptide_df.drop(["internal_id", "fit_weight"], axis=1)\
                               .rename({"external_id" : "pep_id"}, axis=1)
    
    def _predict_single_sample(self, internal_ids, sample_name):
        try:
            model = self._model_dict[sample_name]
            pred = model.predict(
                self._peptide_df.learned_rt.values[internal_ids].reshape(*self._target_shape)
            )
            
            return pred.flatten()
            
        except KeyError:
            return np.nan()

    def predict(self, pep_ids, sample_names):
        pred_df = self._peptide_df.join(
            pd.DataFrame({"sample_name" : sample_names.flatten()}, index=pep_ids.flatten()),
            on="external_id"
            )
        
        pred_df = pred_df.groupby("sample_name").apply(
            lambda df: pd.DataFrame(
                {
                     "predicted_rt": self._predict_single_sample(
                          df.internal_id,
                          df.sample_name.iloc[0]
                )},
                index = pd.Index(df.external_id, name="pep_id")
            ),
        ).reset_index()
        
        return pred_df
    
    def apply(self, X, pep_ids, sample_names):
        print("Applying learned models to new data...")
        result = pd.DataFrame({"pep_id" : np.unique(pep_ids)})
        result["occurence"] = np.zeros(result.shape[0])
        result["learned_rt"] = np.zeros(result.shape[0])
        result = result.set_index("pep_id")
        
        tmp_df = pd.DataFrame(dict(
            pep_id=pep_ids.flatten(),
            sample_name=sample_names.flatten(),
            rt = X.flatten(), 
            )
        ).sort_values(["sample_name", "pep_id"])
        sample_bounds = np.arange(tmp_df.shape[0])[~tmp_df.sample_name.duplicated()]
        
        for ind in range(len(sample_bounds)):
            front = sample_bounds[ind]
            if ind < len(sample_bounds) - 1:
                back = sample_bounds[ind + 1]
            else:
                back = tmp_df.shape[0]
                
            sample = tmp_df.sample_name.values[front]
            model = self._model_dict[sample]
            rts = tmp_df.rt.iloc[front:back].values
            pred = model.invert(np.zeros_like(rts), rts)
            
            result.occurence.loc[tmp_df.pep_id.iloc[front:back]] += 1
            result.learned_rt.loc[tmp_df.pep_id.iloc[front:back]] += pred
        
        result.learned_rt /= result.occurence
        print("Average spearman correlation: {:.3f}".format(
            tmp_df.join(result, on = "pep_id")\
                  .groupby("sample_name")[["rt", "learned_rt"]]\
                  .corr("spearman")\
                  .iloc[0::2,-1]\
                  .mean()
        ))
                
        return result.reset_index()

class RTCalculator:
    def __init__(self, model, seed, **kwargs):
        self.seed = seed
        assert model in ("linear", "isotonic")
        if model == "linear":
            base_regressor = AugmentedLinearRegression
            br_kws = {}
            self.target_shape = (-1, 1)
        else:
            base_regressor = AugmentedIsotonicRegression
            br_kws = {"out_of_bounds" : "clip"}
            self.target_shape = (-1,)

        self.model = RTAligner(
            base_regressor, 
            br_kws={},
            **kwargs
        )

    def _load_data(self, path):
        print("Loading Peptides...")
        peptides = pd.read_csv(os.path.join(path, "peptides.csv"))
        peptides = peptides.drop(["hits", "learned_rt", "test_error"],
                                 axis=1, errors="ignore")

        print("Loading PSMs...")
        psms = pd.read_csv(os.path.join(path, "psms.csv"),
                           dtype={"precursor_charge": str})
        psms = psms.loc[:, ["id", "pep_id", 
                            "sample_name",
                            "score", "label",
                            "qvalue", "scan_rt"]]

        return peptides, psms

    def _filter_data(self, peptides, psms):
        score_cutoff = peptides[peptides.qvalue < 0.01].score.min()
        peptides = peptides[
            np.logical_and(
                peptides.label == "target",
                peptides.score > score_cutoff
            )]
        psms = psms[
            np.logical_and(
                psms.label == "target",
                psms.score > score_cutoff
            )]
        psms = peptides.id.rename("pep_id")\
                       .to_frame()\
                       .join(psms.set_index("pep_id"),
                             on = "pep_id",
                             how = "left")
        psms = psms.sort_values("qvalue")\
                   .groupby(["pep_id", "sample_name"])\
                   .scan_rt\
                   .first()\
                   .reset_index()

        file_count = psms.groupby("sample_name")\
                         .pep_id\
                         .count()\
                         .rename("total_peptides")

        psms = psms.join(file_count, on="sample_name", how="right")
        psms = psms[psms.total_peptides > 100]

        return psms

    def _partition_psms(self, psms):
        peptide_count = psms.groupby("pep_id")\
                            .scan_rt\
                            .count()\
                            .rename("hits")
        psms = psms.join(peptide_count, on="pep_id")
        high_count_psms = psms[psms.hits >= 5]
        low_count_psms = psms[psms.hits < 5]

        return high_count_psms, low_count_psms

    def _split_and_scale(self, detections):
        train, test = train_test_split(
            detections,
            test_size=.2,
            stratify = detections.pep_id
        )

        train_mean = train.groupby("sample_name")\
                          .scan_rt\
                          .mean()\
                          .rename("train_mean")
        train_std = train.groupby("sample_name")\
                         .scan_rt\
                         .std()\
                         .rename("train_std")

        train = train.join(train_mean, on="sample_name")\
                     .join(train_std, on="sample_name")
        train["scaled_rt"] = (train.scan_rt - train.train_mean)/train.train_std

        test = test.join(train_mean, on="sample_name")\
                   .join(train_std, on="sample_name")
        test["scaled_rt"] = (test.scan_rt - test.train_mean)/test.train_std

        return train, test

    def _evaluate(self, test):
        print("Evaluating on test data...")
        test_corr = test.groupby("sample_name")[["scaled_rt", "predicted_rt"]]\
                        .corr().iloc[0::2,-1].mean()
        quantiles = np.linspace(.1, .9, 9)
        errors = test.groupby("sample_name")\
                     .abs_scaled_error\
                     .median()\
                     .quantile(quantiles)
        print("Pearson's correlation:\t{:.4f}".format(test_corr))
        print("Quantiles of scaled absolute error:")
        for q, v in zip(quantiles, errors):
            print("\t{:.2f}:\t{:.3f}".format(q, v))

    def process_path(self, path):
        np.random.seed(self.seed)
        peptides, psms = self._load_data(path)
        psms = self._filter_data(peptides, psms)
        high_count_detections, low_count_detections = self._partition_psms(psms)
        train, test = self._split_and_scale(high_count_detections)

        self.model.fit(
            X = train.scaled_rt.values.reshape(self.target_shape), 
            pep_ids = train.pep_id.values.reshape(self.target_shape), 
            sample_names = train.sample_name.values.reshape(self.target_shape)
        )
        test = test.join(
            self.model.predict(
                pep_ids = test.pep_id.values.reshape(self.target_shape),
                sample_names = test.sample_name.values.reshape(self.target_shape)
            ).set_index(["sample_name", "pep_id"]),
            on=["sample_name", "pep_id"]
        )
        test["abs_scaled_error"] = np.abs(test.scaled_rt - test.predicted_rt)
        self._evaluate(test)

        low_count_detections = low_count_detections.join(
            train.groupby("sample_name").first()[["train_mean", "train_std"]],
            on = "sample_name", how = "inner"
        )
        low_count_detections["scaled_rt"] = (
            low_count_detections.scan_rt - low_count_detections.train_mean
            )/low_count_detections.train_std
        retention_times = pd.concat([
            self.model.get_retention(),
            self.model.apply(
                X = low_count_detections.scaled_rt.values.reshape(self.target_shape), 
                pep_ids = low_count_detections.pep_id.values.reshape(self.target_shape), 
                sample_names = low_count_detections.sample_name.values.reshape(self.target_shape)
            )
        ])

        lower_bound = retention_times.learned_rt.quantile(0.001)
        upper_bound = retention_times.learned_rt.quantile(.999)
        retention_times.learned_rt -= lower_bound
        retention_times.learned_rt *= 100./(upper_bound - lower_bound)
        test.abs_scaled_error *= 100./(upper_bound - lower_bound)

        peptides = peptides.join(retention_times.set_index("pep_id"), 
                                 on = "id")\
                           .join(test.groupby("pep_id")\
                                     .abs_scaled_error\
                                     .median(),\
                                 on = "id")\
                           .rename({"occurence" : "hits",
                                    "abs_scaled_error" : "test_error"},
                                   axis=1)
        peptides.to_csv(os.path.join(path, "peptides.csv"))
