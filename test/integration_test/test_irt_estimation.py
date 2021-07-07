import unittest
import numpy as np
import scipy as sc
import pandas as pd
import matplotlib.pyplot as plt
from integration_engine.rt import *


def random_linear(x, sd):
    coefs = [np.random.randn(),
             np.abs(np.random.randn())]
    x = coefs[0] + coefs[1] * x
    x = x + sd*np.random.randn(*x.shape)
    return x

def random_nonlinear(x, sd):
    coefs = [np.random.randn(),
             np.abs(np.random.randn() + 2),
             0.5*np.abs(np.random.randn())]
    x = coefs[0] + x + coefs[2] * (x**3)
    x = x + sd*np.random.randn(*x.shape)
    return x


class TestRegressors(unittest.TestCase):
    def test_linear_regression(self):
        # Build data
        np.random.seed(98115)
        test_x = np.random.randn(1000, 1)
        test_y = random_linear(test_x, sd=.1)

        # Test model
        model = AugmentedLinearRegression()
        model.fit(test_x, test_y)
        inverted_y = model.invert(test_x, test_y)

        var = sum((test_x - inverted_y)**2)/(test_x.shape[0] - 1)
        self.assertTrue(np.isclose(var, .1, rtol=0., atol=.01))
        
        grad_x = model.get_input_grad(test_x, test_y)
        self.assertTrue(np.isclose(np.mean(grad_x), 0.))

    def test_isotonic_regression(self):
        # Build data
        np.random.seed(98115)
        test_x = np.random.randn(1000)
        test_y = random_nonlinear(test_x, sd=.1)

        # Test model
        model = AugmentedIsotonicRegression()
        model.fit(test_x, test_y)
        inverted_y = model.invert(test_x, test_y)

        var = sum((test_x - inverted_y)**2)/(test_x.shape[0] - 1)
        self.assertTrue(var < .01)

        grad_x = model.get_input_grad(test_x, test_y)
        self.assertTrue(np.isclose(np.mean(grad_x), 0.))


class TestAlignmentClass(unittest.TestCase):
    def test_model_initialization(self):
        # Build data
        np.random.seed(98115)
        test_x = np.random.randn(1000)
        n = test_x.shape[0]
        data = pd.concat([
            pd.DataFrame({"sample_name" : ["sample_1"]*n,
                          "pep_id" : np.arange(n),
                          "rt" : test_x}),
            pd.DataFrame({"sample_name" : ["sample_2"]*n,
                          "pep_id" : np.arange(n),
                          "rt" : test_x})
            ])

        # Test initialization
        print()
        model = RTAligner(AugmentedLinearRegression)
        model._initialize_params(X=data.rt.values,
                                 pep_ids=data.pep_id.values,
                                 sample_names=data.sample_name.values)
        self.assertTrue(
            np.all(model._peptide_df.internal_id == model._peptide_df.external_id)
        )
        self.assertTrue(
            np.isclose(
                np.mean(test_x - model._peptide_df.learned_rt.values.flatten()),
                0.
            )
        )
                                 

    def test_linear_alignment(self):
        # Build data
        np.random.seed(98115)
        n = 1000
        overlap = 250
        test_x = np.random.randn(n)

        data_list = [pd.DataFrame({"sample_name" : ["sample_0"]*n,
                                   "pep_id" : np.arange(n),
                                   "rt" : random_linear(test_x, .1)})]
        for sample_ind in range(1, 11):
            select = np.random.choice(np.arange(n),
                                      overlap,
                                      replace=False)
            data_list.append(
                 pd.DataFrame({"sample_name" : [f"sample_{sample_ind}"]*overlap,
                               "pep_id" : np.arange(n)[select],
                               "rt" : random_linear(test_x[select], .1)})
            )

        data = pd.concat(data_list)

        # Get initialized corr
        print()
        model = RTAligner(AugmentedLinearRegression)
        model._initialize_params(X=data.rt.values,
                                 pep_ids=data.pep_id.values,
                                 sample_names=data.sample_name.values)
        initial_corr = sc.stats.spearmanr(
                           model.get_retention().learned_rt.values.flatten(),
                           test_x
                       )[0]

        # Train model by descent
        print()
        model = RTAligner(AugmentedLinearRegression,
                          method="descent",
                          lr=1e-1,
                          max_epochs=10)
        model.fit(X=data.rt.values.reshape(-1, 1),
                  pep_ids=data.pep_id.values.reshape(-1, 1),
                  sample_names=data.sample_name.values.reshape(-1, 1))
        corr = sc.stats.spearmanr(
                   model.get_retention().learned_rt.values.flatten(),
                   test_x
               )[0]
        self.assertTrue(corr > initial_corr)

        # Train model by minimizing
        print()
        model = RTAligner(AugmentedLinearRegression,
                          method="min",
                          max_epochs=5)
        model.fit(X=data.rt.values.reshape(-1, 1),
                  pep_ids=data.pep_id.values.reshape(-1, 1),
                  sample_names=data.sample_name.values.reshape(-1, 1))
        corr = sc.stats.spearmanr(
                   model.get_retention().learned_rt.values.flatten(),
                   test_x
               )[0]
        self.assertTrue(corr > initial_corr)


    def test_nonlinear_alignment(self):
        # Build data
        np.random.seed(98115)
        n = 1000
        overlap = 250
        test_x = np.random.randn(n)

        data_list = [pd.DataFrame({"sample_name" : ["sample_0"]*n,
                                   "pep_id" : np.arange(n),
                                   "rt" : random_linear(test_x, .1)})]
        for sample_ind in range(1, 11):
            select = np.random.choice(np.arange(n), 
                                      overlap, 
                                      replace=False)
            random_func = np.random.choice([random_linear, random_nonlinear])
            data_list.append(
                 pd.DataFrame({"sample_name" : [f"sample_{sample_ind}"]*overlap,
                               "pep_id" : np.arange(n)[select],
                               "rt" : random_func(test_x[select], .1)})
            )

        data = pd.concat(data_list)

        # Get initialized corr
        print()
        model = RTAligner(AugmentedLinearRegression)
        model._initialize_params(X=data.rt.values,
                                 pep_ids=data.pep_id.values,
                                 sample_names=data.sample_name.values)
        initial_corr = sc.stats.spearmanr(
                           model.get_retention().learned_rt.values.flatten(),
                           test_x
                       )[0]

        # Train model by descent
        print()
        model = RTAligner(AugmentedIsotonicRegression,
                          method="descent",
                          lr=1e-1,
                          max_epochs=10)
        model.fit(X=data.rt.values,
                  pep_ids=data.pep_id.values,
                  sample_names=data.sample_name.values)
        corr = sc.stats.spearmanr(
                   model.get_retention().learned_rt.values.flatten(),
                   test_x
               )[0]
        self.assertTrue(corr > initial_corr)


class TestRTCalculator(unittest.TestCase):
    def test_detection_split(self):
        # Build data
        np.random.seed(98115)
        test_x = np.random.randn(1000)

        data_list = []
        for nsamples, nvar in zip([4, 1], [1000, 250]): 
            for sample_ind in range(nsamples):
                data_list.append(
                    pd.DataFrame({"sample_name" : [f"sample_{sample_ind}"]*nvar,
                                  "pep_id" : np.arange(nvar),
                                  "scan_rt" : random_linear(test_x[:nvar], .1)})
                    )

        data = pd.concat(data_list)

        calculator = RTCalculator("linear", 98115)
        high_count, low_count = calculator._partition_psms(data)
        self.assertEqual(high_count.shape[0], 1250)
        self.assertEqual(low_count.shape[0], 3000)
