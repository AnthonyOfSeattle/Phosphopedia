import os
import numpy as np
import pandas as pd

class FDRCalculator:
    def _calculate_fdr(self, score_array, label_array):
        # Arguments should be sorted such that best scores come first
        # This function is very similar to that provided by Percolator
        qvalues = []
    
        score_array = np.atleast_1d(score_array)
        label_array = np.atleast_1d(label_array)
        ntargets = np.cumsum(label_array == "target")
        ndecoys = np.arange(1, ntargets.shape[0] + 1) - ntargets
        qvalues = (ndecoys + 1)/ntargets
        for ind in np.arange(qvalues.shape[0])[::-1]:
             if ind - 1 >= 0 and score_array[ind] == score_array[ind-1]:
                 qvalues[ind - 1] = qvalues[ind]
    
             if ind + 1 < qvalues.shape[0]:
                 qvalues[ind] = min(qvalues[ind], qvalues[ind + 1])
    
        return qvalues

    def _process_psms(self, path):
        print("Loading PSMs...")
        psms = pd.read_csv(os.path.join(path, "psms.csv"),
                           dtype={"precursor_charge": str})
        
        print("Sorting PSMs by score...")
        psms.sort_values("score", inplace=True, ascending=False)

        print("Calculating PSM FDR...")
        psms["qvalue"] = self._calculate_fdr(
            psms.score.values, psms.label.values
        )

        print("Sorting PSMs by ID...")
        psms.sort_values("id", inplace=True, ascending=True)
        
        print("Dumping PSMs")
        psms.to_csv(os.path.join(path, "psms.csv"), index=False)

    def _process_peptides(self, path):
        print("Loading Peptides...")
        peptides = pd.read_csv(os.path.join(path, "peptides.csv"))

        print("Sorting Peptides by score...")
        peptides.sort_values("score", inplace=True, ascending=False)

        print("Calculating Peptide FDR...")
        peptides["qvalue"] = self._calculate_fdr(
            peptides.score.values, peptides.label.values
        )

        print("Sorting PSMs by ID...")
        peptides.sort_values("id", inplace=True, ascending=True)

        print("Dumping Peptides...")
        peptides.to_csv(os.path.join(path, "peptides.csv"), index=False)

    def _process_sites(self, path):
        print("Loading and annotating Phosphosites...")
        prot_data = pd.read_csv(os.path.join(path, "proteins.csv"), index_col="id")
        sites = pd.read_csv(os.path.join(path, "sites.csv"))
        sites = sites.join(prot_data, on="prot_id")

        print("Sorting Phosphosites on score...")
        sites.sort_values("score", inplace=True, ascending=False)

        print("Calculating Phosphosite FDR...")
        sites["qvalue"] = self._calculate_fdr(
            sites.score.values, sites.label.values
        )

        print("Dumping Phosphosites...")
        sites.drop(["label", "accession"], axis=1).to_csv(os.path.join(path, "sites.csv"), index=False)

    def process_path(self, path):
        print("Calculating the FDR for: " + path)
        self._process_psms(path)
        self._process_peptides(path)
        self._process_sites(path)

