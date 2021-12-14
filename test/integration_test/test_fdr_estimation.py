import unittest
import tempfile
import os
import numpy as np
import pandas as pd
from integration_engine.fdr import *


# Define truth for all tests
SCORE_ARRAY = np.arange(20)[::-1]
LABEL_ARRAY = np.array(["target", "target", "target", "target", "target",
                        "target", "target", "decoy", "target", "target",
                        "target", "decoy", "target", "target", "decoy",
                        "target", "target", "target", "decoy", "decoy"])
TRUE_QVALUES = np.array([0.1428571, 0.1428571, 0.1428571, 0.1428571, 0.1428571,
                         0.1428571, 0.1428571, 0.2000000, 0.2000000, 0.2000000,
                         0.2000000, 0.2500000, 0.2500000, 0.2500000, 0.2666667,
                         0.2666667, 0.2666667, 0.2666667, 0.3333333, 0.4000000])


def get_fake_psm_data():
    psm_data = pd.DataFrame({"score" : SCORE_ARRAY,
                             "label" : LABEL_ARRAY})
    psm_data["id"] = np.random.choice(np.arange(SCORE_ARRAY.shape[0]),
                                      SCORE_ARRAY.shape[0],
                                      replace=False)
    for colname in ["sample_name", "scan_number", "scan_rt",
                    "precursor_mz", "precursor_charge", "pep_id"]:
        psm_data[colname] = np.nan

    psm_data = psm_data.sort_values("id")\
                       .loc[:, ["id", "sample_name", "scan_number",
                                "scan_rt", "precursor_mz", "precursor_charge",
                                "score", "pep_id", "label"]]
    return psm_data


def get_fake_peptide_data():
    peptide_data = pd.DataFrame({"score" : SCORE_ARRAY,
                                 "label" : LABEL_ARRAY})
    peptide_data["id"] = np.random.choice(np.arange(SCORE_ARRAY.shape[0]),
                                          SCORE_ARRAY.shape[0],
                                          replace=False)
    peptide_data["sequence"] = np.nan

    peptide_data = peptide_data.sort_values("id")\
                               .loc[:, ["id", "sequence", "score", "label"]]
    return peptide_data


def get_fake_site_data():
    protein_data = pd.DataFrame({"id" : [1, 2],
                                 "accession" : ["sp|target_protein",
                                                "sp|decoy_protein"],
                                 "label" : ["target", "decoy"]})


    site_data = pd.DataFrame({"score" : SCORE_ARRAY,
                              "prot_id" : [{"target" : 1, "decoy" : 2}[l]
                                           for l in LABEL_ARRAY]})
    site_data["position"] = np.random.choice(np.arange(SCORE_ARRAY.shape[0]),
                                             SCORE_ARRAY.shape[0],
                                             replace=False)
    site_data["residue"] = np.nan

    site_data = site_data.sort_values(["prot_id", "position"])\
                         .loc[:, ["prot_id", "position", "residue", "score"]]
    return protein_data, site_data


class TestFDRCalculator(unittest.TestCase):
    def test_qvalue_calculation(self):
        calc = FDRCalculator()
        qvalues = calc._calculate_fdr(SCORE_ARRAY, LABEL_ARRAY)
        self.assertTrue(np.all(np.isclose(TRUE_QVALUES, qvalues))) 

    def test_folder_processing(self):

        calc = FDRCalculator()
        with tempfile.TemporaryDirectory() as temp_path:
            # Test psm calculations
            psm_file = os.path.join(temp_path, "psms.csv")
            psms = get_fake_psm_data()
            psms.to_csv(psm_file, index=False)

            calc._process_psms(temp_path)

            psms = pd.read_csv(psm_file)
            psms = psms.sort_values("score", ascending=False)
            self.assertTrue(np.all(np.isclose(TRUE_QVALUES, psms.qvalue)))
            
            # Test peptide calculations
            peptide_file = os.path.join(temp_path, "peptides.csv")
            peptides = get_fake_peptide_data()
            peptides.to_csv(peptide_file, index=False)

            calc._process_peptides(temp_path)

            peptides = pd.read_csv(peptide_file)
            peptides = peptides.sort_values("score", ascending=False)
            self.assertTrue(np.all(np.isclose(TRUE_QVALUES, peptides.qvalue)))

            # Test peptide calculations
            protein_file = os.path.join(temp_path, "proteins.csv")
            site_file = os.path.join(temp_path, "sites.csv")
            proteins, sites = get_fake_site_data()
            proteins.to_csv(protein_file, index=False)
            sites.to_csv(site_file, index=False)

            calc._process_sites(temp_path)

            sites = pd.read_csv(site_file)
            sites = sites.sort_values("score", ascending=False)
            self.assertTrue(np.all(np.isclose(TRUE_QVALUES, sites.qvalue)))

