import os
import numpy as np
import pandas as pd

class DatabaseBuild:
    def __init__(self, build_path, fdr_cutoff = 1e-2):
        self.build_path = build_path
        self.fdr_cutoff = fdr_cutoff
        self.psms = self._load_psms()
        self.peptides, self.peptide_modifications = self._load_peptides()
        self.proteins, self.peptide_protein = self._load_proteins()
        self.sites = self._load_sites()

    def _load_psms(self):
        psms = pd.read_csv(os.path.join(self.build_path,
                                        "psms.csv"), 
                           na_values="None")
        psms = psms[np.logical_and(psms.qvalue <= self.fdr_cutoff,
                                   psms.label == "target")]

        return psms

    def _load_peptides(self):
        peptides = pd.read_csv(os.path.join(self.build_path,
                                            "peptides.csv"), 
                               na_values="None")
        peptides = peptides[np.logical_and(peptides.qvalue <= self.fdr_cutoff,
                                           peptides.label == "target")]

        peptide_modifications =  pd.read_csv(os.path.join(self.build_path, 
                                                          "peptide_modifications.csv"),
                                             na_values="None")
        peptide_modifications = peptides[["id"]].rename({"id": "pep_id"}, axis=1)\
                                                .join(peptide_modifications.set_index("pep_id"),
                                                      how="inner",
                                                      on="pep_id")
        return peptides, peptide_modifications

    def _load_proteins(self):
        peptide_protein = pd.read_csv(os.path.join(self.build_path,
                                                   "peptide_protein.csv"),
                                      na_values="None")
        peptide_protein = self.peptides[["id"]].rename({"id": "pep_id"}, axis=1)\
                                               .join(peptide_protein.set_index("pep_id"),
                                                     on="pep_id")

        proteins = pd.read_csv(os.path.join(self.build_path,
                                            "proteins.csv"),
                               na_values="None")
        proteins = peptide_protein.prot_id\
                                  .drop_duplicates()\
                                  .rename("id")\
                                  .to_frame()\
                                  .join(proteins.set_index("id"),
                                        on="id")

        return proteins, peptide_protein

    def _load_sites(self):
        sites = pd.read_csv(os.path.join(self.build_path,
                                         "sites.csv"),
                            na_values="None")
        sites = sites[sites.qvalue <= self.fdr_cutoff]

        sites = self.proteins.id\
                             .rename("prot_id")\
                             .to_frame()\
                             .join(sites.set_index("prot_id"),
                                   how="inner",
                                   on="prot_id")

        return sites

