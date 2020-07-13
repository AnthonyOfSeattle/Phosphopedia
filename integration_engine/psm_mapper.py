import re
import os
import pickle
import numpy as np
import pandas as pd

class PSMMapper:

    @staticmethod
    def initialize(protein_groups, group_number):
        global PROTEIN_GROUPS
        PROTEIN_GROUPS = protein_groups

        global GROUP_NUMBER
        GROUP_NUMBER = group_number

    @staticmethod
    def gather_data(scan_info_path, psm_path, localization_path):
        scan_info = pd.read_csv(scan_info_path, sep="\t").set_index("scan")
        localizations = pd.read_csv(localization_path, sep="\t").set_index("Scan")
        psm_scores = pd.read_csv(psm_path, sep="\t",
                                 usecols=["scan",
                                          "percolator score",
                                          "percolator q-value",
                                          "percolator PEP",
                                          "protein id"]).set_index("scan", drop=False)
        if PROTEIN_GROUPS and GROUP_NUMBER >= 0:
            select = [PROTEIN_GROUPS.get(ref_set.split(",")[0].lstrip("decoy_"), -1) == GROUP_NUMBER
                      for ref_set in psm_scores["protein id"]]
            psm_scores = psm_scores[select]
        localized_scores = psm_scores.join(localizations, how="left")\
                                     .join(scan_info, how="left")
        localized_scores = localized_scores[~localized_scores.PepScore.isna()]
        return localized_scores

    @staticmethod
    def yield_mod(seq):
        for ind, mod in enumerate(re.finditer("[A-Zn](\[[^A-Z]+\])?", seq), 1):
            mod_mass = re.search("(?<=\[)([^A-Z]*)(?=\])", mod.group())
            if mod_mass is not None:
                # Subtract 1 if the first character is an n-terminal designation
                yield mod.group()[0], ind - int(seq[0] == "n"), float(mod_mass.group())

    @staticmethod
    def create_entries(row, label, basename):
        # Puts all modifications in a dictionary to be jsonified
        mods = {}
        phospho_scores = iter(zip(row["Ascores"].split(";"), str(row["AltSites"]).split(";")))
        for mod_res, mod_ind, mod_mass in PSMMapper.yield_mod(row["LocalizedSequence"]):
            score = np.inf
            alt_pos = []
            if np.isclose(79.966331, mod_mass, rtol=0, atol=1):
                score, alt_pos = next(phospho_scores)
                score = float(score)
                if alt_pos and alt_pos != "nan":
                  alt_pos = [int(i) for i in alt_pos.split(",")]

            mods[mod_ind] = {"score": score, "residue" : mod_res, "mass": mod_mass, "alt_pos": alt_pos}

        psm = dict(sample_name = basename,
                   scan_number = row["scan"],
                   scan_rt = row["rt"],
                   precursor_mz = row["precursor"],
                   precursor_charge = row["charge"],
                   label = label,
                   base_sequence =  re.sub("[^A-Z]+", "", row["LocalizedSequence"]),
                   score = row["percolator score"],
                   qvalue = row["percolator q-value"],
                   pep = row["percolator PEP"],
                   modifications = pickle.dumps(mods))
        psm["proteins"] = ["|".join(seq.split("|")[:-1]) for seq in row["protein id"].split(",")]
        psm["hash_value"] = hash(psm["base_sequence"])

        return psm

    @staticmethod
    def run(file_tuple):
        localized_scores = PSMMapper.gather_data(*file_tuple)
        psm_entries = localized_scores.apply(PSMMapper.create_entries, axis=1,
                                             label=os.path.split(file_tuple[1])[1].split(".")[2],
                                             basename=os.path.split(file_tuple[1])[1].split(".")[0]).tolist()
        return psm_entries
