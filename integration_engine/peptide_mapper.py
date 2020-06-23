import pickle
import numpy as np
from itertools import groupby
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

class PeptideMapper:

    @staticmethod
    def establish_connection(db_path):
        global SESSION
        SESSION = sessionmaker(bind=create_engine(db_path))()

    @staticmethod
    def get_localized_modifications(records):
        localized = {}

        # Gather all modifications from all peptides
        all_modifications = []
        [all_modifications.extend(psm["modifications"].items()) for psm in records]

        # Set positions in localized map to be the highest scores
        for pos, mod in all_modifications:
            localized[pos] = max(localized.get(pos, 0.),
                                 mod["score"])

            #if mod["score"] < 13:
            #    [localized.setdefault(int(alt_pos), 0.) for alt_pos in mod.alternative_positions.split(",")]

        return localized

    @staticmethod
    def build_peptide(psm, localized_ptms):
        tokenized_seq = list(psm["base_sequence"])
        peptide_mods = {pos: dict(info) for pos, info in psm["modifications"].items()}
        for pos in list(peptide_mods.keys()):
            mod = peptide_mods.pop(pos)

            if mod["score"] < 13:
                for alt_pos in mod["alt_pos"]:
                    if localized_ptms.get(alt_pos, 0.) > mod["score"] and not peptide_mods.get(alt_pos, False):
                        pos = alt_pos

            peptide_mods[pos] = mod
            if pos == 0:
                tokenized_seq[0] = "n" + "[{:.0f}]".format(mod["mass"]) + tokenized_seq[0]
            else:
                tokenized_seq[pos - 1] = tokenized_seq[pos - 1] + "[{:.0f}]".format(mod["mass"])

        peptide = dict(psm_id = psm["id"],
                       label = psm["label"],
                       base_sequence = psm["base_sequence"],
                       sequence = "".join(tokenized_seq),
                       score = psm["score"],
                       modifications = pickle.dumps(peptide_mods))

        return peptide

    @staticmethod
    def collapse_peptides(peptide_group):
        psm_ids = []
        best_peptide = next(peptide_group[1])
        psm_ids.append(best_peptide["psm_id"])
        for pep in peptide_group[1]:
            psm_ids.append(pep["psm_id"])
            if pep["score"] > best_peptide["score"]:
                best_peptide = pep

        best_peptide["psm_ids"] = psm_ids
        return best_peptide

    @staticmethod
    def merge_psms(records):
        localized_ptms = PeptideMapper.get_localized_modifications(records)
        unmerged_peptides = [PeptideMapper.build_peptide(psm, localized_ptms) for psm in records]
        unmerged_peptides.sort(key=lambda p: p["sequence"])
        peptide_groups = groupby(unmerged_peptides, key=lambda p: p["sequence"])
        return list(map(PeptideMapper.collapse_peptides, peptide_groups))

    @staticmethod
    def run(records):
        for r in records:
            r["modifications"] = pickle.loads(r["modifications"])

        if len(np.unique([r["base_sequence"] for r in records])) > 1:
            print("Hash collision!")

        merged_psms = PeptideMapper.merge_psms(records)

        return merged_psms
