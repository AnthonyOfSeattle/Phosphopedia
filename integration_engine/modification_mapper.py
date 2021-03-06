import pickle
import numpy as np

class ModificationMapper:

    @staticmethod
    def get_mod_pos(pep_seq, mods, prot_seq, is_decoy=False):
        prot_mod_pos = []
        start = prot_seq.find(pep_seq)
        for pep_pos in mods:
            if pep_pos == 0:
                pos = start

            elif is_decoy:
                pos = start + (len(pep_seq) - pep_pos - 1)

            else:
                pos = start + (pep_pos - 1)

            prot_mod_pos.append((pos, prot_seq[pos]))

        return prot_mod_pos

    @staticmethod
    def run(prot):
        mod_scores = {}
        for peptide in prot["peptide_list"]:
            seq, score, mods = peptide
            mods = pickle.loads(mods)
            mods = {pos : info for pos, info in mods.items() if info["residue"] in "STY"}
            if prot["label"] == "decoy":
                seq = seq[:-1][::-1] + seq[-1]
                prot_mod_pos = ModificationMapper.get_mod_pos(seq, mods, prot["seq"], True)
            else:
                prot_mod_pos = ModificationMapper.get_mod_pos(seq, mods, prot["seq"])

            for site in prot_mod_pos:
                if mod_scores.get(site, -np.inf) < score:
                    mod_scores[site] = score

        return prot["id"], mod_scores
