import re
import numpy as np
from numpy import isclose
from pyteomics.pepxml import PepXML
from pyteomics import mass
std_aa_mass = mass.std_aa_mass

class AbstractExtractor:
    def __init__(self, score_string=None):
        # Storage 
        self.match_list = None
        self.match = None
        self.score_string = score_string

    def _get_score(self):
        try:
            return self.match[self.score_string]
        except KeyError:
            return None
    
    def _initialize_results(self):
        nmatches = len(self.match_list)
        fields = ["scans", "scores", "charge_states",
                  "peptides", "mod_positions", "mod_masses"]
        self.results = {}
        for f in fields:
            self.results[f] = [None] * nmatches

    def extract(self, entry):
        self.entry = entry
        self.match_list = self._get_matches()

        self._initialize_results()
        for ind, self.match in enumerate(self.match_list):
            self.results["scans"][ind] = self._get_scan()
            self.results["scores"][ind] = self._get_score()
            self.results["charge_states"][ind] = self._get_charge()
            self.results["peptides"][ind] = self._get_peptide()
            self.results["mod_positions"][ind], self.results["mod_masses"][ind] = self._get_mod_info()

        return self.results

class PepXMLExtractor(AbstractExtractor):
    def _get_matches(self):
        try:
            return self.entry["search_hit"]
        except KeyError:
            return []

    def _get_scan(self):
        try:
            return int(self.entry["start_scan"])
        except KeyError:
            return -1

    def _get_charge(self):
        try:
            return int(self.entry["assumed_charge"])
        except KeyError:
            return 0

    def _get_score(self):
        try:
            return float(self.match["search_score"][self.score_string])
        except KeyError:
            return None

    def _get_peptide(self):
        try:
            return self.match["peptide"]
        except KeyError:
            return ""

    def _get_mod_info(self):
        try:
            mod_list = self.match["modifications"]
            pos = np.zeros(len(mod_list), dtype=np.int32)
            mass = np.zeros(len(mod_list), dtype=np.float32)
            for ind, mod in enumerate(mod_list):
                pos[ind] = int(mod["position"])
                mass[ind] = float(mod["mass"])
            return pos, mass

        except KeyError:
            return np.zeros(0, dtype=np.int32), np.zeros(0, dtype=np.float32)


class ModificationParser:
    def __init__(self, spec_file_name, 
                 match_file_name, 
                 match_file_format,
                 mod_masses={},
                 masses_absolute=False,
                 constant_mods=['C'], 
                 score_string=None,
                 score_threshold=None,
                 score_higher_better=False,
                 score_func=None):
        self.spec_file_name = spec_file_name

        self.mod_data = []
        self.constant_mods = constant_mods
        if not masses_absolute:
            self.mod_masses = self._make_absolute_masses(mod_masses)
        else:
            self.mod_masses = mod_masses
        
        assert all([self.mod_masses.get(const, 0.) > 0. for const in self.constant_mods]), \
               "All constant mods should have non-zero specified masses"

        self.score_threshold = score_threshold
        self.score_higher_better = score_higher_better
        self.score_func = score_func
        if match_file_format == "pepXML":
            self.reader = PepXML(match_file_name)
            self.extractor = PepXMLExtractor(score_string)
        else:
            raise ValueError("{} not supported at this time.".format(file_format))

    def _make_absolute_masses(self, mod_masses):
        absolute_masses = {}
        for res, mass in mod_masses.items():
            absolute_masses[res] = std_aa_mass.get(res, 0.) + mass
        return absolute_masses

    def extract(self):
        self.mod_data = [e for e in self.reader.map(self.extractor.extract) if len(e['peptides']) > 0]
        self.mod_data = sorted(self.mod_data, key=lambda e: e["scans"][0])

    def _correct_mass(self, res, pos, mass):
        mod_tol = 1.5
        std_mass = std_aa_mass.get(res, np.inf)
        mod_mass = self.mod_masses.get(res, np.inf)
        nmod_mass = self.mod_masses.get('n', np.inf)
        if pos == -1 and isclose(mass, nmod_mass,
                                 rtol=0., atol=mod_tol):
            return (res,), (-100,), (nmod_mass,)

        elif pos == 1 and isclose(mass, std_mass + nmod_mass,
                                  rtol=0., atol=mod_tol):
            return (res,), (-100,), (nmod_mass,)

        elif pos == 1 and isclose(mass, mod_mass + nmod_mass,
                                  rtol=0., atol=mod_tol):
            return ('n', res), (-100, pos), (nmod_mass, mod_mass)

        elif isclose(mass, mod_mass,
                     rtol=0., atol=mod_tol):
            return (res,), (pos,), (mod_mass,)

        else:
            raise ValueError("Unrecognized mod on {} at position {} with mass: {}".format(res, pos, mass))

    def _format_mods(self, peptide, positions, masses):
        handles_constant_mods = False
        mod_strings = []
        for pos, mass in zip(positions, masses):
            res = 'n' if pos == 0 else peptide[pos - 1]
            res, pos, mass = self._correct_mass(res, pos, mass)
            for p, m in zip(pos, mass):
                mod_strings.append("{}={:.6f}".format(p if p == -100 else p-1, m))
        return ",".join(mod_strings)

    def _passes_scoring(self, score):
        if self.score_threshold is None:
            return True

        elif score is None:
            return False

        else:
            return (self.score_threshold - score) * (-1 ** self.score_higher_better) > 0
 
    def _generate_tsv_entries(self):
        for match in self.mod_data:
            for ind in range(len(match['peptides'])):
                mod_string = self._format_mods(match['peptides'][ind],
                                               match['mod_positions'][ind],
                                               match['mod_masses'][ind])
                score = match['scores'][ind]
                if self.score_func is not None and score is not None:
                    score = self.score_func(score)

                if len(mod_string) == 0 or not self._passes_scoring(score):
                    continue

                fields = [self.spec_file_name,
                          match['scans'][ind],
                          match['charge_states'][ind],
                          score,
                          match['peptides'][ind],
                          mod_string]
                yield "\t".join([str(e) for e in fields])
                                     

    def to_tsv(self, filename, sep="\t"):
        with open(filename, "w") as dest:
            dest.write("\t".join(["srcFile",
                                  "scanNum",
                                  "charge",
                                  "PSMscore",
                                  "peptide",
                                  "modSites"]) + "\n")
            for entry in self._generate_tsv_entries():
                dest.write(entry + "\n")
            

parser = ModificationParser("09143.mzML",
                            "data/09143.target.pep.xml",
                            "pepXML",
                            mod_masses = {"n": 42.010565,
                                          "M": 15.9949,
                                          "S": 79.966331,
                                          "T": 79.966331,
                                          "Y": 79.966331,
                                          "C": 57.021464},
                            score_string="percolator_PEP",
                            score_threshold=0.,
                            score_higher_better=True,
                            score_func=lambda x: 1-x)
parser.extract()
parser.to_tsv("./09143.input")
