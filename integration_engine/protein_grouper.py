import sys
import re
import json
import numpy as np
import pandas as pd
from itertools import product, cycle, combinations, chain
from collections import deque
from multiprocessing import Pool

class PercolatorReader:
    @staticmethod
    def run(psm_path, fdr_filter):
        psm_map = pd.read_csv(
            psm_path, sep="\t",
            usecols=["scan", "percolator q-value", "protein id"]
            ).set_index("scan")
        psm_map = psm_map[psm_map["percolator q-value"] < fdr_filter]

        protein_groups = []
        for prot_list in psm_map["protein id"]:
            prot_list = prot_list.split(",")
            #prot_list = ["|".join(p.split("|")[:-1]) for p in prot_list]
            isdecoy = np.array([p.startswith("decoy_") for p in prot_list])
            if not (np.all(isdecoy) or np.all(~isdecoy)):
                continue

            protein_groups.append(prot_list)

        return protein_groups

class PercolatorProteinGrouper:
    def __init__(self, n_groups, n_workers, chunk_size, fdr_filter=1.):
        self.n_groups = n_groups
        self.n_workers = n_workers
        self.chunk_size = int(chunk_size)
        self.fdr_filter = fdr_filter
        self.protein_matches = {}
        self.group_map = {}

    def _match_proteins(self, groups):
        for ref_list in groups:
            for ref in ref_list:
                self.protein_matches.setdefault(ref, set())
                self.protein_matches[ref].update(ref_list)

    def _get_clique(self, ref):
        clique = set()
        ref_stack = [ref]
        while ref_stack:
            ref = ref_stack.pop()
            for ref in self.protein_matches[ref]:
                if ref not in clique:
                    clique.add(ref)
                    ref_stack.append(ref)

        return clique

    def _create_groups(self):
        references = self.protein_matches.keys()
        group_count = {group : 0 for group in range(self.n_groups)}
        group_cycle = cycle(group_count.keys())
        for ref in references:
            if ref not in self.group_map:
                clique = self._get_clique(ref)

                group = next(group_cycle)
                this_count = group_count[group]
                next_count = group_count[(group + 1) % self.n_groups]
                while this_count > next_count:
                    group = next(group_cycle)
                    this_count = group_count[group]
                    next_count = group_count[(group + 1) % self.n_groups]

                self.group_map.update({ref : group for ref in clique})
                group_count[group] += len(clique)

    def group(self, files):
        workers = Pool(self.n_workers)
        for chunk_start in range(0, len(files), self.chunk_size):
            chunk_end = min(len(files), chunk_start + self.chunk_size)
            print("Working on files {} through {}".format(chunk_start, chunk_end),
                  flush=True)

            protein_groups = list(chain.from_iterable(
                workers.starmap(PercolatorReader.run,
                                [(f, self.fdr_filter)
                                 for f in files[chunk_start:chunk_end]])
                ))
            self._match_proteins(protein_groups)

        self._create_groups()

    def to_json(self, file_name, **kwargs):
        with open(file_name, "w") as dest:
            json.dump(self.group_map, dest, **kwargs) 

