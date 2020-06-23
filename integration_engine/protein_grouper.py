import sys
import re
import json
import numpy as np
import pandas as pd
from itertools import product, cycle

class InputError(Exception):
    def __init__(self, message):
        self.message = message


class Digester:
    def __init__(self, enzyme, n_missed = 2, min_len=6, max_len=60):
        # This should be all that you need to change to add a new enzyme
        # enzyme names <- use lowercase
        # amino acids <- use UPPERCASE
        enzyme_regex = {"trypsin" : "[KR](?!P)",
                        "lysc" : "K"}
        self.cut_site = re.compile(enzyme_regex[enzyme.lower()])
        self.n_missed = n_missed
        self.min_len = min_len
        self.max_len = max_len

    def _build_missed_cleavages(self, seq_list, group_size):
        if group_size > len(seq_list):
            return []

        args = [iter(seq_list) for i in range(group_size)]
        for ind, it in enumerate(args):
            for i in range(ind):
                next(it)

        return ["".join(peps) for peps in zip(*args)]

    def _size_select(self, peptides):
        select = [len(p) >= self.min_len and len(p) <= self.max_len for p in peptides]
        return peptides[select]

    def digest(self, sequence):
        sequence_digest = [[]]

        last = 0
        for pep in re.finditer(self.cut_site, sequence):
            sequence_digest[0].append(sequence[ last:pep.span()[1] ])
            last = pep.span()[1]

        if last != len(sequence):
            sequence_digest[0].append(sequence[ last: ])
        
        for i in range(self.n_missed):
            sequence_digest.append(self._build_missed_cleavages(sequence_digest[0], i + 2))
        
        sequence_digest = self._size_select(np.concatenate(sequence_digest))
        return sequence_digest


class PeptideExtractor:
    def __init__(self, file_name, enzyme):
        self.file_name = file_name
        self.digester = Digester(enzyme)
        self.peptide_list = []
        self.id_list = []

    def _process_fasta_entry(self, sequence_buffer):
        # Check if proper fasta sequence
        if (">" != sequence_buffer[0][0] or
                len(sequence_buffer[0]) == 1):
            raise InputError("FASTA record lacking header")

        if len(sequence_buffer) == 1:
            raise InputError("FASTA record lacking sequence")

        sequence_id = sequence_buffer[0].split()[0][1:]

        peptides = self.digester.digest("".join(sequence_buffer[1:]))

        self.peptide_list.extend(peptides)
        self.id_list.extend([sequence_id] * len(peptides))

    def zip(self):
        return zip(self.id_list, self.peptide_list)

    def extract(self):
        with open(self.file_name) as source:
            sequence_buffer = []
            for line in source:
                line = line.rstrip()
                if line[0] == ">" and len(sequence_buffer) > 0:
                    self._process_fasta_entry(sequence_buffer)
                    sequence_buffer = [line]

                else:
                    sequence_buffer.append(line)

            self._process_fasta_entry(sequence_buffer)


class ProteinGrouper:
    def __init__(self, extractor, n_groups, n_iter=10):
        self.extractor = extractor
        self.n_groups = n_groups
        self.n_iter = n_iter
        
        self.inverted_index = None
        self.protein_matches = None
        self.group_map = None
        self.group_var = np.inf

    def _create_inverted_index(self):
        inverted_index = {}
        for ref, pep in self.extractor.zip():
            inverted_index.setdefault(pep, set())
            inverted_index[pep].add(ref)

        return inverted_index

    def _match_proteins(self):
        protein_matches = {}
        for ref_set in self.inverted_index.values():
            for ref1, ref2 in product(ref_set, ref_set):
                protein_matches.setdefault(ref1, set())
                protein_matches[ref1].add(ref2)

        return protein_matches

    def _create_groups(self):
        references = list(self.protein_matches.keys())
        np.random.shuffle(references)

        group_cycle = cycle(range(self.n_groups))
        group_map = {}
        for ref in references:
            ref_set = self.protein_matches[ref]
            group = None
            for ref in ref_set:
                group = group_map.get(ref, None)
    
            if group is None:
                group = next(group_cycle)
    
            group_map.update({ref : group for ref in ref_set})

        _, group_counts = np.unique(list(group_map.values()), return_counts=True)
        return group_map, np.var(group_counts)

    def group(self):
        self.extractor.extract()
        self.inverted_index = self._create_inverted_index()
        self.protein_matches = self._match_proteins()

        for i in range(self.n_iter):
            group_map, group_var = self._create_groups()
            if group_var < self.group_var:
                self.group_var = group_var
                self.group_map = group_map

    def to_json(self, file_name, **kwargs):
        with open(file_name, "w") as dest:
            json.dump(self.group_map, dest, **kwargs)

if __name__ == "__main__":
    fasta_name = sys.argv[1]
    destination = sys.argv[2]
    n_groups = int(sys.argv[3])
    pg = ProteinGrouper(
            PeptideExtractor(fasta_name, "trypsin"),
            n_groups
         )
    pg.group()
    pg.to_json(destination)
