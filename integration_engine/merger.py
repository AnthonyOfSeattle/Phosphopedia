import os
import numpy as np
import pandas as pd

class IdCounter:
    def __init__(self):
        self.psm_count = 0
        self.protein_count = 0


class SubintegrationMerger:
    def __init__(self, dest):
        self.counter = IdCounter()
        self.dest = dest
        self.first = True

    def _dump(self, df, dest):
        if self.first:
            df.to_csv(dest, index=False)
        else:
            df.to_csv(dest, index=False, header=False, mode="a")

    def _process_sites(self, source):
        data = pd.read_csv(source)
        data.prot_id += self.counter.protein_count
        self._dump(data, os.path.join(self.dest, "sites.csv"))

    def _process_peptide_protein_links(self, source):
        data = pd.read_csv(source)
        data.pep_id += self.counter.psm_count
        data.prot_id += self.counter.protein_count
        self._dump(data, os.path.join(self.dest, "peptide_protein.csv"))

    def _process_peptide_modifications(self, source):
        data = pd.read_csv(source)
        data.pep_id += self.counter.psm_count
        self._dump(data, os.path.join(self.dest, "peptide_modifications.csv"))

    def _process_proteins(self, source):
        data = pd.read_csv(source)
        data.id += self.counter.protein_count
        self._dump(data, os.path.join(self.dest, "proteins.csv"))
        self.counter.protein_count += data.id.max()

    def _process_peptides(self, source):
        data = pd.read_csv(source)
        data.id += self.counter.psm_count
        self._dump(data, os.path.join(self.dest, "peptides.csv"))

    def _process_psms(self, source):
        data = pd.read_csv(source)
        data.id += self.counter.psm_count
        data.pep_id += self.counter.psm_count
        self._dump(data, os.path.join(self.dest, "psms.csv"))
        self.counter.psm_count = data.id.max()

    def process(self, source_list):
        source_list = np.sort(np.asarray(source_list))

        for source_path in source_list:
            # Easier for count tracking to process aux tables first
            self._process_sites(
                os.path.join(source_path, "sites.csv")
            )

            self._process_peptide_protein_links(
                os.path.join(source_path, "peptide_protein.csv")
            )

            self._process_peptide_modifications(
                os.path.join(source_path, "peptide_modifications.csv")
            )

            self._process_proteins(
                os.path.join(source_path, "proteins.csv")
            )

            self._process_peptides(
                os.path.join(source_path, "peptides.csv")
            )

            self._process_psms(
                os.path.join(source_path, "psms.csv")
            )

            self.first = False
