import numpy as np

class ProteinCoverageAnalyzer:
    @staticmethod
    def get_peptide_ranges(peptides, seq):
        peptide_ranges = []
        for upep in peptides:
            start = seq.find(upep)
            while start >= 0:
                peptide_ranges.append((start, start + len(upep)))
                start = seq.find(upep, start + 1)
        return sorted(peptide_ranges)

    @staticmethod
    def calculate_coverage(ranges):
        last_start = 0
        last_end = 0
        total_coverage = 0
        for start, end in ranges:
            if start >= last_end:
                total_coverage += last_end - last_start
                last_start = start
            last_end = max(last_end, end)
        total_coverage += last_end - last_start

        return total_coverage

    @staticmethod
    def run(prot):
        # Deal with peptide sequences and reverse if decoys
        unique_peptides = np.unique(prot['peptide_list'])
        if prot["label"] == "decoy":
            unique_peptides = np.array([pep[:-1][::-1] + pep[-1] for pep in unique_peptides])

        peptide_ranges = ProteinCoverageAnalyzer.get_peptide_ranges(unique_peptides, prot["seq"])
        total_coverage = ProteinCoverageAnalyzer.calculate_coverage(peptide_ranges)
        covered_proportion = total_coverage / len(prot["seq"])
        return prot["id"], covered_proportion
