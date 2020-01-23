import os
import re
import time
import shutil
import pickle
import numpy as np
import pandas as pd
from glob import glob
from multiprocessing import Pool
from itertools import groupby
from numpy import inf

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine, Column, Integer, Numeric, String, JSON, BLOB, update, delete
from sqlalchemy.sql import text
from sqlalchemy.orm import sessionmaker, joinedload

#############################
#                           #
# Temporary database schema #
#                           #
#############################

TempBase = declarative_base()

class PSM(TempBase):

    __tablename__ = "psms"

    hash_value = Column(Integer, index=True)

    id = Column(Integer, primary_key=True)
    sample_name = Column(String)
    scan_number = Column(Integer)
    label = Column(String)
    base_sequence = Column(String)
    score = Column(Numeric(asdecimal=False))
    qvalue = Column(Numeric(asdecimal=False))
    pep = Column(Numeric(asdecimal=False))
    modifications = Column(BLOB)

class Peptide(TempBase):

    __tablename__ = "peptides"

    psm_id = Column(Integer, primary_key = True) #index=True)

    #id = Column(Integer, primary_key=True)
    label = Column(String)
    base_sequence = Column(String)
    sequence = Column(String)
    score = Column(Numeric(asdecimal=False))
    qvalue = Column(Numeric(asdecimal=False))
    pep = Column(Numeric(asdecimal=False))
    modifications = Column(BLOB)

class PSMProtein(TempBase):

    __tablename__ = "psms_proteins"

    prot_id = Column(String, index=True)
    psm_id = Column(Integer, index=True)
    id = Column(Integer, primary_key=True)

class Protein(TempBase):
    
    __tablename__ = "proteins"

    id = Column(Integer, primary_key=True)
    accession = Column(String)
    coverage = Column(Numeric(asdecimal=False))

##########################
#                        #
# PSM processing classes #
#                        #
##########################

class PSMMapper:

    @staticmethod
    def gather_data(psm_path, localization_path):
        localizations = pd.read_csv(localization_path, sep="\t").set_index("Scan")
        psm_scores = pd.read_csv(psm_path, sep="\t",
                                 usecols=["scan",
                                          "percolator score",
                                          "percolator q-value",
                                          "percolator PEP",
                                          "protein id"]).set_index("scan", drop=False)
        localized_scores = psm_scores.join(localizations, how="left")
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
                                             label=os.path.split(file_tuple[0])[1].split(".")[2],
                                             basename=os.path.split(file_tuple[0])[1].split(".")[0]).tolist()
        return psm_entries

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
        best_peptide = next(peptide_group[1])
        for pep in peptide_group[1]:
            if pep["score"] > best_peptide["score"]:
                best_peptide = pep

        return best_peptide

    @staticmethod
    def merge_psms(records):
        localized_ptms = PeptideMapper.get_localized_modifications(records)
        unmerged_peptides = [PeptideMapper.build_peptide(psm, localized_ptms) for psm in records]
        unmerged_peptides.sort(key=lambda p: p["sequence"])
        peptide_groups = groupby(unmerged_peptides, key=lambda p: p["sequence"])
        return list(map(PeptideMapper.collapse_peptides, peptide_groups))

    @staticmethod
    def run(key):
        objects = SESSION.query(PSM).filter(PSM.hash_value == key).all()
        records = [o.__dict__ for o in objects] 
        for r in records:
            r["modifications"] = pickle.loads(r["modifications"])

        if len(np.unique([r["base_sequence"] for r in records])) > 1:
            print("Hash collision!")

        merged_psms = PeptideMapper.merge_psms(records)
        return merged_psms


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
        if prot["acc"].startswith("decoy"):
            unique_peptides = np.array([pep[:-1][::-1] + pep[-1] for pep in unique_peptides])

        peptide_ranges = ProteinCoverageAnalyzer.get_peptide_ranges(unique_peptides, prot["seq"])
        total_coverage = ProteinCoverageAnalyzer.calculate_coverage(peptide_ranges)
        return prot["id"], total_coverage

###########################
#                         #
# Statistical calculation # 
#                         #
###########################

def calculate_fdr(score_array, label_array):
    # Arguments should be sorted such that best scores come first
    # This function is very similar to that provided by Percolator
    qvalues = []

    score_array = np.atleast_1d(score_array)
    label_array = np.atleast_1d(label_array)
    ntargets = np.cumsum(label_array == "target")
    ndecoys = np.arange(1, ntargets.shape[0] + 1) - ntargets
    qvalues = (ndecoys + 1)/ntargets
    for ind in np.arange(qvalues.shape[0])[::-1]:
         if ind - 1 >= 0 and score_array[ind] == score_array[ind-1]:
             qvalues[ind - 1] = qvalues[ind]

         if ind + 1 < qvalues.shape[0]:
             qvalues[ind] = min(qvalues[ind], qvalues[ind + 1])

    return qvalues

##########################
#                        #
# Main integration class # 
#                        #
##########################

class IntegrationManager:
    def __init__(self, temp_path, nworkers, write_buffer_size = 1e5, overwrite=True):
        self.nworkers = nworkers
        self.temp_path = os.path.join(temp_path, "integration_engine")

        if overwrite and os.path.exists(self.temp_path):
            shutil.rmtree(self.temp_path)

        if not os.path.exists(self.temp_path):
            os.makedirs(self.temp_path, mode=0o700)

        self.write_buffer_size = write_buffer_size
        self.db_path = "sqlite:///" + os.path.join(self.temp_path, "psm.db")
        self.engine = create_engine(self.db_path)
        self.session = sessionmaker(bind=self.engine)()
        TempBase.metadata.create_all(self.engine)

        self.prot_acc_map = {}

    def extract_protein_links(self, unormalized_proteins):
        psm_protein_links = []
        for psm_id, protein_acc_list in unormalized_proteins:
            for protein_acc in protein_acc_list:
                self.prot_acc_map.setdefault(protein_acc, len(self.prot_acc_map) + 1)
                psm_protein_links.append(
                    dict(psm_id = psm_id, prot_id = self.prot_acc_map[protein_acc])
                )

        return psm_protein_links

    def map_psms(self, psm_files, localization_files):
        workers = Pool(self.nworkers)
        n_total_psms = 0
        for ind in range(0, len(psm_files), 50):
            mapped_results = workers.map(
                PSMMapper.run, zip(
                    psm_files[ind:min(len(psm_files), ind + 50)], 
                    localization_files[ind:min(len(psm_files), ind + 50)]
                )
            )
            psms = np.concatenate(mapped_results)
            [p.setdefault("id", pind) for pind, p in enumerate(psms, n_total_psms + 1)]
            n_total_psms += psms.shape[0] 
            psm_protein_links = self.extract_protein_links(
                [(p["id"], p.pop("proteins")) for p in psms]
            )
            n_total_psms += psms.shape[0]

            self.session.bulk_insert_mappings(PSM, psms)
            self.session.bulk_insert_mappings(PSMProtein, psm_protein_links)
            self.session.commit()

        proteins = [dict(id = prot_id, accession = acc) for acc, prot_id in self.prot_acc_map.items()]
        self.session.bulk_insert_mappings(Protein, proteins)
        self.session.commit()

    def map_peptides(self):
        hash_values = [q[0] for q in self.session.query(PSM.hash_value).distinct().all()]
        workers = Pool(self.nworkers, PeptideMapper.establish_connection, (self.db_path,))
        mapped_results = workers.map(PeptideMapper.run, hash_values)
        peptides = np.concatenate(mapped_results)
        self.session.bulk_insert_mappings(Peptide, peptides)
        self.session.commit()

    def infer_protein_coverage(self, db_file):
        # Read in proteins
        cur_acc = ""
        cur_seq = ""
        protein_sequence_map = {}
        with open(db_file, "r") as src:
            for line in src:
                line = line.split()
                if line[0][0] == ">":
                    if cur_seq:
                        protein_sequence_map[cur_acc] = cur_seq
                    cur_acc = "|".join((line[0].lstrip(">")).split("|")[:-1])
                    cur_seq = ""
                else:
                    cur_seq += line[0]
            protein_sequence_map[cur_acc] = cur_seq

        matches = (self.session.query(Protein.id, Protein.accession, Peptide.base_sequence)
                       .join(PSMProtein, Protein.id == PSMProtein.prot_id)
                       .join(Peptide, PSMProtein.psm_id == Peptide.psm_id)
                       .order_by(Protein.id)).all()
        grouped_matches = groupby(matches, lambda m : (m[0], m[1]))
        protein_dicts = [dict(id = prot_id, acc = prot_acc, 
                              peptide_list = [m[2] for m in match_iter], 
                              seq = protein_sequence_map[prot_acc.lstrip("decoy_")])
                         for (prot_id, prot_acc), match_iter in grouped_matches]
        workers = Pool(self.nworkers)
        mapped_results = workers.map(ProteinCoverageAnalyzer.run, protein_dicts)

        for prot_id, coverage in mapped_results:
            self.session.execute(update(Protein).where(Protein.id == prot_id).values(coverage = coverage))
        self.session.commit()

    def drop_low_coverage(self):
        n_connections = self.session.query(PSMProtein).count()
        self.session.execute(text(
            """
            WITH ranked_coverage AS (
              SELECT pp.prot_id as prot_id, pp.psm_id as psm_id,
                     RANK() OVER(
                       PARTITION BY pp.psm_id
                       ORDER BY pro.coverage, pro.accession
                     ) rank
              FROM psms_proteins pp
              LEFT JOIN proteins pro
              ON pp.prot_id = pro.id
              ORDER BY pp.prot_id, pp.psm_id
            )
            INSERT INTO psms_proteins (psm_id, prot_id)
            SELECT psm_id, prot_id
            FROM ranked_coverage
            WHERE rank = 1;
            """
        ))
        self.session.execute(delete(PSMProtein).where(PSMProtein.id <= n_connections))
        self.session.commit()

    def update_peptide_fdr(self):
        peptide_scores = self.session.query(Peptide.id, Peptide.score, Peptide.label).all()
        #peptide_scores = self.session.query(PSM.id, PSM.score, PSM.label, PSM.qvalue).all()
        peptide_scores.sort(key=lambda pep: -pep[1])

        score_array = np.array([pep[1] for pep in peptide_scores])
        label_array = np.array([pep[2] for pep in peptide_scores])
        print(np.sum(calculate_fdr(score_array, label_array) < 0.01))

if __name__ == "__main__":
    print("Finding files")
    percolator_files = glob("../../test/percolator/*/*/*.percolator.*.psms.txt")
    percolator_files.sort(
        key=lambda path: (os.path.split(path)[1].split(".")[0], os.path.split(path)[1].split(".")[2])
    )
    ascore_files = glob("../../test/ascore/*/*/*.ascore.*.txt")
    ascore_files.sort(
        key=lambda path: (os.path.split(path)[1].split(".")[0], os.path.split(path)[1].split(".")[2])
    )

    print("Processing")
    manager = IntegrationManager(os.getenv("TMPDIR"), 8, overwrite=1)

    t0 = time.time()
    manager.map_psms(percolator_files[:100], ascore_files[:100])
    print("PSMs took {} seconds".format(time.time() - t0))

    t0 = time.time()
    manager.map_peptides()
    print("Peptides took {} seconds".format(time.time() - t0))

    t0 = time.time()
    db_file = "../../test/config/sp_iso_HUMAN_4.9.2015_UP000005640.fasta" 
    manager.infer_protein_coverage(db_file)
    print("Coverage estimation took {} seconds".format(time.time() - t0))

    t0 = time.time()
    manager.drop_low_coverage()
    print("High coverage selection took {} seconds".format(time.time() - t0))

    t0 = time.time()
    #manager.update_peptide_fdr()
    print("Peptide FDR update took {} seconds".format(time.time() - t0))
