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
from sqlalchemy import create_engine, Column, Integer, Numeric, String, JSON, BLOB
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

    id = Column(Integer, primary_key=True)
    psm_id = Column(Integer)
    label = Column(String)
    base_sequence = Column(String)
    sequence = Column(String)
    score = Column(Numeric(asdecimal=False))
    qvalue = Column(Numeric(asdecimal=False))
    pep = Column(Numeric(asdecimal=False))
    modifications = Column(BLOB)

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
                                          "percolator PEP"]).set_index("scan", drop=False)
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

    def map_psms(self, psm_files, localization_files):
        workers = Pool(self.nworkers)
        for ind in range(0, len(psm_files), 50):
            mapped_results = workers.map(
                PSMMapper.run, zip(
                    psm_files[ind:min(len(psm_files), ind + 50)], 
                    localization_files[ind:min(len(psm_files), ind + 50)]
                )
            )
            psms = np.concatenate(mapped_results)
            self.session.bulk_insert_mappings(PSM, psms)
            self.session.commit()

    def map_peptides(self):
        hash_values = [q[0] for q in self.session.query(PSM.hash_value).distinct().all()]
        workers = Pool(self.nworkers, PeptideMapper.establish_connection, (self.db_path,))
        mapped_results = workers.map(PeptideMapper.run, hash_values)
        peptides = np.concatenate(mapped_results)
        self.session.bulk_insert_mappings(Peptide, peptides)
        self.session.commit()

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
    manager.map_psms(percolator_files, ascore_files)
    print("PSMs took {} seconds".format(time.time() - t0))

    t0 = time.time()
    manager.map_peptides()
    print("Peptides took {} seconds".format(time.time() - t0))
