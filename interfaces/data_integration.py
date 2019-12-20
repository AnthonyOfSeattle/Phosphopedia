import sys
import warnings
from multiprocessing import Process, Manager
from sqlalchemy.exc import OperationalError, SAWarning
from queue import Empty
from os import getpid

from database_schema import *
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, joinedload
from itertools import groupby
from numpy import inf

class PSMMerger:
    def __init__(self, db_path, 
                 inputq, outputq,
                 term_sig, hold_sig):
        self.db_path=db_path
        self.inputq = inputq
        self.outputq = outputq
        self.term_sig = term_sig
        self.hold_sig = hold_sig

    @staticmethod
    def get_records(seq, session):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=SAWarning)
            records = session.query(PSM).join(PSMModification).filter(PSM.base_sequence == seq).all()
        return records

    @staticmethod
    def get_localized_modifications(records):
        localized = {}

        all_modifications = []
        [all_modifications.extend(psm.modifications) for psm in records]

        for mod in all_modifications:
            localized[mod.position] = max(localized.get(mod.position, 0.),
                                          mod.localization_score)
            if mod.localization_score < 13:
                [localized.setdefault(int(pos), 0.) for pos in mod.alternative_positions.split(",")]

        return localized

    @staticmethod
    def build_peptide(psm, localized_ptms):
        tokenized_seq = list(psm.base_sequence)
        occupied_pos = {m.position : True for m in psm.modifications}
        peptide_mods = []
        for mod in psm.modifications:
            mod = PeptideModification(
                position = mod.position,
                residue = mod.residue,
                mass = mod.mass,
                localization_score = mod.localization_score,
                alternative_positions = mod.alternative_positions
            )

            if mod.localization_score < 13:
                occupied_pos.pop(mod.position)
                alternative_positions = [int(pos) for pos in mod.alternative_positions.split(",")]
                for pos in alternative_positions:
                    if localized_ptms[pos] > mod.localization_score and not occupied_pos.get(pos, False):
                        mod.position = pos
                        mod.localization_score = localized_ptms[pos]

            occupied_pos.setdefault(mod.position, True)
            if mod.position == 0:
                tokenized_seq[0] = "n" + "[{:.0f}]".format(mod.mass) + tokenized_seq[0]
            else:
                tokenized_seq[mod.position - 1] = tokenized_seq[mod.position - 1] + "[{:.0f}]".format(mod.mass)

            peptide_mods.append(mod)

        new_peptide = Peptide(
            label = psm.label,
            sequence = "".join(tokenized_seq),
            score = psm.score
        )
        new_peptide.modifications = peptide_mods
        return new_peptide

    @staticmethod
    def collapse_peptides(peptide_group):
        best_peptide = next(peptide_group[1])
        for pep in peptide_group[1]:
            best_mods = []
            for modl, modr in zip(best_peptide.modifications, pep.modifications):
                if modl.localization_score == inf:
                    best_mods.append(modl)
                if modl.localization_score < modr.localization_score:
                    best_mods.append(modr)
                elif (modl.localization_score == modr.localization_score
                      and len(modl.alternative_positions) < len(modr.alternative_positions)):
                    best_mods.append(modr)
                else:
                    best_mods.append(modl)

            if pep.score > best_peptide.score:
                best_peptide.score = pep.score

            best_peptide.modifications = best_mods

        return best_peptide

    @staticmethod
    def merge_psms(records):
        localized_ptms = PSMMerger.get_localized_modifications(records)
        unmerged_peptides = [PSMMerger.build_peptide(psm, localized_ptms) for psm in records]
        unmerged_peptides.sort(key=lambda p: p.sequence)
        peptide_groups = groupby(unmerged_peptides, key=lambda p: p.sequence)
        return list(map(PSMMerger.collapse_peptides, peptide_groups))

    def run(self):
        process_engine = create_engine(self.db_path)
        process_session = sessionmaker(bind=process_engine)()

        while not self.term_sig.is_set():
            try:
                seq = self.inputq.get(True, 0.01)
                self.hold_sig.clear()
                records = PSMMerger.get_records(seq, process_session)
                merged_records = PSMMerger.merge_psms(records)
                for mr in merged_records:
                    self.outputq.put(mr)

            except Empty:
                self.hold_sig.set()

class IntegrationManager:
    def __init__(self, db_path, nmergers, batch_size):
        # Establish database connection
        self.db_path = db_path
        self.engine = create_engine(self.db_path)
        self.session = sessionmaker(bind=self.engine)()
        self.batch_size = batch_size

        # Establish process management
        self.nmergers = nmergers
        self.mergers = []
        self.manager = Manager()

        self.sequenceq = self.manager.Queue(batch_size * 2)
        self.mergedq = self.manager.Queue(batch_size * 10)

        self.term_sig = self.manager.Event()
        self.hold_sigs = [self.manager.Event() for m in range(nmergers)]

    def gather_sequences(self):
        print("-- Gathering sequences to parse --")
        return [q[0] for q in self.session.query(PSM.base_sequence).distinct().all()]

    def start_mergers(self):
        print("-- Creating and starting {} mergers --".format(self.nmergers))
        for i in range(self.nmergers):
            m = PSMMerger(self.db_path, self.sequenceq, self.mergedq,
                          self.term_sig, self.hold_sigs[i])
            self.mergers.append(Process(target=m.run))

        [m.start() for m in self.mergers]

    def flush_to_database(self):
        results = []
        while True:
            try:
                results.append(self.mergedq.get(True, 1))
            except Empty:
                break

        print("-- Got {} results --".format(len(results)))

    def shutdown(self):
        print("-- Attempting shutdown -- ")
        self.term_sig.set()
        for m in self.mergers:
            m.join()

        self.manager.shutdown()

    def run(self):
        sequences = self.gather_sequences()[:10000]
        self.start_mergers()

        for batchn, ind in enumerate(range(0, len(sequences), self.batch_size)):
            print("-- Processing batch {} --".format(batchn))
            for seq in sequences[ind:min(len(sequences), ind + self.batch_size)]:
                self.sequenceq.put(seq)

            processes_held = [False]
            while not all(processes_held) and not self.sequenceq.empty():
                processes_held = [h.is_set() for h in self.hold_sigs]

            print("-- Flushing batch {} --".format(batchn))
            self.flush_to_database()

        self.shutdown()


if __name__ == "__main__":
    try:
        manager = IntegrationManager("sqlite:////net/villen/vol2/users/valenta4/Phosphopedia/test/phosphopedia.db", 16, 1000)
        manager.run()
    except (KeyboardInterrupt, SystemExit):
        print('interrupt signal received')
        sys.exit(1)
