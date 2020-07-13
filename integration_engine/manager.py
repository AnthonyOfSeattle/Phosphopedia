import csv
import json

from multiprocessing import Pool
from itertools import count, chain
from sqlalchemy import create_engine, update, delete, and_, func
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.sql import text

from .fdr import *
from .psm_mapper import *
from .peptide_mapper import *
from .protein_coverage_analyzer import *
from .modification_mapper import *
from .database_schema import *

class SubintegrationManager:
    def __init__(self, nworkers = 8, file_chunk_size = 1e3, record_chunk_size = 1e4):
        # Processing parameters
        self.nworkers = nworkers
        self.file_chunk_size = int(file_chunk_size)
        self.record_chunk_size = int(record_chunk_size)

        # Initialize in memory database
        self.engine = create_engine("sqlite://")
        #self.engine = create_engine(self.db_path)
        #self.session = sessionmaker(bind=self.engine)()
        TempBase.metadata.create_all(self.engine)

        # Storage for protein accessions
        self.prot_acc_map = {}
        self.protein_sequence_map = {}

    def _get_session(self):
        Session = scoped_session(sessionmaker(bind=self.engine))
        return Session

    def _extract_protein_links(self, unormalized_proteins):
        psm_protein_links = []
        for psm_id, protein_acc_list in unormalized_proteins:
            for protein_acc in protein_acc_list:
                self.prot_acc_map.setdefault(protein_acc, len(self.prot_acc_map) + 1)
                psm_protein_links.append(
                    dict(psm_id = psm_id, prot_id = self.prot_acc_map[protein_acc])
                )

        return psm_protein_links

    def map_psms(self, scan_info_files, psm_files, localization_files, group_file = None, group_number = -1):
        # Initialize workers with group information
        if group_file is not None:
            with open(group_file, "r") as group_src:
                protein_groups = json.load(group_src)
        else:
            protein_groups = {}

        workers = Pool(self.nworkers,
                       PSMMapper.initialize,
                       [protein_groups, group_number])

        # Iterate through files in chunks
        psm_counter = count()
        session = self._get_session()
        for chunk_start in range(0, len(psm_files), self.file_chunk_size):
            chunk_end = min(len(psm_files), chunk_start + self.file_chunk_size)
            print("Working on files {} through {}".format(chunk_start, chunk_end),
                  flush=True)

            psms = list(chain.from_iterable(
                workers.map(
                   PSMMapper.run, zip(
                       scan_info_files[chunk_start:chunk_end],
                       psm_files[chunk_start:chunk_end],
                       localization_files[chunk_start:chunk_end]
                   )
                )
            ))

            [p.setdefault("id", pind) for pind, p in zip(psm_counter, psms)]
            psm_protein_links = self._extract_protein_links(
                [(p["id"], p.pop("proteins")) for p in psms]
            )

            session.execute(PSM.__table__.insert(), psms)
            session.execute(PSMProtein.__table__.insert(), psm_protein_links)

        print("Got this many PSMS:", session.query(PSM).count())
        # Write protein entries for each PSM
        proteins = [dict(id = prot_id, accession = acc) for acc, prot_id in self.prot_acc_map.items()]
        for prot in proteins:
            if prot["accession"].startswith("decoy_"):
                prot["label"] = "decoy"
                prot["accession"] = prot["accession"].lstrip("decoy_")
            else:
                prot["label"] = "target"

        session.execute(Protein.__table__.insert(), proteins)
        session.commit()
        session.remove()

    def map_peptides(self):
        db_path = "sqlite://"
        workers = Pool(self.nworkers, PeptideMapper.establish_connection, ("sqlite://",))

        mapped_results = []
        session = self._get_session()
        hash_values = [q[0] for q in session.query(PSM.hash_value).distinct().order_by(PSM.hash_value).all()]
        for chunk in range(0, len(hash_values), self.record_chunk_size):
            lower_bound = hash_values[chunk]
            upper_bound = hash_values[min(chunk + self.record_chunk_size - 1, len(hash_values) - 1)]
            print("Working on peptides {} through {}".format(
                chunk, min(chunk + self.record_chunk_size - 1, 
                           len(hash_values) - 1)
            ))

            # Query a chunk of psms
            psm_objects = session.query(PSM).filter(
                and_(PSM.hash_value >= lower_bound, PSM.hash_value <= upper_bound)
            ).order_by(PSM.hash_value).all()
            psms = [p.__dict__ for p in psm_objects]

            # Group psms and analyze
            psms = [list(group) for key, group in groupby(psms, key = lambda p : p["hash_value"])]
            mapped_results.extend(workers.map(PeptideMapper.run, psms))

        # Flush peptides 
        peptides = list(chain.from_iterable(mapped_results))
        psm_to_pep = []
        for pep in peptides:
            for psm_id in pep.pop("psm_ids"):
                psm_to_pep.append({"pep_id": pep["psm_id"],
                                   "psm_id": psm_id})

        session.execute(PSMPeptide.__table__.insert(), psm_to_pep)
        session.execute(Peptide.__table__.insert(), peptides)
        session.commit()

    def read_in_fasta(self, db_file):
        # Read in proteins
        cur_acc = ""
        cur_seq = ""
        with open(db_file, "r") as src:
            for line in src:
                line = line.split()
                if line[0][0] == ">":
                    if cur_seq:
                        self.protein_sequence_map[cur_acc] = cur_seq
                    cur_acc = "|".join((line[0].lstrip(">")).split("|")[:-1])
                    cur_seq = ""
                else:
                    cur_seq += line[0]
            self.protein_sequence_map[cur_acc] = cur_seq

    def infer_protein_coverage(self):
        session = self._get_session()
        matches = (session.query(Protein.id,
                                 Protein.accession,
                                 Protein.label,
                                 Peptide.base_sequence)
                       .join(PSMProtein, Protein.id == PSMProtein.prot_id)
                       .join(Peptide, PSMProtein.psm_id == Peptide.psm_id)
                       .order_by(Protein.id)).all()

        grouped_matches = groupby(matches, lambda m : (m[0], m[1], m[2]))
        protein_dicts = [dict(id = group_info[0],
                              acc = group_info[1],
                              label = group_info[2],
                              peptide_list = [m[3] for m in match_iter],
                              seq = self.protein_sequence_map[group_info[1]])
                         for group_info, match_iter in grouped_matches]

        workers = Pool(self.nworkers)
        mapped_results = workers.map(ProteinCoverageAnalyzer.run, protein_dicts)
        update_dict = {prot_id : coverage for prot_id, coverage in mapped_results}

        max_prot_id = session.query(func.max(Protein.id)).all()[0][0]
        for chunk in range(0, max_prot_id, 1000):
            prot_objects = (session.query(Protein)
                               .filter(Protein.id.between(chunk, chunk + 1000))
                               .all())

            for prot in prot_objects:
                prot.coverage = update_dict.get(prot.id, 0.0)

            session.flush()

        session.commit()
        session.remove()

    def drop_low_coverage(self):
        session = self._get_session()
        n_connections = session.query(PSMProtein).count()
        session.execute(text(
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
        session.execute(delete(PSMProtein).where(PSMProtein.id <= n_connections))
        session.commit()
        session.remove()

    def map_modifications(self):
        session = self._get_session()
        matches = (session.query(Protein.id,
                                 Protein.accession,
                                 Protein.label,
                                 Peptide.base_sequence,
                                 Peptide.score,
                                 Peptide.modifications)
                       .join(PSMProtein, Protein.id == PSMProtein.prot_id)
                       .join(Peptide, PSMProtein.psm_id == Peptide.psm_id)
                       .order_by(Protein.id)).all()

        grouped_matches = groupby(matches, lambda m : (m[0], m[1], m[2]))
        protein_dicts = [dict(id = group_info[0],
                              acc = group_info[1],
                              label = group_info[2],
                              peptide_list = [m[3:] for m in match_iter],
                              seq = self.protein_sequence_map[group_info[1]])
                         for group_info, match_iter in grouped_matches]

        workers = Pool(self.nworkers)
        mapped_results = workers.map(ModificationMapper.run, protein_dicts)

        ptms = []
        for prot_id, ptm_dict in mapped_results:
            for pos, score in ptm_dict.items():
                ptms.append(
                    dict(prot_id = prot_id, position=pos, score=score)
                )
        session.execute(PTM.__table__.insert(), ptms)
        session.commit()
        session.remove()

    def update_peptide_fdr(self):
        session = self._get_session()
        peptide_scores = session.query(Peptide.psm_id, Peptide.score, Peptide.label).all()
        peptide_scores.sort(key=lambda pep: -pep[1])

        score_array = np.array([pep[1] for pep in peptide_scores])
        label_array = np.array([pep[2] for pep in peptide_scores])
        fdr_array = calculate_fdr(score_array, label_array)
        update_dict = {pep[0] : fdr for pep, fdr in zip(peptide_scores, fdr_array)}

        max_psm_id = session.query(func.max(Peptide.psm_id)).all()[0][0]
        for chunk in range(0, max_psm_id, 1000):
            pep_objects = (session.query(Peptide)
                               .filter(Peptide.psm_id.between(chunk, chunk + 1000))
                               .all())

            for pep in pep_objects:
                pep.qvalue = update_dict[pep.psm_id]

            session.flush()

        session.commit()
        session.remove()

    def update_ptm_fdr(self):
        session = self._get_session()
        ptm_scores = (session.query(PTM.prot_id,
                                    PTM.position,
                                    PTM.score,
                                    Protein.label)
                          .join(Protein, PTM.prot_id == Protein.id)).all()
        ptm_scores.sort(key=lambda ptm: -ptm[2])

        score_array = np.array([ptm[2] for ptm in ptm_scores])
        label_array = np.array([ptm[3] for ptm in ptm_scores])
        fdr_array = calculate_fdr(score_array, label_array)
        update_dict = {(ptm[0], ptm[1]) : fdr for ptm, fdr in zip(ptm_scores, fdr_array)}

        max_prot_id = session.query(func.max(PTM.prot_id)).all()[0][0]
        for chunk in range(0, max_prot_id, 1000):
            ptm_objects = (session.query(PTM)
                               .filter(PTM.prot_id.between(chunk, chunk + 1000))
                               .all()
            )

            for ptm in ptm_objects:
                ptm.qvalue = update_dict[(ptm.prot_id, ptm.position)]

            session.flush()

        session.commit()
        session.remove()

    def dump(self, path = ""):
        session = self._get_session()

        # Dump psms
        psm_list = (session.query(PSM.id,
                                  PSM.sample_name,
                                  PSM.scan_number,
                                  PSM.scan_rt,
                                  PSM.precursor_mz,
                                  PSM.precursor_charge,
                                  PSM.score,
                                  PSMPeptide.pep_id,
                                  PSM.qvalue,
                                  PSM.label)
                        .join(PSMPeptide, PSM.id == PSMPeptide.psm_id)
                   ).all()

        psm_list = [",".join([str(e) for e in psm]) + "\n" for psm in psm_list]
        with open(os.path.join(path, "psms.csv"), 'w') as dest:
            dest.write(",".join(["id", "sample_name", 
                                 "scan_number", "scan_rt",
                                 "precursor_mz", "precursor_charge",
                                 "score", "pep_id", 
                                 "qvalue", "label"]) + "\n")
            dest.writelines(psm_list)

        # Dump peptides and associated PTMs
        # The setup here allows PTMs to be associated at the peptide level
        pep_mod_list = (session.query(Peptide.psm_id,
                                      Peptide.sequence,
                                      Peptide.score,
                                      Peptide.qvalue,
                                      Peptide.label,
                                      Peptide.modifications)
                        ).all()

        pep_list = []
        mod_list = []
        for entry in pep_mod_list:
            pep_list.append(entry[:-1])
            mod = pickle.loads(entry[-1])
            for pos, info in mod.items():
                mod_list.append((entry[0], pos, info["score"], info["residue"], info["mass"]))

        pep_list = [",".join([str(e) for e in pep]) + "\n" for pep in pep_list]
        with open(os.path.join(path, "peptides.csv"), 'w') as dest:
            dest.write(",".join(["id", "sequence", "score", "qvalue", "label"]) + "\n")
            dest.writelines(pep_list)

        mod_list = [",".join([str(e) for e in mod]) + "\n" for mod in mod_list]
        with open(os.path.join(path, "peptide_modifications.csv"), 'w') as dest:
            dest.write(",".join(["pep_id", "pos", "score", "residue", "mass", "label"]) + "\n")
            dest.writelines(mod_list)

        # Dump peptide to protein links
        junction_list = (session.query(PSMProtein.psm_id,
                                       PSMProtein.prot_id)
                            .join(Peptide, PSMProtein.psm_id == Peptide.psm_id)
                        ).all()

        junction_list = [",".join([str(e) for e in link]) + "\n" for link in junction_list]
        with open(os.path.join(path, "peptide_protein.csv"), 'w') as dest:
            dest.write(",".join(["pep_id", "prot_id"]) + "\n")
            dest.writelines(junction_list)

        # Dump proteins
        prot_list = (session.query(Protein.id,
                                   Protein.accession,
                                   Protein.label)
                         .join(PTM, Protein.id == PTM.prot_id)
                         .group_by(Protein.id, Protein.accession)
                   ).all()

        prot_list = [",".join([str(e) for e in prot]) + "\n" for prot in prot_list]
        with open(os.path.join(path, "proteins.csv"), 'w') as dest:
            dest.write(",".join(["id", "accession", "label"]) + "\n")
            dest.writelines(prot_list)

        # Dump modified sites
        site_list = (session.query(PTM.prot_id,
                                   PTM.position,
                                   PTM.score,
                                   PTM.qvalue)
                         .join(Protein, PTM.prot_id == Protein.id)
                   ).all()

        site_list = [",".join([str(e) for e in site]) + "\n" for site in site_list]
        with open(os.path.join(path, "sites.csv"), 'w') as dest:
            dest.write(",".join(["prot_id", "position", "score", "qvalue"]) + "\n")
            dest.writelines(site_list)

        session.remove()
