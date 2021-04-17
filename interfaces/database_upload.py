import os
import re
import numpy as np
import pandas as pd
from itertools import product
from sqlalchemy import create_engine, update
from sqlalchemy.orm import sessionmaker, scoped_session
from database_schema import *

class DatabaseBuild:
    def __init__(self, build_path):
        components = [
            "psms",
            "peptides",
            "peptide_modifications",
            "peptide_protein",
            "proteins",
            "sites"
        ]

        for c in components:
            path = os.path.join(build_path, c + ".csv")
            setattr(
                self, c,
                pd.read_csv(path, na_values="None")
            )

class Uploader:
    def __init__(self,
                 host, 
                 build_path, 
                 fasta_path, 
                 annotation_path="",
                 buffer_size=1e6):
        self.engine = create_engine(host)
        self.buffer_size = int(buffer_size)
        PhosphopediaBase.metadata.create_all(self.engine)

        self.build = DatabaseBuild(build_path)
        self.fasta_path = fasta_path
        self.annotation_path = annotation_path

        self.psms = None
        self.peptides = None
        self.peptide_sites = None
        self.proteins = None
        self.sites = None
        self.pathways = None
        self.pathway_sites = None

    def _convert_psms(self):
        self.psms = self.build.psms[
            np.logical_and(
                 self.build.psms.qvalue <= 0.01,
                 self.build.psms.label == "target"
            )
        ].loc[:,
            ["id", "pep_id", "qvalue",
             "sample_name", "scan_number", 
             "precursor_charge",
             "precursor_mz"]
        ]

        self.psms.columns = ["idPSM", "idPeptide", "fdr",
                             "sampleName", "scanNum", 
                             "pCharge", "pMz"]

    def _get_psm_info(self):
        psm_info = self.build.psms[
            np.logical_and(
                 self.build.psms.qvalue <= 0.01,
                 self.build.psms.label == "target"
            ) 
        ]

        psm_info = pd.DataFrame(
            product(
                np.unique(psm_info.reset_index().pep_id),
                np.arange(2, 7)
            ), columns = ["pep_id", "precursor_charge"]
        ).join(
            psm_info[
                ~psm_info.precursor_charge.isnull()
            ].groupby(["pep_id", "precursor_charge"]).id.count(),
            on = ["pep_id", "precursor_charge"],
            how = "left"
        ).fillna(0.)
        
        psm_info["precursor_charge"] = [
            "z{}".format(int(z)) for z in psm_info["precursor_charge"]
        ]
        psm_info = psm_info.pivot(index="pep_id", columns="precursor_charge")
        psm_info.columns = psm_info.columns.get_level_values(1)

        return psm_info
         
    def _convert_peptides(self):
        self.peptides = self.build.peptides[
            np.logical_and(
                self.build.peptides.qvalue <= 0.01,
                self.build.peptides.label == "target"
            )
        ].loc[:,
            ["id", "qvalue", "sequence", 
             "learned_rt", "hits", "test_error"]
        ]
        self.peptides.columns = ["idPeptide", "fdr", "sequence",
                                 "iRT", "nRTExamples", "errorRT"]

        self.peptides = self.peptides.join(self._get_psm_info(), on="idPeptide")

    def _read_in_fasta(self):
        # Read in proteins
        entries = []
        cur_entry = {}
        with open(self.fasta_path, "r") as src:
            for line in src:
                if line.startswith(">"):
                    if cur_entry.get("sequence", False):
                        entries.append(cur_entry)
                        cur_entry = {}
                    # Remove sp| and make accession/reference
                    uniprot = re.search("(?<=\>).+?(?=\s)", line).group()
                    uniprot = uniprot.lstrip("sp|")
                    uniprot = uniprot.split("|")
                    cur_entry["accession"] = uniprot[0]
                    cur_entry["reference"] = uniprot[-1]

                    # Pull out description of protein
                    desc = re.search("(?<=\s).+$", line)
                    cur_entry["description"] = desc if desc is None else desc.group()
                else:
                    cur_entry["sequence"] = cur_entry.get("sequence", "") + line.rstrip()
            entries.append(cur_entry)

        return pd.DataFrame.from_records(entries)

    def _convert_proteins(self):
        self.proteins = pd.DataFrame(
                {"prot_id" : self.peptides.idPeptide\
                                          .to_frame()\
                                          .join(self.build.peptide_protein\
                                                          .set_index("pep_id"),
                                                on="idPeptide"
                                          ).prot_id\
                                          .unique()
                }).join(
                    self.build.proteins[
                        self.build.proteins.label == "target"
                    ].drop("label", axis=1)\
                     .set_index("id"),
                    on="prot_id",
                    how="inner"
                )
        self.proteins.accession = self.proteins.accession.str.replace(".*\|", "")
        self.proteins = self.proteins.join(
                            self._read_in_fasta().set_index("accession"), 
                            on="accession"
                        )
        self.proteins.columns = ["idProtein", "accession", 
                                 "reference", "description", 
                                 "sequence"]

    def _convert_sites(self):
        # Build table of FDR corrected sites
        site_info = self.build.proteins.join(
            self.build.sites.reset_index().set_index("prot_id"),
            on = "id"
        ).set_index("id")

        site_info = site_info[
            np.logical_and(
                site_info.qvalue < 0.01,
                site_info.label == "target"
            )
        ]
        self.sites = self.proteins.idProtein.to_frame().join(
            site_info,
            on="idProtein",
            how="right"
        ).loc[:, ["index", "idProtein", "position", "residue", "qvalue"]]
        self.sites = self.sites.rename({"index" : "idSite"}, axis=1)
        self.sites.position += 1
        self.sites.idSite += 1

        # Build table linking peptide evidence to sites
        self.peptide_sites = self.peptides.loc[:,["idPeptide", "sequence"]].join(
            self.build.peptide_protein.set_index("pep_id"),
            on="idPeptide"
            ).join(
                self.proteins.set_index("idProtein"),
                on="prot_id",
                how="inner",
                lsuffix="_peptide",
                rsuffix="_protein"
            )
        self.peptide_sites["protein_pos"] = self.peptide_sites.apply(
            lambda row: row.sequence_protein.find(
                re.sub("[^A-Z]", "", row.sequence_peptide
            )),
            axis = 1
        )
        self.peptide_sites = self.peptide_sites[self.peptide_sites.protein_pos > -1]
        self.peptide_sites = self.peptide_sites.join(
            self.build.peptide_modifications.loc[
                self.build.peptide_modifications.mass == 80.0,
                ["pep_id", "pos", "score"]
            ].set_index("pep_id"),
            on = "idPeptide"
        )
        self.peptide_sites["localized_pos"] = (
            self.peptide_sites.protein_pos + self.peptide_sites.pos
            )
        self.peptide_sites = self.peptide_sites.join(
            self.sites.set_index(["idProtein", "position"]),
            on = ["prot_id", "localized_pos"],
            how = "inner"
            )
        self.peptide_sites = self.peptide_sites.loc[
            :, ["idPeptide", "idSite", "score"]
        ]
        self.peptide_sites = self.peptide_sites.rename({"score" :"ascore"}, axis=1)
        self.peptide_sites.ascore.values[np.isinf(self.peptide_sites.ascore)] = 1000.

        # Grade sites
        cutoff = 13
        select = self.peptide_sites.ascore > cutoff
        self.sites = self.sites.join(
            self.peptide_sites[select]\
                .groupby("idSite")\
                .idPeptide\
                .count()\
                .rename("grade"), 
            on = "idSite"
            )
        letter_grade = np.zeros(self.sites.shape[0], dtype=str)
        letter_grade[self.sites.grade > 1] = "A"
        letter_grade[self.sites.grade == 1] = "B"
        letter_grade[self.sites.grade.isna()] = "C"
        self.sites.grade = letter_grade

    def _annotate_pathways(self):
        annotations = pd.read_csv(self.annotation_path).loc[
            :, ["signature", "site.uniprot"]
        ]
     
        annotations.columns = ["name", "accession"]
        annotations["position"] = annotations.accession.str\
                                             .replace(".*;[A-Z]", "")\
                                             .astype(np.int)
        annotations["accession"] = annotations.accession.str.replace(";.*", "")

        self.pathways = pd.DataFrame(
            {"name" : annotations.name\
                                 .sort_values()\
                                 .unique()}
            ).reset_index()\
             .rename({"index" : "idPathway"},
                     axis=1)
     
        self.pathway_sites = self.pathways.join(
                annotations.set_index("name"),
                on = "name"
            ).join(
                self.proteins.set_index("accession")\
                             .idProtein,
                on = "accession",
                how = "inner"
           ).join(
                self.sites.set_index(["idProtein", "position"])\
                          .idSite,
                on = ["idProtein", "position"],
                how="inner"
            ).sort_values(
                ["idPathway", "idSite"]
            ).loc[:, ["idPathway", "idSite"]]

    def convert(self):
        self._convert_psms()
        self._convert_peptides()
        self._convert_proteins()
        self._convert_sites()

        # Remove sequence from proteins
        self.proteins = self.proteins.drop("sequence", axis=1)
        if self.annotation_path:
            self._annotate_pathways()

    def upload(self):
        table_names = ["psms", "peptides", "peptide_sites",
                       "proteins", "sites", "pathways", 
                       "pathway_sites"]
        table_objects = [PSM, Peptide, PeptideSite,
                         Protein, Site, Pathway,
                         PathwaySite]

        session = scoped_session(sessionmaker(bind=self.engine))
        for name, obj in zip(table_names, table_objects):
            table = getattr(self, name)
            if table is not None:
                print("Uploading: " + name)
                session.query(obj).delete()
                session.commit()
                for begin in range(0, table.shape[0], self.buffer_size):
                    end = min(begin + self.buffer_size, table.shape[0])
                    sub_table = table.iloc[begin:end, :]
                    session.execute(obj.__table__.insert(),
                                    sub_table.to_dict("records"))
                    session.commit()
                    print("Uploaded {} of {}".format(end, table.shape[0]))

        session.remove()

if __name__ == "__main__":
    host = "sqlite:///phosphopedia.db"
    base_path = "/net/villen/vol2/users/valenta4/Phosphopedia/phosphopedia_yeast/"
    build_dir = base_path + "integration_engine"
    fasta_path = base_path + "config/UP000002311_saccharomyces_cerevisiae_2020_03_22.fasta"
    annotation_path = None

    uploader = Uploader(host, build_dir, fasta_path) 
    uploader.convert()
    uploader.upload()
