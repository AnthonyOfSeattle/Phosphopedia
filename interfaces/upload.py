import re
import numpy as np
import pandas as pd
from itertools import product
from .backend import DatabaseBackend
from .containers import DatabaseBuild
from .database_schema import *

class BuildUploader:
    def __init__(self,
                 database_path,
                 build_path,
                 fasta_path,
                 annotation_path=None,
                 buffer_size=1e6):
        self.database = DatabaseBackend(database_path)
        self.database.initialize_database()
        self.buffer_size = int(buffer_size)

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
        psms = self.build.psms.loc[:, ["id", 
                                       "pep_id",
                                       "qvalue",
                                       "sample_name",
                                       "scan_number",
                                       "precursor_charge",
                                       "precursor_mz"]]

        psms.columns = ["idPSM", "idPeptide", "fdr",
                        "sampleName", "scanNum",
                        "pCharge", "pMz"]

        return psms

    def _get_charge_counts(self):
        charge_counts = self.psms[~self.psms.pCharge.isnull()]\
                            .groupby(["idPeptide", "pCharge"])\
                            .idPSM\
                            .count()\
                            .rename("nDetections")

        charge_counts = pd.DataFrame(product(self.psms.idPeptide.unique(),
                                             np.arange(2, 7)), 
                                     columns = ["idPeptide", "pCharge"])\
                                     .join(charge_counts,
                                           how="left",
                                           on=["idPeptide", "pCharge"])\
                                     .fillna(0)
        
        charge_counts["pCharge"] = [
            "z{}".format(int(z)) for z in charge_counts["pCharge"]
        ]

        charge_counts = charge_counts.pivot(index="idPeptide", columns="pCharge")
        charge_counts.columns = charge_counts.columns.get_level_values(1)

        return charge_counts

    def _convert_peptides(self):
        peptides = self.build.peptides.loc[:, ["id",
                                               "qvalue",
                                               "sequence",
                                               "learned_rt",
                                               "hits", 
                                               "test_error"]
                                          ]

        peptides.columns = ["idPeptide", "fdr", "sequence",
                            "iRT", "nRTExamples", "errorRT"]

        peptides = peptides.join(self._get_charge_counts(), on="idPeptide")
        
        return peptides

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
        proteins = self.build.proteins.drop("label", axis=1)
        proteins.accession = proteins.accession.str.replace(".*\|", "", regex=True)
        proteins = proteins.join(
                       self._read_in_fasta().set_index("accession"),
                       on="accession"
                       )
        proteins.columns = ["idProtein", "accession",
                            "reference", "description",
                            "sequence"]
        return proteins

    def _convert_sites(self):
        sites = self.build.sites.set_index("prot_id")
        sites["idSite"] = np.arange(sites.shape[0], dtype=int) + 1
        sites["position"] += 1
        sites = self.proteins\
                    .idProtein\
                    .to_frame()\
                    .join(sites,
                          on="idProtein",
                          how="right"
                          )\
                    .loc[:, ["idSite", "idProtein", 
                             "position", "residue",
                             "qvalue"]]\
                    .rename({"qvalue" : "fdr"}, axis=1)

        # To remove, explicit filter for STY
        sites = sites[sites.residue.str.contains("[STY]", regex=True)]
        return sites

    def _link_peptides_to_sites(self):
        # Map peptide to start in protein
        peptide_sites = self.peptides.loc[:,["idPeptide", "sequence"]]\
                                     .join(self.build.peptide_protein.set_index("pep_id"),
                                           on="idPeptide")\
                                     .join(self.proteins.set_index("idProtein").sequence,
                                           on="prot_id",
                                           how="inner",
                                           lsuffix="_peptide",
                                           rsuffix="_protein"
                                          )
        peptide_sites.sequence_peptide = peptide_sites.sequence_peptide\
                                                      .str\
                                                      .replace("[^A-Z]", "", regex=True)
        peptide_sites.sequence_peptide = peptide_sites.sequence_peptide\
                                                      .str\
                                                      .replace("[IL]", "[IL]", regex=True)
        peptide_sites["protein_pos"] = peptide_sites.apply(
            lambda row: re.search(row.sequence_peptide, row.sequence_protein).start(),
            axis = 1
        )

        # Load modifications on peptide and find true modification sites
        peptide_sites = peptide_sites.join(
            self.build.peptide_modifications.loc[
                self.build.peptide_modifications.mass == 80.0,
                ["pep_id", "pos", "score"]
            ].set_index("pep_id"),
            on = "idPeptide"
        )
        peptide_sites["localized_pos"] = (
            peptide_sites.protein_pos + peptide_sites.pos
            )

        # Match peptide sites to real sites
        # To Fix: Should be right join
        peptide_sites = peptide_sites.join(
            self.sites.set_index(["idProtein", "position"]),
            on = ["prot_id", "localized_pos"],
            how = "inner"
            )
        peptide_sites = peptide_sites.loc[
            :, ["idPeptide", "idSite", "score"]
        ]
        peptide_sites = peptide_sites.rename({"score" :"ascore"}, axis=1)
        peptide_sites.ascore.values[np.isinf(peptide_sites.ascore)] = 1000.

        return peptide_sites

    def _reduce_proteins(self):
        proteins = self.sites.idProtein\
                             .drop_duplicates()\
                             .to_frame()\
                             .join(self.proteins.set_index("idProtein"),
                                   on="idProtein")

        return proteins

    def _reduce_peptides(self):
        peptides = self.peptide_sites.idPeptide\
                                     .drop_duplicates()\
                                     .to_frame()\
                                     .join(self.peptides.set_index("idPeptide"),
                                           on="idPeptide")

        return  peptides

    def _reduce_psms(self):
        psms = self.peptides.idPeptide\
                            .to_frame()\
                            .join(self.psms.set_index("idPeptide"),
                                  on="idPeptide")

        return psms

    def _annotate_pathways(self):
        annotations = pd.read_csv(self.annotation_path).iloc[:, :3]
        annotations.columns = ["name", "accession", "position"]

        pathways = annotations.name\
                              .drop_duplicates()\
                              .reset_index()\
                              .rename({"index" : "idPathway"},
                                      axis=1)

        pathway_sites = pathways.join(annotations.set_index("name"),
                                      on="name")\
                                .join(self.proteins.set_index("accession"),
                                      how="inner",
                                      on="accession")\
                                .join(self.sites.set_index(["idProtein", "position"]).idSite,
                                      how="inner",
                                      on=["idProtein", "position"])\
                                .loc[:,["idPathway", "idSite"]]

        pathways = pathway_sites.idPathway\
                                .drop_duplicates()\
                                .to_frame()\
                                .join(pathways.set_index("idPathway"),
                                      on="idPathway")

        return pathways, pathway_sites    
                                 
    def convert(self):
        # Load individual tables
        self.psms = self._convert_psms()
        self.peptides = self._convert_peptides()
        self.proteins = self._convert_proteins()
        self.sites = self._convert_sites()
        self.peptide_sites = self._link_peptides_to_sites()

        # Reduce tables with sites as reference
        self.proteins = self._reduce_proteins()
        self.peptides = self._reduce_peptides()
        self.psms = self._reduce_psms()
    
        # Remove sequence from proteins
        self.proteins = self.proteins.drop("sequence", axis=1)

        # Annotate pathways only if provided
        if self.annotation_path is not None:
            self.pathways, self.pathway_sites = self._annotate_pathways()

    def upload(self):
        table_names = ["psms", "peptides", "peptide_sites",
                       "proteins", "sites", "pathways",
                       "pathway_sites"]
        table_objects = [PSM, Peptide, PeptideSite,
                         Protein, Site, Pathway,
                         PathwaySite]

        for name, obj in zip(table_names, table_objects):
            table = getattr(self, name)
            if table is None:
                continue

            print("Dropping from: " + name)
            drop_query = self.database.session.query(obj).delete
            self.database.safe_run(drop_query)
            print("Uploading: " + name)
            for begin in range(0, table.shape[0], self.buffer_size):
                end = min(begin + self.buffer_size, table.shape[0])
                sub_table = table.iloc[begin:end, :]
                self.database.safe_add_bulk(obj, sub_table.to_dict("records"))
                print("==> {} of {}".format(end, table.shape[0]))
