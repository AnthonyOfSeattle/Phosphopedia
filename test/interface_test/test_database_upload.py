import unittest
import os
import tempfile
import numpy as np
from interfaces import schema, DatabaseBuild, BuildUploader


class TestBuildUpload(unittest.TestCase):
    def test_filtering(self):
        """Test loading and filtering of integration engine database build"""

        sample_build_path = os.path.join(os.getcwd(),
                                         "test",
                                         "interface_test",
                                         "sample_build")

        fdr_cutoff = 0.1
        db_build = DatabaseBuild(sample_build_path, fdr_cutoff)

        # Check main table filters
        self.assertTrue(np.all(db_build.psms.qvalue <= fdr_cutoff))
        self.assertTrue(np.all(db_build.psms.label == "target"))
        self.assertTrue(np.all(db_build.peptides.qvalue <= fdr_cutoff))
        self.assertTrue(np.all(db_build.peptides.label == "target"))
        self.assertTrue(np.all(db_build.sites.qvalue <= fdr_cutoff))
        self.assertTrue(np.all(db_build.proteins.label == "target"))

        # Check table connections
        self.assertFalse(np.union1d(db_build.peptides.id,
                                    db_build.peptide_modifications.pep_id
                                   ).shape[0] > db_build.peptides.shape[0])
        self.assertFalse(np.union1d(db_build.proteins.id,
                                    db_build.sites.prot_id
                                   ).shape[0] > db_build.proteins.shape[0])

    def test_conversion_and_upload(self):
        """Test that all data is converted to its new form correctly"""

        # Database must be on disk due to NullPool for engine creation
        with tempfile.TemporaryDirectory() as temp_path:
            test_db_path = "sqlite:///" + temp_path + "/phosphopedia.db"
            sample_build_path = os.path.join(
                                    os.getcwd(),
                                    "test",
                                    "interface_test",
                                    "sample_build"
                                )
            sample_fasta_path = os.path.join(
                                    sample_build_path,
                                    "UP000002311_saccharomyces_cerevisiae_2020_03_22.fasta"
                                )
            sample_annotation_path = os.path.join(
                                         sample_build_path,
                                         "sample_pathways.csv"
                                     )
            database = BuildUploader(test_db_path,
                                     sample_build_path,
                                     sample_fasta_path,
                                     annotation_path = sample_annotation_path)
            database.convert()
            database.upload()            

            # Check psm statistics are correct
            for charge in range(2, 7):
                select = database.psms.pCharge == charge
                expected_detections = database.peptides\
                                              .join(database.psms[select]\
                                                            .set_index("idPeptide")\
                                                            .idPSM,
                                                    how="inner",
                                                    on="idPeptide")
                self.assertEqual(int(database.peptides["z" + str(charge)].sum()),
                                 expected_detections.shape[0])

            # Check that main filtering worked
            self.assertEqual(np.intersect1d(database.sites.idProtein,
                                            database.proteins.idProtein).shape[0],
                             np.union1d(database.sites.idProtein,
                                        database.proteins.idProtein).shape[0])

            self.assertEqual(np.intersect1d(database.sites.idSite,
                                            database.peptide_sites.idSite).shape[0],
                             np.union1d(database.sites.idSite,
                                        database.peptide_sites.idSite).shape[0])

            self.assertEqual(np.intersect1d(database.peptide_sites.idPeptide,
                                            database.peptides.idPeptide).shape[0],
                             np.union1d(database.peptide_sites.idPeptide,
                                        database.peptides.idPeptide).shape[0])

            self.assertEqual(np.intersect1d(database.peptides.idPeptide,
                                            database.psms.idPeptide).shape[0],
                             np.union1d(database.peptides.idPeptide,
                                        database.psms.idPeptide).shape[0])

            # Check that annotations work
            self.assertEqual(database.pathways.name.str.contains("POSITIVE").sum(), 3)
            self.assertEqual(database.pathways.name.str.contains("NEGATIVE").sum(), 0)
            self.assertEqual(database.pathway_sites.shape[0], 9)

