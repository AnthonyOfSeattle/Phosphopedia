import unittest
import tempfile
import ppx
import os
from interfaces import schema, update


class TestDatasetManager(unittest.TestCase):
    def test_local_update(self):
        """Load in raw files saved in a local drive"""

        with tempfile.TemporaryDirectory() as temp_path:
            test_db_path = "sqlite:///" + temp_path + "/phosphopedia.db"

            # Add fake raw files to database
            nfiles = 100
            for ind in range(nfiles):
                open(temp_path + f"/file_{ind}.raw", "w").close()

            print()
            manager = update.DatasetManager(test_db_path)
            manager.add_datasets([temp_path])

            # Make sure all files added
            sample_query = manager.database\
                                  .session\
                                  .query(schema.Sample)\
                                  .order_by(schema.Sample.sampleName)
            sample_entries = manager.database.safe_run(sample_query.all)
            self.assertEqual(len(sample_entries), nfiles)

            # Add more fake raw files
            nfiles = 150
            for ind in range(nfiles):
                open(temp_path + f"/file_{ind}.raw", "w").close()

            print()
            manager.add_datasets([temp_path])

            # Make sure dataset only added once and new files added
            dataset_query = manager.database\
                                   .session\
                                   .query(schema.Dataset)\
                                   .filter(schema.Dataset.title == temp_path.split("/")[-1])
            dataset_entries = manager.database.safe_run(dataset_query.all)
            self.assertEqual(len(dataset_entries), 1)

            sample_query = manager.database\
                                  .session\
                                  .query(schema.Sample)\
                                  .order_by(schema.Sample.sampleName)
            sample_entries = manager.database.safe_run(sample_query.all)
            self.assertEqual(len(sample_entries), nfiles)

    def test_remote_update(self):
        """Test remote connection to ProteomeXchange"""

        prev_wd = os.getcwd()
        with tempfile.TemporaryDirectory() as temp_path:
            os.chdir(temp_path)
            test_db_path = "sqlite:///" + temp_path + "/phosphopedia.db"

            print()
            accession = "PXD001492"
            manager = update.DatasetManager(test_db_path)
            manager.add_datasets([accession])

            # Make sure all files added
            project = ppx.find_project(accession, local=temp_path)
            nfiles = len(project.remote_files("*.raw"))
            sample_query = manager.database\
                                  .session\
                                  .query(schema.Sample)\
                                  .order_by(schema.Sample.sampleName)
            sample_entries = manager.database.safe_run(sample_query.all)
            self.assertEqual(len(sample_entries), nfiles)

            os.chdir(prev_wd)
