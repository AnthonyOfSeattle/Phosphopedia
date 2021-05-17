import unittest
import tempfile
import os
from interfaces import schema, managers


class TestSampleManager(unittest.TestCase):
    def test_error(self):
        """Add error code to files"""

        with tempfile.TemporaryDirectory() as temp_path:
            test_db_path = "sqlite:///" + temp_path + "/phosphopedia.db"

            # Add fake raw files to database
            nfiles = 10
            for ind in range(nfiles):
                open(temp_path + f"/file_{ind}.raw", "w").close()

            print()
            dataset_manager = managers.DatasetManager(test_db_path)
            dataset_manager.add_datasets([temp_path])

            sample_manager = managers.SampleManager(test_db_path)
            error_file_id = sample_manager.lookup_id("LOC000001", "file_0")
            sample_manager.add_error(error_file_id, "testError")
            
            manifest = sample_manager.get_sample_manifest()
            self.assertEqual(manifest.shape[0], nfiles - 1)

    def test_flags(self):
        """Add flags to files"""

        with tempfile.TemporaryDirectory() as temp_path:
            test_db_path = "sqlite:///" + temp_path + "/phosphopedia.db"

            # Add fake raw files to database
            nfiles = 10
            for ind in range(nfiles):
                open(temp_path + f"/file_{ind}.raw", "w").close()

            print()
            dataset_manager = managers.DatasetManager(test_db_path)
            dataset_manager.add_datasets([temp_path])

            sample_manager = managers.SampleManager(test_db_path)
            flag_file_id = sample_manager.lookup_id("LOC000001", "file_1")
            sample_manager.add_flag(flag_file_id, "testFlag")

            manifest = sample_manager.get_incomplete_samples("testFlag")
            self.assertEqual(manifest.shape[0], nfiles - 1)
