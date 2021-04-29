import unittest
import os
import tempfile
from multiprocessing.pool import Pool
from interfaces import backend, schema


def new_dataset(ind, database):
    """Utility function to add single fake dataset to database"""

    str_ind = str(ind)
    accession = "LOC" + "0"*(6-len(str_ind)) + str_ind
    test_entry = schema.Dataset(accession=accession,
                                title="Test " + str_ind)
    database.safe_add(test_entry)


class TestDatabaseInit(unittest.TestCase):
    def test_init(self):
        """Create database and check that all required tables are present"""

        # Database must be on disk due to NullPool for engine creation
        with tempfile.TemporaryDirectory() as temp_path:
            test_db_path = "sqlite:///" + temp_path + "/phosphopedia.db"
            database = backend.DatabaseBackend(test_db_path)
            database.initialize_database()

            engine = database.create_engine()
            for name, cls in schema.__dict__.items():
                if isinstance(cls, type) and hasattr(cls, "__tablename__"):
                    self.assertTrue(
                        engine.dialect.has_table(engine, cls.__tablename__)
                        )

    def test_update(self):
        """Test database updates in a single process"""

        with tempfile.TemporaryDirectory() as temp_path:
            test_db_path = "sqlite:///" + temp_path + "/phosphopedia.db"
            database = backend.DatabaseBackend(test_db_path)
            database.initialize_database()

            new_dataset(1, database)

            # Get entry back
            session = database.create_session()
            query = session.query(schema.Dataset)
            retrieved_entry = database.safe_run(query.one)
            self.assertEqual(int(retrieved_entry.accession[3:]), 1)

            session.remove()

    def test_update_parallel(self):
        """Test database updates over multiple processes"""

        with tempfile.TemporaryDirectory() as temp_path:
            test_db_path = "sqlite:///" + temp_path + "/phosphopedia.db"
            database = backend.DatabaseBackend(test_db_path)
            database.initialize_database()

            with Pool() as pool:
                pool.starmap(new_dataset, [(ind, database) for ind in range(100)])

            # Get all entries back
            session = database.create_session()
            query = session.query(schema.Dataset).order_by(schema.Dataset.accession)
            retrieved_entries = database.safe_run(query.all)
            for ind, entry in enumerate(retrieved_entries):
                self.assertEqual(int(entry.accession[3:]), ind)

            session.remove()
