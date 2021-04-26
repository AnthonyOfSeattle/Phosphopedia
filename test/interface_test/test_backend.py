import unittest
import os
import tempfile
from interfaces import backend, schema

class TestDatabaseInit(unittest.TestCase):
    def test_init(self):
        """Create database and check that all required tables are present"""

        # Database must be on disk due to NullPool for engine creation
        with tempfile.TemporaryDirectory() as temp_path:
            test_db_path = "sqlite:///" + temp_path + "/phosphopedia.db"
            database = backend.DatabaseBackend(test_db_path)
            database.initialize_database()

            engine = database._get_engine()
            for name, cls in schema.__dict__.items():
                if isinstance(cls, type) and hasattr(cls, "__tablename__"):
                    self.assertTrue(
                        engine.dialect.has_table(engine, cls.__tablename__)
                        )

