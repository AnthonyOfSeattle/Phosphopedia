import time
from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.exc import OperationalError
from .schema import *


class DatabaseBackend:
    def __init__(self, database_path, backoffs=100):
        self.database_path = database_path
        self.engine = create_engine(self.database_path, poolclass=NullPool)
        self.session = scoped_session(sessionmaker(bind=self.engine))
        self.backoffs = 100

    def __del__(self):
        self.session.remove()

    def _backoff(self):
        yield True
        for attempt in range(self.backoffs):
            time.sleep(0.1 * 2 ** attempt)
            yield True

        yield False

    def initialize_database(self):
        with self.engine.connect() as connection:
            try:
                PhosphopediaBase.metadata.create_all(
                    connection
                )
            except OperationalError:
                print("Database locked, ignoring initialization check.")

    def safe_add(self, obj):
        while self._backoff():
            try:
                self.session.add(obj)
                self.session.commit()
                self.session.remove()
                return True
            except OperationalError:
                self.session.rollback()

        return False

    def safe_run(self, func):
        while self._backoff():
            try:
                return func()
            except OperationalError:
                continue

        OperationalError()
