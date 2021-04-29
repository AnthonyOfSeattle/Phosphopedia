import time
from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.exc import OperationalError
from .schema import *


class DatabaseBackend:
    def __init__(self, database_path, backoffs=100):
        self.database_path = database_path
        self.backoffs = 100

    def _backoff(self):
        yield True
        for attempt in range(self.backoffs):
            time.sleep(0.1 * 2 ** attempt)
            yield True

        yield False

    def create_engine(self):
        return create_engine(self.database_path, poolclass=NullPool)

    def create_session(self):
        return scoped_session(sessionmaker(bind=self.create_engine()))

    def initialize_database(self):
        with self.create_engine().connect() as connection:
            try:
                PhosphopediaBase.metadata.create_all(
                    connection
                )
            except OperationalError:
                print("Database locked, ignoring initialization check.")

    def safe_add(self, obj):
        session = self.create_session()
        while self._backoff():
            try:
                session.add(obj)
                session.commit()
                session.remove()
                return True
            except OperationalError:
                session.rollback()

        session.remove()
        return False

    def safe_run(self, func):
        while self._backoff():
            try:
                return func()
            except OperationalError:
                continue

        OperationalError()
