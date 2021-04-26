from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.exc import OperationalError
from .schema import *

SAFE_QUERY_ITERS = 1000

def safe_get_all(query):
    for i in map(time.sleep, np.random.uniform(0, 5, size=SAFE_QUERY_ITERS)):
        try:
            result = query.all()
            break
        except OperationalError:
            continue

    return result

def safe_get_one(query):
    for i in map(time.sleep, np.random.uniform(0, 5, size=SAFE_QUERY_ITERS)):
        try:
            result = query.one()
            break
        except OperationalError:
            continue

    return result

def safe_update(session, obj):
    for i in map(time.sleep, np.random.uniform(0, 5, size=SAFE_QUERY_ITERS)):
        try:
            session.add(obj)
            session.commit()
            break
        except OperationalError:
            session.rollback()

class DatabaseBackend:
    def __init__(self, database_path):
        self.database_path = database_path

    def _get_engine(self):
        return create_engine(self.database_path, poolclass=NullPool)

    def _get_session(self):
        Session = scoped_session(sessionmaker(bind=self._get_engine()))
        return Session

    def initialize_database(self):
        with self._get_engine().connect() as connection:
            try:
                PhosphopediaBase.metadata.create_all(
                    connection
                )
            except OperationalError:
                print("Database locked, ignoring initialization check.")
