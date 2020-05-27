import os
import re
import time
import glob
import urllib3
import json
from sqlalchemy.exc import OperationalError
from sqlalchemy.pool import NullPool
from sqlalchemy.orm import scoped_session
from database_schema import *

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

class DatabaseInterface:
    def __init__(self, database_path):
        self.database_path = database_path

    def __get_engine(self):
        return create_engine(self.database_path, poolclass=NullPool)

    def __get_session(self):
        Session = scoped_session(sessionmaker(bind=self.__get_engine()))
        return Session

    def __get_pride_dataset(self, dataset, include_string):
        http = urllib3.PoolManager()
        session = self.__get_session() 

        if session.query(Dataset).filter(Dataset.accession == dataset).count():
            print("==> {} already present in database, so it was ignored".format(dataset))
            session.remove()
            return

        # Issue GET request to PRIDE for dataset meta
        url = 'http://www.ebi.ac.uk/pride/ws/archive/project/' + dataset
        req = http.request('GET', url)
        dataset_meta = json.loads(req.data.decode("utf8"))
        dataset_entry = Dataset(
            accession=dataset_meta["accession"],
            title=dataset_meta["title"]
        )
        print("==> {} found, gathering samples".format(dataset))

        # Issue GET request to PRIDE for samples
        url = 'http://www.ebi.ac.uk/pride/ws/archive/file/list/project/' + dataset
        req = http.request('GET', url)
        sample_list = [sample for sample in json.loads(req.data.decode("utf8"))["list"] if sample["fileType"] == "RAW"]

        for sample in sample_list:
            if re.search(include_string, sample["fileName"]) is not None:
               dataset_entry.samples.append(
                   Sample(
                       parentDataset = dataset_entry.accession,
                       sampleName = os.path.splitext(sample["fileName"])[0],
                       fileName = sample["fileName"],
                       fileSize = sample["fileSize"],
                       fileLocation = sample["downloadLink"]
                   )
               )

        try:
            session.add(dataset_entry)
            session.commit()
            print("==> SUCCESS: {} added to database.".format(dataset))
        except DatabaseError:
            session.rollback()
            print("==> DatabaseError caused {} not to be added".format(dataset)) 
        finally:
            session.remove()

    def __get_local_dataset(self, dataset_path):
        # Check if data is in database by name
        session = self.__get_session()
        dataset_title = os.path.split(dataset_path)[1]
        if session.query(Dataset).filter(Dataset.title == dataset_title).count():
           print("==> {} already present in database, so it was ignored".format(dataset_title))
           session.remove()
           return

        n_previous_local = session.query(Dataset).filter(Dataset.accession.like("LOC%")).count()
        local_index = str(n_previous_local + 1)
        dataset_accession = "LOC" + "0"*(6 - len(local_index)) + local_index

        dataset_entry = Dataset(
           accession=dataset_accession,
           title=dataset_title
        )
        print("==> {} found, gathering samples".format(dataset_title))

        sample_list = glob.glob(
            os.path.join(dataset_path, "*.raw")
        )

        for sample in sample_list:
            dataset_entry.samples.append(
                Sample(
                    parentDataset = dataset_entry.accession,
                    sampleName = os.path.splitext(os.path.split(sample)[1])[0],
                    fileName = os.path.split(sample)[1],
                    fileSize = os.stat(sample).st_size,
                    fileLocation = sample
                )
            )

        try:
            session.add(dataset_entry)
            session.commit()
            print("==> SUCCESS: {} added to database.".format(dataset_title))
        except OperationalError:
            session.rollback()
            print("==> DatabaseError caused {} not to be added".format(dataset_title))
            raise
        finally:
            session.remove()
         
    def initialize_database(self):
        with self.__get_engine().connect() as connection:
            try:
                PhosphopediaBase.metadata.create_all(
                    connection
                )
            except OperationalError:
                print("Database locked, ignoring initialization check.")

    def update_datasets(self, dataset_list):
        for dataset in dataset_list:
            if isinstance(dataset, list):
                include_string = dataset[1]
                dataset = dataset[0]

            print("==> Attempting to add {} to database".format(dataset))

            if dataset.upper().startswith("PXD"):
                self.__get_pride_dataset(dataset, include_string)
            elif os.path.isdir(dataset):
                self.__get_local_dataset(dataset)
            else:
                raise ValueError("{} is not a directory or proteome exchange id".format(dataset))

    def get_sample_manifest(self):
        session = self.__get_session()
        query = session.query(Sample.id,
                                           Sample.sampleName,
                                           Sample.parentDataset).\
                                     outerjoin(Error, Sample.id == Error.sampleId).\
                                     filter(Error.errorCode == None)

        q = safe_get_all(query)
        session.remove()
        return pd.DataFrame(
            q, columns = ["id", "sampleName", "parentDataset"]
        )

    def get_sample(self, sample_id = None, sample_name = None):
        session = self.__get_session()
        query = session.query(Sample)

        if sample_id is not None:
            query = query.filter(Sample.id == sample_id)

        elif sample_name is not None:
            query = query.filter(Sample.sampleName == sample_name)

        else:
            raise ValueError("Input sample id or name")

        q = safe_get_one(query)
        session.remove()
        return q

    def get_sample_params(self, sample_id = None, sample_name = None):
        session = self.__get_session()
        if sample_id is None and sample_name is not None:
            query = session.query(Sample.id).filter(Sample.sampleName == sample_name)
            sample_id = query.one()[0]

        elif sample_id is None and sample_name is None:
            raise ValueError("Input sample id or name")

        params = session.query(Parameters).filter(Parameters.sampleId == sample_id).one()
        session.remove()
        return params

    def __safe_get_one(self, query):
        for i in map(time.sleep, np.random.uniform(0, 5, size=SAFE_QUERY_ITERS)):
            try:
                result = query.one()
                break
            except OperationalError:
                continue

        return result

    def __safe_update(self, obj):
        session = self.__get_session()
        for i in map(time.sleep, np.random.uniform(0, 5, size=SAFE_QUERY_ITERS)):
            try:
                session.add(obj)
                session.commit()
                break
            except OperationalError:
                session.rollback()

        session.remove()

    def add_sample_params(self, sample_id = None, sample_name = None, 
                          ms1_analyzer = None, ms2_analyzer = None):
        session = self.__get_session()
        if sample_id is None and sample_name is not None:
            query = session.query(Sample.id).filter(Sample.sampleName == sample_name)
            sample_id = self.__safe_get_one(query)[0]

        elif sample_id is None and sample_name is None:
            raise ValueError("Input sample id or name")

        try:
            params = self.__safe_get_one(session.query(Parameters).filter(Parameters.sampleId == sample_id))
            params.ms1Analyzer = ms1_analyzer if ms1_analyzer is not None else params.ms1Analyzer
            params.ms2Analyzer = ms2_analyzer if ms2_analyzer is not None else params.ms2Analyzer
        except Exception as e:
            print(e)
            params = Parameters(sampleId = sample_id,
                                ms1Analyzer = ms1_analyzer,
                                ms2Analyzer = ms2_analyzer)

        safe_update(session, params)
        session.remove()
 

    def add_error(self, sample_id = None, sample_name = None, error_code = None):
        session = self.__get_session()
        if sample_id is None and sample_name is not None:
            query = session.query(Sample.id).filter(Sample.sampleName == sample_name)
            sample_id = self.__safe_get_one(query)[0]

        elif sample_id is None and sample_name is None:
            raise ValueError("Input sample id or name")

        error = Error(sampleId = sample_id,
                      errorCode = error_code)

        safe_update(session, error)
        session.remove()

