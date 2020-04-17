import os
import re
import glob
import urllib3
import json
from sqlalchemy.exc import OperationalError
from database_schema import *

class DatabaseInterface:
    def __init__(self, database_path):
        self.database_path = database_path

    def __get_engine(self):
        return create_engine(self.database_path)

    def __get_session(self):
        Session = sessionmaker(bind=self.__get_engine())
        return Session()

    def __get_pride_dataset(self, dataset, include_string):
        http = urllib3.PoolManager()
        session = self.__get_session() 

        if session.query(Dataset).filter(Dataset.accession == dataset).count():
            print("==> {} already present in database, so it was ignored".format(dataset))
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
                       accession = dataset_entry.accession,
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

    def __get_local_dataset(self, dataset_path):
        # Check if data is in database by name
        session = self.__get_session()
        dataset_title = os.path.split(dataset_path)[1]
        if session.query(Dataset).filter(Dataset.title == dataset_title).count():
           print("==> {} already present in database, so it was ignored".format(dataset_title))
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
                    accession = dataset_entry.accession,
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
         
    def initialize_database(self):
        try:
            PhosphopediaBase.metadata.create_all(
                self.__get_engine()
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
