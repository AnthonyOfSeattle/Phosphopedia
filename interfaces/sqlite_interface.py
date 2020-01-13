import os
import glob
import re
import sqlite3
import urllib3
import json
import pandas as pd

class SQLiteInterface:
    """
    Performs common pipeline queries with easy to use interface.
    """
    def __init__(self, database_name):
        self.connection = sqlite3.connect(database_name)
        self.cursor = self.connection.cursor()
        self.http = urllib3.PoolManager()

    def __del__(self):
        self.connection.close()

    def __init_dataset_table(self):
        table_exists = self.cursor.execute(
            """
            SELECT count(name) FROM sqlite_master
            WHERE type='table' and name='dataset';
            """
        ).fetchone()[0]

        if not table_exists:
            self.cursor.execute(
                """
                CREATE TABLE dataset (
                    accession TEXT NOT NULL PRIMARY KEY,
                    title TEXT NOT NULL
                );
                """
            )

        self.connection.commit()

    def __init_sample_table(self):
        table_exists = self.cursor.execute(
            """
            SELECT count(name) FROM sqlite_master
            WHERE type='table' and name='sample';
            """
        ).fetchone()[0]

        if not table_exists:
            self.cursor.execute(
                """
                CREATE TABLE sample (
                    sampleID INTEGER PRIMARY KEY,
                    sampleName TEXT NOT NULL,
                    fileName TEXT NOT NULL,
                    fileSize INTEGER NOT NULL,
                    fileLocation TEXT NOT NULL,
                    ms1Analyzer TEXT,
                    ms2Analyzer TEXT,
                    errorCode TEXT,
                    parentDataset TEXT NOT NULL,
                    FOREIGN KEY(parentDataset) REFERENCES dataset(accession)
                );
                """
            )

        self.connection.commit()

    def __get_pride_dataset_meta(self, dataset):
        # Check if data is in database by accession
        dataset_exists = self.cursor.execute(
            """
            SELECT count(accession) FROM dataset
            WHERE accession='{}';
            """.format(dataset)
        ).fetchone()[0]

        if dataset_exists:
            raise sqlite3.Error

        # Issue GET request to PRIDE
        url = 'http://www.ebi.ac.uk/pride/ws/archive/project/' + dataset
        req = self.http.request('GET', url)
        dataset_meta = json.loads(req.data.decode("utf8"))

        self.cursor.execute(
            """
            INSERT INTO dataset (
                accession, title
            ) 
            VALUES (
                '{accession}', '{title}'
            );
            """.format(**dataset_meta)
        )

    def __get_local_dataset_meta(self, dataset):
        # Check if data is in database by name
        dataset_title = os.path.split(dataset)[1]

        dataset_exists = self.cursor.execute(
            """
            SELECT count(accession) FROM dataset
            WHERE title='{}';
            """.format(dataset_title)
        ).fetchone()[0]

        if dataset_exists:
            raise sqlite3.Error

        n_previous_local = self.cursor.execute(
            """
            SELECT count(accession) FROM dataset
            WHERE accession LIKE 'LOC%';
            """
        ).fetchone()[0]

        local_index = str(n_previous_local + 1)
        dataset_accession = "LOC" + "0"*(6 - len(local_index)) + local_index 

        self.cursor.execute(
            """
            INSERT INTO dataset (
                accession, title
            ) 
            VALUES (
                '{accession}', '{title}'
            );
            """.format(accession=dataset_accession, title=dataset_title)
        )

    def __get_pride_dataset_samples(self, dataset, include_string="."):
        # Issue GET request to PRIDE
        url = 'http://www.ebi.ac.uk/pride/ws/archive/file/list/project/' + dataset
        req = self.http.request('GET', url)
        sample_list = json.loads(req.data.decode("utf8"))["list"]

        tuples = []
        for sample in sample_list:
            if sample["fileType"] == "RAW":
                if re.search(include_string, sample["fileName"]) is None:
                    continue

                sample["parentDataset"] = dataset
                sample["sampleName"] = os.path.splitext(sample["fileName"])[0]
                tuples.append((sample["sampleName"], sample["fileName"], 
                               sample["fileSize"], sample["downloadLink"],
                               sample["parentDataset"]))
                print(sample["sampleName"])

        statement = """INSERT INTO sample (
                       sampleName, fileName, fileSize, 
                       fileLocation, parentDataset
                       ) VALUES (?, ?, ?, ?, ?);""" 
        self.cursor.executemany(statement, tuples)

    def __get_local_dataset_samples(self, dataset):
        dataset_title = os.path.split(dataset)[1]
        accession = self.cursor.execute(
            """
            SELECT accession FROM dataset
            WHERE title='{}';
            """.format(dataset_title)
        ).fetchone()[0]

        file_list = glob.glob(
            os.path.join(dataset, "*.raw")
        )

        tuples = []
        for path in file_list:
            tuples.append((
                os.path.splitext(os.path.split(path)[1])[0],
                os.path.split(path)[1],
                os.stat(path).st_size,
                path,
                accession
            ))  

        statement = """INSERT INTO sample (
                       sampleName, fileName, fileSize, 
                       fileLocation, parentDataset
                       ) VALUES (?, ?, ?, ?, ?);"""
        self.cursor.executemany(statement, tuples)

    def update_datasets(self, dataset_list):
        self.__init_dataset_table()
        self.__init_sample_table()

        for dataset in dataset_list:
            if isinstance(dataset, list):
                include_string = dataset[1]
                dataset = dataset[0]

            print("==> Attempting to add {} to database".format(dataset))
            try: 
                if dataset.upper().startswith("PXD"):
                   self.__get_pride_dataset_meta(dataset)
                   self.__get_pride_dataset_samples(dataset, include_string)
                   self.connection.commit()
                elif os.path.isdir(dataset):
                   self.__get_local_dataset_meta(dataset)
                   self.__get_local_dataset_samples(dataset)
                   self.connection.commit()
                else:
                   print("    FAIL: {} is not a directory or proteome exchange id".format(dataset))
                   raise ValueError
                print("    SUCCESS: {} added to database.".format(dataset))
            except:
                print("    IGNORE: {} not added to database.".format(dataset))
                self.connection.rollback()

        return self

    def get_sample_manifest(self):
        return pd.read_sql_query("""SELECT sampleName, fileLocation, accession 
                                    FROM sample LEFT JOIN dataset ON parentDataset=accession;""", 
                                 self.connection)

    def update_samples(self, column_name, new_value, where_clause=''):
        statement = """
                    UPDATE Sample
                    SET {} = {}
                    {};
                    """.format(column_name, new_value, where_clause)
        self.cursor.execute(statement)
        self.connection.commit()

    def query_samples(self, column_names, where_clause=''):

        if isinstance(column_names, str):
            column_names = [column_names]

        statement = """
                    SELECT {}
                    FROM Sample
                    {}
                    """.format(",".join(column_names), where_clause)
        return pd.read_sql_query(statement, self.connection)
