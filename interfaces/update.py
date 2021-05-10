import os
import re
import ppx
import glob
from .backend import DatabaseBackend
from .schema import Dataset, Sample


class DatasetManager:
    def __init__(self, database_path):
        self.database = DatabaseBackend(database_path)

    def _add_local_dataset(self, path, filter_str):
        # Check if dataset is already in database
        title = os.path.split(path)[-1]
        dataset_query = self.database.session.query(Dataset).filter(Dataset.title == title)
        if self.database.safe_run(dataset_query.count) == 0:
            print("==> {} not in database, adding now".format(title))
            previous_data_query = self.database.session\
                                               .query(Dataset)\
                                               .filter(Dataset.accession.like("LOC%"))
            n_previous = self.database.safe_run(dataset_query.count)

            index = str(n_previous + 1)
            accession = "LOC" + "0"*(6 - len(index)) + index
            dataset = Dataset(accession=accession, 
                              title=title)

        else: 
            print("==> {} already present in database, updating files".format(title))
            dataset = self.database.safe_run(dataset_query.one)

        # Check for files in path
        file_list = [os.path.split(f)[-1] for f in glob.glob(os.path.join(path, "*.raw"))]
        file_list = [f for f in file_list if re.search(filter_str, f) is not None]
        
        # Add files if not present
        for file_name in file_list:
            sample = Sample(parentDataset = dataset.accession,
                            sampleName = os.path.splitext(file_name)[0],
                            fileType = os.path.splitext(file_name)[1].lstrip(".").lower(),
                            fileName = file_name,
                            fileLocation = path
                           )
            sample_query = self.database.session\
                                        .query(Sample)\
                                        .filter(Sample.parentDataset == sample.parentDataset)\
                                        .filter(Sample.sampleName == sample.sampleName)
            if self.database.safe_run(sample_query.count) == 0:
                dataset.samples.append(sample)

        # Update database
        self.database.safe_add(dataset)

    def _add_remote_dataset(self, accession, filter_str):
        # Initiate ppx call for accesion validation
        cache_path = os.makedirs(os.path.join(os.getcwd(), ".ppx_cache"),
                                 exist_ok=True)
        project = ppx.find_project(accession, local=cache_path)

        # Check if dataset is already in database
        dataset_query = self.database\
                            .session\
                            .query(Dataset)\
                            .filter(Dataset.accession == accession)
        if self.database.safe_run(dataset_query.count) == 0:
            print("==> {} not in database, adding now".format(accession))
            dataset = Dataset(accession=accession,
                              title=accession)

        else:
            print("==> {} already present in database, updating files".format(accession))
            dataset = self.database.safe_run(dataset_query.one)

        # Check for files in remote
        file_list = project.remote_files("*.raw")
        file_list = [f for f in file_list if re.search(filter_str, f) is not None]

        # Add files if not present
        for file_name in file_list:
            sample = Sample(parentDataset = dataset.accession,
                            sampleName = os.path.splitext(file_name)[0],
                            fileType = os.path.splitext(file_name)[1].lstrip(".").lower(),
                            fileName = file_name,
                            fileLocation = "remote"
                           )
            sample_query = self.database.session\
                                        .query(Sample)\
                                        .filter(Sample.parentDataset == sample.parentDataset)\
                                        .filter(Sample.sampleName == sample.sampleName)
            if self.database.safe_run(sample_query.count) == 0:
                dataset.samples.append(sample)

        # Update database
        self.database.safe_add(dataset)

    def add_datasets(self, datasets):
        assert isinstance(datasets, list)

        for entry in datasets:
            if not isinstance(entry, list):
                entry = [entry, "."]

            print("==> Attempting to add {} to database".format(entry[0]))
            if re.search("(^PXD)|(^MSV)", entry[0].upper()) is not None:
                self._add_remote_dataset(*entry)
            elif os.path.isdir(entry[0]):
                self._add_local_dataset(*entry)
            else:
                raise ValueError(
                    "==> Error: {} is not a directory or proteome exchange id".format(entry[0])
                    )

            print("==> Success: {} added to database".format(entry[0]))
