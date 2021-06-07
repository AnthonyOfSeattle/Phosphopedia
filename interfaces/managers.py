import os
import re
import ppx
import glob
import pandas as pd
from sqlalchemy.orm.exc import NoResultFound
from .backend import DatabaseBackend
from .schema import Dataset, Sample, Error, Flag, Parameters


class DatasetManager:
    def __init__(self, database_path):
        self.database = DatabaseBackend(database_path)
        self.database.initialize_database()

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
        file_list.sort()
        
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
                              title=project.title)

        else:
            print("==> {} already present in database, updating files".format(accession))
            dataset = self.database.safe_run(dataset_query.one)

        # Check for files in remote
        file_list = project.remote_files("*.raw")
        file_list = [f for f in file_list if re.search(filter_str, f) is not None]
        file_list.sort()

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


class SampleManager:
    def __init__(self, database_path):
        self.database = DatabaseBackend(database_path)

    def lookup_id(self, parent_dataset, sample_name):
        query = self.database\
                    .session\
                    .query(Sample.id)\
                    .filter(Sample.parentDataset == parent_dataset)\
                    .filter(Sample.sampleName == sample_name)

        return self.database.safe_run(query.one)[0]

    def add_error(self, sample_id, error_code):
        error = Error(sampleId = sample_id,
                      errorCode = error_code)

        self.database.safe_add(error)

    def add_flag(self, sample_id, flag_code):
        flag = Flag(sampleId = sample_id,
                    flagCode = flag_code)

        self.database.safe_add(flag)

    def get_parameters(self, sample_id):
        query = self.database\
                    .session\
                    .query(Parameters)\
                    .filter(Parameters.sampleId == sample_id)

        return self.database.safe_run(query.one)

    def add_mass_analyzers(self, sample_id, ms1_analyzer, ms2_analyzer):
        try:
            params = self.get_parameters(sample_id)
            params.ms1Analyzer = ms1_analyzer
            params.ms2Analyzer = ms2_analyzer
        except NoResultFound as e:
            params = Parameters(sampleId = sample_id,
                                ms1Analyzer = ms1_analyzer,
                                ms2Analyzer = ms2_analyzer)

        self.database.safe_add(params)

    def get_sample(self, sample_id):
        query = self.database\
                    .session\
                    .query(Sample)\
                    .filter(Sample.id == sample_id)

        return self.database.safe_run(query.one)

    def get_sample_manifest(self, exclude_errors=True):
        query = self.database\
                    .session\
                    .query(Sample.id,
                           Sample.sampleName,
                           Sample.parentDataset)

        if exclude_errors:
            query = query.outerjoin(Error, Sample.id == Error.sampleId)\
                         .filter(Error.errorCode == None)

        sample_manifest = pd.DataFrame(
            self.database.safe_run(query.all),
            columns = ["id", "sampleName", "parentDataset"]
        )

        return sample_manifest
        
    def get_incomplete_samples(self, flag_code):
        flag_query = self.database\
                         .session\
                         .query(Flag.sampleId, Flag.flagCode)\
                         .filter(Flag.flagCode == flag_code)

        flag_df = pd.DataFrame(
            self.database.safe_run(flag_query.all),
            columns = ["sampleId", "flagCode"]
        )

        sample_manifest = self.get_sample_manifest()
        sample_manifest = sample_manifest.join(flag_df.set_index("sampleId"),
                                               how="left",
                                               on="id")
        sample_manifest = sample_manifest[sample_manifest.flagCode.isna()]
        sample_manifest = sample_manifest.drop("flagCode", axis=1)        

        return sample_manifest
