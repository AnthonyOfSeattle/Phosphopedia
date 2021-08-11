import os
from sqlalchemy.exc import OperationalError
from .backend import DatabaseBackend
from .schema import Sample, Error


class ErrorCleaner:
    def __init__(self, working_dir, database_path):
        self.active = False
        self.working_dir = working_dir
        self.database = DatabaseBackend(database_path)
        if self._is_error() and self._clean_requested():
            self.active = True

    def _is_error(self):
        error_query = self.database.session.query(Error)
        try:
            error_count = error_query.count()
        except OperationalError:
            error_count = 0

        return error_count > 0

    def _clean_requested(self):
        response = ""
        while response.lower() not in ["y", "n"]:
            response = input("Previously errored files detected, "
                             "would you like to clean and rerun? (y/n) ")

        return response == "y"

    def _del_sample_path(self, dataset, sample):
        print(f"==> Clearing error at: sample/{dataset}/{sample}")
        path = os.path.join(self.working_dir, "samples", dataset, sample)
        for f in os.listdir(path):
            os.remove(os.path.join(path, f))

    def clean(self):
        if not self.active:
            return

        error_query = self.database.session.query(Sample.id, 
                                                  Sample.parentDataset,
                                                  Sample.sampleName,
                                                  Error.errorCode)\
                                           .filter(Sample.id == Error.sampleId)
        errors = self.database.safe_run(error_query.all)
        for e in errors:
            try:
                self._del_sample_path(e[1], e[2])
                self.database.session.query(Error).filter(Error.sampleId == e[0]).delete()
                self.database.session.commit()
            except Exception as err:
                print("====> Something went wrong during clearing")
                self.database.session.rollback()
