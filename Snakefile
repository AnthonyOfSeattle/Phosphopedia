import os
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from snakemake.utils import makedirs, listfiles
from sqlalchemy.exc import OperationalError
include: "interfaces/sqlite_interface.py"
include: "interfaces/database_schema.py"

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
WORKING_DIR = os.getcwd()
DATABASE_PATH = os.path.join(WORKING_DIR, "project.db")

ALCHEMY_PATH = "sqlite:///" + os.path.join(WORKING_DIR, "phosphopedia.db")
try:
    PhosphopediaBase.metadata.create_all(
        create_engine(ALCHEMY_PATH)
    )
except OperationalError:
    print("Database locked, ignoring create.")
print(__name__)
##### load config #####
configfile: "config/config.yaml"

##### target rules #####

localrules: all, clean_pipeline, finish_search, #write_to_database

db_interface = SQLiteInterface(DATABASE_PATH)
rule all:
    input:
        expand("flags/pipeline_flags/{parentDataset}/{sampleName}.pipeline.complete", zip,
               **db_interface.update_datasets(
                     config["datasets"]
                 ).query_samples(["parentDataset", "sampleName"]).sample(100, random_state=0)
              )

rule clean_pipeline:
    input:
        "flags/preprocess_flags/{parentDataset}/{sampleName}.preprocess.complete",
        "flags/search_flags/{parentDataset}/{sampleName}.search.complete",
        #"comet/{parentDataset}/{sampleName}/{sampleName}.comet.target.pep.xml",
        #"percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.target.pep.xml",
        #"ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.target.txt",
        #"ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.decoy.txt"
    output:
        touch("flags/pipeline_flags/{parentDataset}/{sampleName}.pipeline.complete")
    params:
        raw = "samples/{parentDataset}/{sampleName}/{sampleName}.raw",
        mzml = "samples/{parentDataset}/{sampleName}/{sampleName}.mzML"
    shell:
        """
        #echo "Samples retained";
        if [ -f {params.raw} ]; then rm {params.raw}; fi
        if [ -f {params.mzml} ]; then rm {params.mzml}; fi
        """

##### load rules #####
include: "rules/preprocess_data.smk"
include: "rules/database_search.smk"
