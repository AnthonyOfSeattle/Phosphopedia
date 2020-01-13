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

##### load config and update database #####
configfile: "config/config.yaml"

DB_INTERFACE = SQLiteInterface(DATABASE_PATH)
DB_INTERFACE.update_datasets(
    config["datasets"]
)

def generate_sample_manifest():
    return DB_INTERFACE.query_samples(["parentDataset", "sampleName"], "WHERE errorCode IS NULL") #.sample(50, random_state=0)

#ALCHEMY_PATH = "sqlite:///" + os.path.join(WORKING_DIR, "phosphopedia.db")
#try:
#    PhosphopediaBase.metadata.create_all(
#        create_engine(ALCHEMY_PATH)
#    )
#except OperationalError:
#    print("Database locked, ignoring create.")


##### load config #####
configfile: "config/config.yaml"

##### target rules #####

localrules: all, clean_pipeline, finish_search, finalize_preprocessing #write_to_database

def evaluate_samples(wildcards):
    checkpoints.finalize_preprocessing.get()
    return expand(
               "flags/pipeline_flags/{parentDataset}/{sampleName}.pipeline.complete", zip,
               **generate_sample_manifest()
           )
    
rule all:
    input:
        evaluate_samples

rule clean_pipeline:
    input:
        "flags/preprocess_flags/preprocess.complete",
        #"flags/preprocess_flags/{parentDataset}/{sampleName}.preprocess.complete",
        #"flags/search_flags/{parentDataset}/{sampleName}.search.complete",
        "comet/{parentDataset}/{sampleName}/{sampleName}.comet.target.pep.xml",
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
        echo "Samples retained";
        #if [ -f {params.raw} ]; then rm {params.raw}; fi
        #if [ -f {params.mzml} ]; then rm {params.mzml}; fi
        """

##### load rules #####
include: "rules/preprocess_data.smk"
include: "rules/database_search.smk"
