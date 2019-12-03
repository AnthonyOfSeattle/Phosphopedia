import os
from snakemake.utils import makedirs, listfiles
include: "interfaces/sqlite_interface.py"

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
WORKING_DIR = os.getcwd()
DATABASE_PATH = os.path.join(WORKING_DIR, "project.db")

##### load config #####
configfile: "config/config.yaml"

##### target rules #####

localrules: all, clean_pipeline

db_interface = SQLiteInterface(DATABASE_PATH)
rule all:
    input:
        expand("flags/pipeline_flags/{parentDataset}/{sampleName}.pipeline.complete", zip,
               **db_interface.update_datasets(
                     config["datasets"]
                 ).query_samples(["parentDataset", "sampleName"])
              )

rule clean_pipeline:
    input:
        "flags/preprocess_flags/{parentDataset}/{sampleName}.preprocess.complete",
        "comet/{parentDataset}/{sampleName}/{sampleName}.comet.target.pep.xml",
        "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.target.pep.xml",
        "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.txt"
    output:
        touch("flags/pipeline_flags/{parentDataset}/{sampleName}.pipeline.complete")
    params:
        target = "samples/{parentDataset}/{sampleName}/{sampleName}.mzML"
    shell:
        """
        if [ -f {params.target} ]; then rm {params.target}; fi
        """

##### load rules #####
include: "rules/preprocess_data.smk"
include: "rules/database_search.smk"
