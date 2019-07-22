import os
import pandas as pd
from snakemake.utils import makedirs, listfiles

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
WORKING_DIR = os.getcwd()

##### load config #####
configfile: "config/config.yaml"

##### build manifest #####

SAMPLE_MANIFEST = pd.read_csv("sample_manifest", header=0)
SAMPLE_MANIFEST = SAMPLE_MANIFEST.iloc[[1, 5],:]

raw_files = expand("raws/{dataset}/{basename}.raw", zip,
                    dataset=SAMPLE_MANIFEST.dataset,
                    basename=SAMPLE_MANIFEST.basename)

mzml_files = expand("mzmls/{dataset}/{basename}.mzML", zip, 
                    dataset=SAMPLE_MANIFEST.dataset,
                    basename=SAMPLE_MANIFEST.basename)

comet_files = expand("comet/{dataset}/{basename}.pep.xml", zip,
                     dataset=SAMPLE_MANIFEST.dataset,
                     basename=SAMPLE_MANIFEST.basename)

perc_files = expand("percolator/{dataset}/{basename}.target.pep.xml", zip,
                     dataset=SAMPLE_MANIFEST.dataset,
                     basename=SAMPLE_MANIFEST.basename)

SAMPLE_MANIFEST = SAMPLE_MANIFEST.set_index(["dataset", "basename"])

##### target rules #####

rule all:
    input:
        "qc/stats/version_0_stats.csv",
        mzml_files,
        comet_files,
        perc_files
        

##### load rules #####

include: "rules/retrieve_files.smk"
include: "rules/database_search.smk"
include: "rules/generate_report.smk"
