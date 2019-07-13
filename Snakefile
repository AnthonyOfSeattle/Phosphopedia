import os
import pandas as pd
from snakemake.utils import makedirs, listfiles

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
WORKING_DIR = os.getcwd()

##### build manifest #####

SAMPLE_MANIFEST = pd.read_csv("sample_manifest", header=0)
SAMPLE_MANIFEST = SAMPLE_MANIFEST.iloc[:1,:]

raw_files = expand("raws/{dataset}/{basename}.raw", zip,
                    dataset=SAMPLE_MANIFEST.dataset,
                    basename=SAMPLE_MANIFEST.basename)

mzml_files = expand("mzmls/{dataset}/{basename}.mzML", zip, 
                    dataset=SAMPLE_MANIFEST.dataset,
                    basename=SAMPLE_MANIFEST.basename)

##### target rules #####

rule all:
    input:
        "qc/stats/version_0_stats.csv",
        mzml_files,
        

##### load rules #####

include: "rules/retrieve_files.smk"
include: "rules/database_search.smk"
include: "rules/generate_report.smk"
