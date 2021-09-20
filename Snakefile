import re
import os
import time
import shutil
import ppx
import numpy as np
import pandas as pd
from glob import glob
from snakemake.utils import makedirs, listfiles

from interfaces import *
from integration_engine import *

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
WORKING_DIR = os.getcwd()
DATABASE_PATH = "sqlite:///" + os.path.join(WORKING_DIR, "phosphopedia.db")

##### load config #####
configfile: "config/config.yaml"

##### utility functions #####

def add_sample_error(accession, sample_name, error_code):
    manager = SampleManager(DATABASE_PATH)
    sample_id = manager.lookup_id(accession, sample_name)
    manager.add_error(sample_id, error_code)

def add_sample_flag(accession, sample_name, flag_code):
    manager = SampleManager(DATABASE_PATH)
    sample_id = manager.lookup_id(accession, sample_name)
    manager.add_flag(sample_id, flag_code)

##### target rules #####

localrules:
    all,
    integration_checkpoint,
    search_checkpoint,
    finalize_search,
    preprocessing_checkpoint, 
    update_database

def clean_temp_files():
    shutil.rmtree(".ppx_cache", ignore_errors=True)
    shutil.rmtree(".pipeline_flags", ignore_errors=True)

onstart:
    ErrorCleaner(WORKING_DIR, DATABASE_PATH).clean()

onsuccess:
    clean_temp_files()

onerror:
    clean_temp_files()

rule all:
    input:
        ".pipeline_flags/database.updated",
        ".pipeline_flags/preprocess.complete",
        ".pipeline_flags/search.complete",
        ".pipeline_flags/integration.complete",
        ".pipeline_flags/upload.complete"

checkpoint update_database:
    output:
        touch(".pipeline_flags/database.updated")
    run:
        manager = DatasetManager(DATABASE_PATH)
        manager.add_datasets(config["datasets"])


##### load rules #####
include: "rules/preprocess_data.smk"
include: "rules/database_search.smk"
include: "rules/psm_integration.smk"
include: "rules/upload_data.smk"
