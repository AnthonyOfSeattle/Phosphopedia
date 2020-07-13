import os
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from snakemake.utils import makedirs, listfiles
from sqlalchemy.exc import OperationalError

from integration_engine import *
include: "interfaces/sqlite_interface.py"
include: "interfaces/database_interface.py"

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
WORKING_DIR = os.getcwd()
DATABASE_PATH = "sqlite:///" + os.path.join(WORKING_DIR, "phosphopedia.db")

##### load config #####
configfile: "config/config.yaml"

##### target rules #####

localrules: all, clean_pipeline, clean_search, finalize_search, finalize_preprocessing, group_proteins, update_database

rule all:
    input:
        "flags/database_flags/database_updated.flag",
        "flags/preprocess_flags/preprocess.complete",
        "flags/search_flags/search.complete",
        "flags/integration_flags/integration.complete",

##### load rules #####
include: "rules/preprocess_data.smk"
include: "rules/database_search.smk"
include: "rules/psm_integration.smk"
