#!/bin/bash

echo Running profiling on Velos data
python openms_ascore.py --fragment_mass_tolerance .25 \
                        --max_peptide_length 55 \
                        --n_to_profile 10000 \
                        --profile_file_name data/09143.profile.tsv \
                        data/09143.{mzML,pep.xml,ascore.txt}

echo Running profiling on QExactive data
python openms_ascore.py --fragment_mass_tolerance .25 \
                        --max_peptide_length 55 \
                        --n_to_profile 10000 \
                        --profile_file_name data/q04348.profile.tsv \
                        data/q04348.{mzML,pep.xml,ascore.txt}

