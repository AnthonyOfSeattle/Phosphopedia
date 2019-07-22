#! /bin/bash

if [ "$1" == "+" ]; then
    snakemake --use-conda --shadow-prefix /tmp --directory test
else
    snakemake -n --use-conda --shadow-prefix /tmp --directory test
fi
