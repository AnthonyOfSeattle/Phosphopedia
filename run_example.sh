#!/bin/bash

# Local Version:
# This will run the data sequentially without much fuss
# on the small dataset included it could still take a while.
snakemake --use-conda \
          --cores 1 \
          --latency-wait 120 \
          --directory example/

# Distributed Version:
# This is my usual command for running our larger pipelines
# snakemake --use-conda \
#           --latency-wait 120 \
#           --jobs 8 \
#           --max-jobs-per-second .5 \
#           --cluster-config example_cluster.json \
#           --cluster "qsub -P {cluster.project} -p {cluster.priority} -pe serial {cluster.ncores} -l {cluster.resources} -o {cluster.eNo} -j y" \
#           --directory example/

