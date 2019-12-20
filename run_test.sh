#! /bin/bash

if [ "$1" == "+" ]; then
    snakemake --use-conda \
              --latency-wait 60 \
              --jobs 16 \
              --cluster-config cluster.json \
              --cluster "qsub -P {cluster.project} -p {cluster.priority} -pe serial {cluster.ncores} -l {cluster.resources} -o {cluster.eNo} -j y" \
              --directory test
else
    snakemake -n --directory test
fi
