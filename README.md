## Production pipeline for the Phosphopedia database

Integration of post-translational modification (PTM) mass spectrometry data from public repositories is 
important for building targeted assays, but can be difficult due the scale of data analysis needed.
Since 2016, **Phosphopedia** has housed ~1000 phosphoproteomic samples from internal and public sources.
This repository holds the pipeline that we are using to update this dataset to meet the demands of modern research.
It contains a Snakemake pipeline which automates data retrieval, database search, localization scoring, and psm integration.
If you would like to use this pipeline yourself, it can be downloaded and used with Conda and Snakemake to build
personal databases for your modification of choice.

#### Getting started

The pipeline can effectively be run locally, on a cluster, or distributed on a cloud service such as AWS.
However, it has not been tested on the last case currently.
The only requirements are a Conda distribution installed locally and on the worker nodes,
as well as Snakemake installed in the base Conda environment.
Upon running the pipeline, Snakemake will automatically install the environments located in the env/ folder
into a .snakemake sub-directory of the directory in which you want to do data analysis.

You can verify that Snakemake is up and running using the example configuration.
This contains an example config file which shows most of the available parameters
and points at a small yeast phosphoproteomics dataset which is easy to analyze.
By running a dry run with Snakemake from the Phosphopedia directory, 
you can see the general outline of the pipeline:

```
snakemake -n --directory example/
```

The example ProteomXchange repository contains 12 runs, 
which is fairly simple to run even locally and see what the usual pipeline output would be.
This can be achieved by running the run_example.sh script with the --conda flag.
An example for cluster computation is commented out in run_example.sh as well.

```
snakemake --cores 1 --conda --directory example/
```

#### Future development:

The pipeline can currently be used to generate an FDR controlled final list of
PSMs, peptides, and sites for just about any PTM you can describe a mass for.
In the coming months, we are expecting to add several features:

 - Global dataset post-processing for iRT and charge state distributions.
 - An installable command-line interface that can easily be called in any directory.
 - An integrated web interface to make local browsing easier.
