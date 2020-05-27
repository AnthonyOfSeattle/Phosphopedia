checkpoint update_database:
    output:
        touch("flags/database_flags/database_updated.flag")
    run:
        database = DatabaseInterface(DATABASE_PATH)
        database.initialize_database()
        database.update_datasets(
            config["datasets"]
        )

rule obtain_data:
    input:
        ancient("flags/database_flags/database_updated.flag")
    output:
        "samples/{parentDataset}/{sampleName}/{sampleName}.raw"
    params:
        prefix = lambda wildcards: wildcards.parentDataset[:3],
        location = lambda wildcards : DatabaseInterface(DATABASE_PATH).get_sample(sample_name = wildcards.sampleName).fileLocation,
        download_error = "samples/{parentDataset}/{sampleName}/download.error"
    group:
        "preprocess"
    shell:
        """
        trap "touch {output} {params.download_error}; exit 0" ERR

        echo {params.location}
        if [ "{params.prefix}" = "LOC" ] 
        then
          ln -s {params.location} {output}
        else if [ "{params.prefix}" = "PXD" ]
        then
          wget -O {output} {params.location}
          fi
        fi
        """

rule raw_to_mzml:
    input:
        ancient("samples/{parentDataset}/{sampleName}/{sampleName}.raw")
    output:
        "samples/{parentDataset}/{sampleName}/{sampleName}.mzML"
    log:
        "samples/{parentDataset}/{sampleName}/{sampleName}.raw_to_mzml.log.txt"
    conda:
        SNAKEMAKE_DIR + "/envs/raw_parser.yaml"
    benchmark:
        "benchmarks/raw_to_mzml/{parentDataset}/{sampleName}.benchmark.txt"
    params:
        conversion_error = "samples/{parentDataset}/{sampleName}/conversion.error"
    group:
        "preprocess"
    shell:
        """
        trap "touch {output} {params.conversion_error}; exit 0" ERR

        ThermoRawFileParser.sh -i {input} \
                               -b {output} \
                               -f 2 \
        &> {log}
        """

def update_mass_analyzer(input_file, sample_name, database):
    from pyteomics.mzml import MzML

    mz = MzML(input_file)

    ms1_analyzer = ""
    ms2_analyzer = ""

    for scan in mz:
        analyzer = scan["scanList"]["scan"][0]["filter string"][:4]
        if scan["ms level"] == 1:
            ms1_analyzer = analyzer
        elif scan["ms level"] == 2:
            ms2_analyzer = analyzer

        if ms1_analyzer and ms2_analyzer:
            break

    database.add_sample_params(sample_name = sample_name,
                               ms1_analyzer = ms1_analyzer,
                               ms2_analyzer = ms2_analyzer)

rule file_checks:
    input:
        ancient("samples/{parentDataset}/{sampleName}/{sampleName}.mzML")
    output:
        touch("flags/preprocess_flags/{parentDataset}/{sampleName}.preprocess.complete")
    params:
        download_error = "samples/{parentDataset}/{sampleName}/download.error",
        conversion_error = "samples/{parentDataset}/{sampleName}/conversion.error"
    group:
        "preprocess"
    run:
        import os
        database = DatabaseInterface(DATABASE_PATH)
        if os.path.isfile(params.download_error) or os.path.isfile(params.conversion_error):
            database.add_error(sample_name = wildcards.sampleName, error_code = "preprocessError")
        else:
            update_mass_analyzer(input[0], wildcards.sampleName, database)

def evaluate_samples_to_preprocess(wildcards):
    checkpoints.update_database.get(**wildcards)
    return expand(
        "flags/preprocess_flags/{parentDataset}/{sampleName}.preprocess.complete",
        zip, **DatabaseInterface(DATABASE_PATH).get_sample_manifest() #.sample(100, random_state = 0)
    )

checkpoint finalize_preprocessing:
    input:
        evaluate_samples_to_preprocess
    output:
        touch("flags/preprocess_flags/preprocess.complete")
