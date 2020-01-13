rule promise_data:
    input:
        ancient("project.db")
    group:
        "preprocess"
    output:
        temp(touch("samples/{parentDataset}/{sampleName}/{sampleName}.promise"))

rule obtain_data:
    input:
        ancient("samples/{parentDataset}/{sampleName}/{sampleName}.promise")
    output:
        "samples/{parentDataset}/{sampleName}/{sampleName}.raw"
    params:
        prefix = lambda wildcards: wildcards.parentDataset[:3],
        location = lambda wildcards : DB_INTERFACE.query_samples(
            "fileLocation", "WHERE parentDataset='{}' AND sampleName='{}'".format(wildcards.parentDataset, wildcards.sampleName)
        ).iloc[0]["fileLocation"],
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

def update_errorcode(sample_name):
    interface = SQLiteInterface(
            os.path.join(WORKING_DIR, "project.db")
        )
    interface.update_samples("errorCode", "'PreprocessError'",
                             "WHERE sampleName='{}'".format(sample_name))

def update_mass_analyzer(input_file, sample_name):
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

    interface = SQLiteInterface(
        os.path.join(WORKING_DIR, "project.db")
    )
    interface.update_samples("ms1Analyzer", "'{}'".format(ms1_analyzer),
                             "WHERE sampleName='{}'".format(sample_name))
    interface.update_samples("ms2Analyzer", "'{}'".format(ms2_analyzer),
                             "WHERE sampleName='{}'".format(sample_name))

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
        if os.path.isfile(params.download_error) or os.path.isfile(params.conversion_error):
            update_errorcode(wildcards.sampleName)
        else:
            update_mass_analyzer(input[0], wildcards.sampleName)

checkpoint finalize_preprocessing:
    input:
        expand(
            "flags/preprocess_flags/{parentDataset}/{sampleName}.preprocess.complete", zip,
            **generate_sample_manifest()
        )
    output:
        touch("flags/preprocess_flags/preprocess.complete")
