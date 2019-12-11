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
        location = lambda wildcards : db_interface.query_samples(
            "fileLocation", "WHERE parentDataset='{}' AND sampleName='{}'".format(wildcards.parentDataset, wildcards.sampleName)
        ).iloc[0]["fileLocation"]
    group:
        "preprocess"
    shell:
        """
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
    group:
        "preprocess"
    shell:
        """
        ThermoRawFileParser.sh -i {input} \
                               -b {output} \
                               -f 2 \
        &> {log}
        """

rule finalize_preprocessing:
    input:
        ancient("samples/{parentDataset}/{sampleName}/{sampleName}.mzML")
    output:
        touch("flags/preprocess_flags/{parentDataset}/{sampleName}.preprocess.complete")
    group:
        "preprocess"
    run:
        from pyteomics.mzml import MzML

        mz = MzML(input[0])

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
                                 "WHERE sampleName='{}'".format(wildcards.sampleName))
        interface.update_samples("ms2Analyzer", "'{}'".format(ms2_analyzer), 
                                 "WHERE sampleName='{}'".format(wildcards.sampleName))
