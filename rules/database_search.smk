rule raw_to_mzml:
    input:
        "raws/{dataset}/{basename}.raw"
    output:
        "mzmls/{dataset}/{basename}.mzML"
    log:
        "logs/raw_to_mzml/{dataset}/{basename}.log"
    shell:
        """
        {{ time \
        mono {SNAKEMAKE_DIR}/tools/ThermoRawFileParser/ThermoRawFileParser.exe \
            -i {WORKING_DIR}/{input} \
            -b {WORKING_DIR}/{output} \
            -f 2 \
        ; }} &> {WORKING_DIR}/{log}
        """
