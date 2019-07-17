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
        ; }} &> {log}
        """

def get_comet_param_file(wildcards):
    comet_kwd = SAMPLE_MANIFEST.loc[(wildcards.dataset, wildcards.basename), "comet"]
    return config["comet"]["param_file"][comet_kwd]


rule comet_search:
    input:
        mzml = "mzmls/{dataset}/{basename}.mzML",
        parameter_file = get_comet_param_file,
        ref = config["comet"]["ref"]
    output:
        "comet_results/{dataset}/{basename}.pep.xml"
    log:
        "logs/comet/{dataset}/{basename}.log"
    params:
        fileroot = "{basename}",
        temp_dir = "/tmp/{dataset}_{basename}_crux_output",
        temp_file = "{basename}.comet.target.pep.xml"
    conda:
        SNAKEMAKE_DIR + "/envs/crux.yaml"
    shell:
        """
        {{ time \
        crux comet --parameter-file {input.parameter_file} \
                   --fileroot {params.fileroot} \
                   --output-dir {params.temp_dir} \
                   --overwrite T \
                   {input.mzml} {input.ref} \
        ; }} &> {log}
        mv {params.temp_dir}/{params.temp_file} {output}
        rm -r {params.temp_dir}
        """
