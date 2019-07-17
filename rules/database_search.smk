rule raw_to_mzml:
    input:
        "raws/{dataset}/{basename}.raw"
    output:
        "mzmls/{dataset}/{basename}.mzML"
    log:
        "logs/raw_to_mzml/{dataset}/{basename}.log"
    conda:
        SNAKEMAKE_DIR + "/envs/raw_parser.yaml"
    shell:
        """
        {{ time \
        ThermoRawFileParser.sh -i {input} \
                               -b {output} \
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
        "comet/{dataset}/{basename}.pep.xml"
    log:
        "logs/comet/{dataset}/{basename}.log"
    params:
        fileroot = "{basename}"
    conda:
        SNAKEMAKE_DIR + "/envs/crux.yaml"
    shadow:
        "minimal"
    shell:
        """
        {{ time \
        crux comet --parameter-file {input.parameter_file} \
                   --fileroot {params.fileroot} \
                   {input.mzml} {input.ref} \
        ; }} &> {log}
        mv crux-output/*.pep.xml {output}
        """
