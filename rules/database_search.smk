from itertools import chain

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

rule percolator_correction:
  input:
    "comet/{dataset}/{basename}.pep.xml"
  output:
    target_mzid = "percolator/{dataset}/{basename}.target.mzid",
    decoy_mzid = "percolator/{dataset}/{basename}.decoy.mzid",
    target_pep_xml = "percolator/{dataset}/{basename}.target.pep.xml",
    decoy_pep_xml = "percolator/{dataset}/{basename}.decoy.pep.xml",
    target_tsv = "percolator/{dataset}/{basename}.target.psms.txt",
    decoy_tsv = "percolator/{dataset}/{basename}.decoys.psms.txt",
    weights = "percolator/{dataset}/{basename}.percolator.weights.txt"
  log:
    "logs/percolator/{dataset}/{basename}.log"
  params:
    fileroot = "{basename}"
  conda:
    SNAKEMAKE_DIR + "/envs/crux.yaml"
  shadow:
    "minimal"
  shell:
    """
    {{ time \
    crux percolator --only-psms T \
                    --output-weights T \
                    --fileroot {params.fileroot} \
                    --pepxml-output T \
                    --mzid-output T \
                    --top-match 1 \
                    --verbosity 40 \
                    {input} \
    ; }} &> {log}
    mv crux-output/*.target.mzid {output.target_mzid}
    mv crux-output/*.decoy.mzid {output.decoy_mzid}
    mv crux-output/*.target.pep.xml {output.target_pep_xml}
    mv crux-output/*.decoy.pep.xml {output.decoy_pep_xml}
    mv crux-output/*.target.psms.txt {output.target_tsv}
    mv crux-output/*.decoy.psms.txt {output.decoy_tsv}
    mv crux-output/*.percolator.weights.txt {output.weights}
    """

def get_ascore_params(wildcards):
    comet_kwd = SAMPLE_MANIFEST.loc[(wildcards.dataset, wildcards.basename), "comet"]
    shared_dict = config["ascore"]["params"]["shared"]
    unique_dict = config["ascore"]["params"][comet_kwd]
    arg_list = ["--{} {}".format(a, v) for a, v in chain(shared_dict.items(), unique_dict.items())]
    return " ".join(arg_list)

rule ascore_localization:
  input:
    mzml = "mzmls/{dataset}/{basename}.mzML",
    pep_xml = "comet/{dataset}/{basename}.pep.xml"
  output:
    "ascore/{dataset}/{basename}.ascore.txt"
  log:
    "logs/ascore/{dataset}/{basename}.log"
  params:
    get_ascore_params
  conda: 
    SNAKEMAKE_DIR + "/envs/openms.yaml"
  shell:
    """
    {{ time \
    python {SNAKEMAKE_DIR}/scripts/openms_ascore.py {params} \
                                                    {input.mzml} \
                                                    {input.pep_xml} \
                                                    {output} \
    ; }} &> {log}
    """
