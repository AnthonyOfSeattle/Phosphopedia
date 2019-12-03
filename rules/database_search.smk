##############################
#                            #
# Database search with comet #
#                            #
##############################

rule comet_param_builder:
    input:
        ancient("flags/preprocess_flags/{parentDataset}/{sampleName}.preprocess.complete")
    output:
        "comet/{parentDataset}/{sampleName}/{sampleName}.comet.params"
    group:
        "comet"
    run:
        ms2_analyzer = SQLiteInterface(
                os.path.join(WORKING_DIR, "project.db")
            ).query_samples(
                "ms2Analyzer", "WHERE parentDataset='{}' AND sampleName='{}'".format(wildcards.parentDataset, wildcards.sampleName)
            ).iloc[0]["ms2Analyzer"]

        param_lines = []
        for key, value in config["comet"]["params"].items():
            if key == ms2_analyzer and isinstance(value, dict):
                param_lines.extend([" = ".join([k, str(v)]) + "\n" for k,v in value.items()])
            elif isinstance(value, dict):
                continue
            else:
                param_lines.append(" = ".join([key, str(value)]) + "\n")

        with open(output[0], "w") as dest:
            dest.writelines(param_lines)

rule comet_search:
    input:
        mzml = "samples/{parentDataset}/{sampleName}/{sampleName}.mzML",
        parameter_file = ancient("comet/{parentDataset}/{sampleName}/{sampleName}.comet.params"),
        ref = "config/" + config["comet"]["ref"]
    output:
        pep_xml = protected("comet/{parentDataset}/{sampleName}/{sampleName}.comet.target.pep.xml"),
        pin = protected("comet/{parentDataset}/{sampleName}/{sampleName}.comet.target.pin")
    params:
        fileroot = "{sampleName}",
        output_dir = "comet/{parentDataset}/{sampleName}"
    conda:
        SNAKEMAKE_DIR + "/envs/crux.yaml"
    benchmark:
        "benchmarks/comet/{parentDataset}/{sampleName}.benchmark.txt"
    group:
        "comet"
    shell:
        """
        crux comet --parameter-file {input.parameter_file} \
                   --fileroot {params.fileroot} \
                   --output-dir {params.output_dir} \
                   --output_txtfile 0 \
                   --output_percolatorfile 1 \
                   --output_pepxmlfile 1 \
                   --overwrite T \
                   {input.mzml} {input.ref}
        """

###############################
#                             #
# PSM scoring with Percolator #
#                             #
###############################

rule percolator_pin_builder:
    input:
        "comet/{parentDataset}/{sampleName}/{sampleName}.comet.target.pin"
    output:
        "percolator/{parentDataset}/{sampleName}/{sampleName}.pin"
    params:
        drop = ["deltCn", "deltLCn"]
    group:
        "percolator"
    script:
        "scripts/clean_comet_pin_files.py"

rule percolator_scoring:
    input:
        ancient("percolator/{parentDataset}/{sampleName}/{sampleName}.pin")
    output:
        pep_xml = protected("percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.target.pep.xml"),
        csv = protected("percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.target.psms.txt")
    params:
        fileroot = "{sampleName}",
        output_dir = "percolator/{parentDataset}/{sampleName}"
    conda:
        SNAKEMAKE_DIR + "/envs/crux.yaml"
    benchmark:
        "benchmarks/percolator/{parentDataset}/{sampleName}.benchmark.txt"
    group:
        "percolator"
    shell:
        """
        crux percolator --fileroot {params.fileroot} \
                        --output-dir {params.output_dir} \
                        --only-psms T \
                        --pepxml-output T \
                        --top-match 1 \
                        --overwrite T \
                        {input}
        """

############################################################
#                                                          #
# Localize phosphorylations in target peptides with Ascore #
#                                                          #
############################################################

rule ascore_param_builder:
    input:
        "comet/{parentDataset}/{sampleName}/{sampleName}.comet.target.pep.xml"
    output:
        "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.params"
    group:
        "ascore"
    run:
        ms2_analyzer = SQLiteInterface(
                os.path.join(WORKING_DIR, "project.db")
            ).query_samples(
                "ms2Analyzer", "WHERE parentDataset='{}' AND sampleName='{}'".format(wildcards.parentDataset, wildcards.sampleName)
            ).iloc[0]["ms2Analyzer"]

        param_lines = []
        for key, value in config["ascore"]["params"].items():
            if key == ms2_analyzer and isinstance(value, dict):
                param_lines.extend([" = ".join([k, str(v)]) + "\n" for k,v in value.items()])
            elif isinstance(value, dict):
                continue
            else:
                param_lines.append(" = ".join([key, str(value)]) + "\n")

        with open(output[0], "w") as dest:
            dest.writelines(param_lines)

rule ascore_localization:
    input:
        mzml = "samples/{parentDataset}/{sampleName}/{sampleName}.mzML",
        pep_xml = "comet/{parentDataset}/{sampleName}/{sampleName}.comet.target.pep.xml",
        parameter_file = ancient("ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.params")
    output:
        protected("ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.txt")
    conda: 
        SNAKEMAKE_DIR + "/envs/ascore.yaml"
    benchmark:
        "benchmarks/ascore/{parentDataset}/{sampleName}.benchmark.txt"
    group:
        "ascore"
    shell:
        """
        python -m pyascore --parameter_file {input.parameter_file} \
                           {input.mzml} \
                           {input.pep_xml} \
                           {output}
        """
