##############################
#                            #
# Database search with comet #
#                            #
##############################

rule comet_param_builder:
    input:
        ancient("flags/preprocess_flags/preprocess.complete")
    output:
        "comet/{parentDataset}/{sampleName}/{sampleName}.comet.params"
    group:
        "comet"
    run:
        params = DatabaseInterface(DATABASE_PATH).get_sample_params(
            sample_name = wildcards.sampleName
        )

        param_lines = []
        for key, value in config["comet"]["params"].items():
            if key == params.ms2Analyzer and isinstance(value, dict):
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
        ref = ancient("config/" + config["comet"]["ref"])
    output:
        pin = "comet/{parentDataset}/{sampleName}/{sampleName}.comet.target.pin"
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
                   --output_pepxmlfile 0 \
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
        ancient("comet/{parentDataset}/{sampleName}/{sampleName}.comet.target.pin")
    output:
        "percolator/{parentDataset}/{sampleName}/{sampleName}.pin"
    params:
        drop = config.get("percolator_pin_builder", dict(drop=[]))["drop"]
    group:
        "percolator"
    script:
        "scripts/clean_comet_pin_files.py"

rule percolator_scoring:
    input:
        ancient("percolator/{parentDataset}/{sampleName}/{sampleName}.pin")
    output:
        "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.target.pep.xml",
        "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.decoy.pep.xml",
        "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.target.psms.txt",
        "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.decoy.psms.txt"
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
        ancient("percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.target.psms.txt")
    output:
        "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.params"
    group:
        "ascore"
    run:
        params = DatabaseInterface(DATABASE_PATH).get_sample_params(
            sample_name = wildcards.sampleName
        )

        param_lines = []
        for key, value in config["ascore"]["params"].items():
            if key == params.ms2Analyzer and isinstance(value, dict):
                param_lines.extend([" = ".join([k, str(v)]) + "\n" for k,v in value.items()])
            elif isinstance(value, dict):
                continue
            else:
                param_lines.append(" = ".join([key, str(value)]) + "\n")

        with open(output[0], "w") as dest:
            dest.writelines(param_lines)

rule ascore_localization:
    input:
        mzml = ancient("samples/{parentDataset}/{sampleName}/{sampleName}.mzML"),
        percolator_ids = ancient("percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.{psmLabel}.psms.txt"),
        parameter_file = ancient("ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.params")
    output:
        "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.{psmLabel}.txt"
    conda:
        SNAKEMAKE_DIR + "/envs/ascore.yaml"
    benchmark:
        "benchmarks/ascore/{parentDataset}/{sampleName}.{psmLabel}.benchmark.txt"
    group:
        "ascore"
    shell:
        """
        python -m pyascore --parameter_file {input.parameter_file} \
                           --ident_file_type percolatorTXT \
                           {input.mzml} \
                           {input.percolator_ids} \
                           {output}
        """


###########################################
#                                         #
# Write output state to finalize searches #
#                                         #
###########################################


rule clean_search:
    input:
        ancient("ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.target.txt"),
        ancient("ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.decoy.txt")
    output:
        touch("flags/search_flags/{parentDataset}/{sampleName}.search.complete")
    params:
        raw = "samples/{parentDataset}/{sampleName}/{sampleName}.raw",
        mzml = "samples/{parentDataset}/{sampleName}/{sampleName}.mzML"
    shell:
        """
        if [ -f {params.raw} ]; then echo rm {params.raw}; fi
        if [ -f {params.mzml} ]; then echo rm {params.mzml}; fi
        """

def evaluate_search_samples(wildcards):
    checkpoints.finalize_preprocessing.get(**wildcards)
    return expand(
               "flags/search_flags/{parentDataset}/{sampleName}.search.complete", zip,
               **DatabaseInterface(DATABASE_PATH).get_sample_manifest()
           )

checkpoint finalize_search:
    input:
        evaluate_search_samples
    output:
        touch("flags/search_flags/search.complete")

