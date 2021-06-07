##############################
#                            #
# Database search with comet #
#                            #
##############################

rule comet_param_builder:
    output:
        "comet/{parentDataset}/{sampleName}/{sampleName}.comet.params"
    group:
        "comet"
    run:
        sample_manager = SampleManager(DATABASE_PATH)
        sample_id = sample_manager.lookup_id(wildcards.parentDataset,
                                             wildcards.sampleName)
        params = sample_manager.get_parameters(sample_id)

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
        parameter_file = "comet/{parentDataset}/{sampleName}/{sampleName}.comet.params",
        ref = "config/" + config["comet"]["ref"]
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
        "comet/{parentDataset}/{sampleName}/{sampleName}.comet.target.pin"
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
        "percolator/{parentDataset}/{sampleName}/{sampleName}.pin"
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
        "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.target.psms.txt"
    output:
        "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.params"
    group:
        "ascore"
    run:
        sample_manager = SampleManager(DATABASE_PATH)
        sample_id = sample_manager.lookup_id(wildcards.parentDataset,
                                             wildcards.sampleName)
        params = sample_manager.get_parameters(sample_id)

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
        percolator_ids = "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.{psmLabel}.psms.txt",
        parameter_file = "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.params"
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
        pyascore --parameter_file {input.parameter_file} \
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

rule finalize_search:
    input:
        "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.target.txt",
        "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.decoy.txt"
    output:
        touch(".pipeline_flags/{parentDataset}/{sampleName}/search.complete")
    run:
        raw = "samples/{parentDataset}/{sampleName}/{sampleName}.raw".format(**wildcards)
        if os.path.isfile(raw):
            os.remove(raw)
        
        mzml = "samples/{parentDataset}/{sampleName}/{sampleName}.mzML".format(**wildcards)
        if os.path.isfile(mzml):
            os.remove(mzml)

        add_sample_flag(wildcards.parentDataset,
                        wildcards.sampleName,
                        "searchComplete")

def evaluate_search_samples(wildcards):
    checkpoints.preprocessing_checkpoint.get(**wildcards)
    return expand(
        ".pipeline_flags/{parentDataset}/{sampleName}/search.complete", zip,
        **SampleManager(DATABASE_PATH).get_incomplete_samples("searchComplete")
        )

checkpoint search_checkpoint:
    input:
        evaluate_search_samples
    output:
        touch(".pipeline_flags/search.complete")

