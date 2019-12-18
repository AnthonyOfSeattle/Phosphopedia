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
        mzml = ancient("samples/{parentDataset}/{sampleName}/{sampleName}.mzML"),
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
        protected("percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.target.pep.xml"),
        protected("percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.decoy.pep.xml"),
        protected("percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.target.psms.txt"),
        protected("percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.decoy.psms.txt")
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
        percolator_ids = "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.{psmLabel}.psms.txt",
        parameter_file = ancient("ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.params")
    output:
        protected("ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.{psmLabel}.txt")
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

####################################
#                                  #
# Write all PSM scores to database #
#                                  #
####################################

rule write_to_database:
    input:
        psm_scores = "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.{psmLabel}.psms.txt",
        localizations = "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.{psmLabel}.txt"
    output:
        touch("flags/search_flags/{parentDataset}/{sampleName}.search.{psmLabel}.write.complete")
    run:
        import re
        import time
        import numpy as np
        import pandas as pd

        localizations = pd.read_csv(input.localizations, sep="\t").set_index("Scan")
        psm_scores = pd.read_csv(input.psm_scores, sep="\t",
                                 usecols=["scan",
                                          "percolator score",
                                          "percolator q-value",
                                          "percolator PEP"]).set_index("scan", drop=False)
        localized_scores = psm_scores.join(localizations, how="left") #.sort_values(by="scan").join(localizations, how="left")
        localized_scores = localized_scores[~localized_scores.PepScore.isna()]
        print(localized_scores.head())
        
        def yield_mod(seq):
            for ind, mod in enumerate(re.finditer("[A-Zn](\[[^A-Z]+\])?", seq), 1):
                mod_mass = re.search("(?<=\[)([^A-Z]*)(?=\])", mod.group())
                if mod_mass is not None:
                    # Subtract 1 if the first character is an n-terminal designation
                    yield mod.group()[0], ind - int(seq[0] == "n"), float(mod_mass.group())

        def create_entries(row):
            psm = PSM(sample_name = wildcards.sampleName,
                      scan_number = row["scan"],
                      psm_label = wildcards.psmLabel,
                      base_sequence =  re.sub("[^A-Z]+", "", row["LocalizedSequence"]),
                      psm_score = row["percolator score"],
                      psm_qvalue = row["percolator q-value"],
                      psm_pep = row["percolator PEP"])

            phospho_scores = iter(zip(row["Ascores"].split(";"), str(row["AltSites"]).split(";")))

            for mod_res, mod_ind, mod_mass in yield_mod(row["LocalizedSequence"]):
                score = np.inf
                alt_pos = ""
                if np.isclose(config["ascore"]["params"]["mod_mass"], mod_mass, rtol=0, atol=1):
                    score, alt_pos = next(phospho_scores)

                psm.modifications.append(Modification(
                    residue=mod_res,
                    position=mod_ind, 
                    mass=mod_mass, 
                    localization_score=score, 
                    alternative_positions=alt_pos
                ))

            return psm

        engine = create_engine(ALCHEMY_PATH)
        session = sessionmaker(bind=engine)()
        entries = localized_scores.apply(create_entries, axis=1).tolist()
        for i in map(time.sleep, np.random.uniform(0, 5, size=100)):
            try:
                session.add_all(localized_scores.apply(create_entries, axis=1).tolist())
                session.commit()
                break
            except OperationalError:
                session.rollback()

rule finish_search:
    input:
        "flags/search_flags/{parentDataset}/{sampleName}.search.target.write.complete",
        "flags/search_flags/{parentDataset}/{sampleName}.search.decoy.write.complete"
    output:
        touch("flags/search_flags/{parentDataset}/{sampleName}.search.complete")
    shell:
        """
        rm {input[0]} {input[1]}
        """
