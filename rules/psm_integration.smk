###################################################
#                                                 #
# Split proteins into mutally disconnected groups #
#                                                 #
###################################################

def get_protein_group_filename():
    base_filename = os.path.split(config["comet"]["ref"])[1]
    group_filename = ".".join(
        base_filename.split(".")[:-1] + ["protein_groups", "json"]
    )

    return os.path.join("integration_engine", group_filename)

rule group_proteins:
    input:
        ref = "config/" + config["comet"]["ref"]
    output:
        group_file = get_protein_group_filename()
    run:
        # List percolator files
        sample_manifest = SampleManager(DATABASE_PATH).get_sample_manifest()
        percolator_files = expand(
            expand(
                "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.{{psmLabel}}.psms.txt", 
                zip, **sample_manifest
            ), psmLabel=["target", "decoy"]
        )
        percolator_files.sort()
        pg = PercolatorProteinGrouper(n_groups=config["integration"]["ngroups"], 
                                      n_workers=8, 
                                      chunk_size=1e3,
                                      fdr_filter=config["integration"]["fdr_filter"])
        pg.group(percolator_files)
        pg.to_json(output.group_file)

########################################
#                                      #
# Integrate PSMs into final detections #
#                                      #
########################################

rule subintegration:
    input:
        ref = "config/" + config["comet"]["ref"],
        group_file = get_protein_group_filename()
    output:
        flag = touch(".pipeline_flags/subintegration/subintegration_{groupNumber}.complete"),
        psms = "integration_engine/subintegration/group_{groupNumber}/psms.csv",
        peptides = "integration_engine/subintegration/group_{groupNumber}/peptides.csv",
        peptide_modifications = "integration_engine/subintegration/group_{groupNumber}/peptide_modifications.csv",
        proteins = "integration_engine/subintegration/group_{groupNumber}/proteins.csv",
        peptide_protein = "integration_engine/subintegration/group_{groupNumber}/peptide_protein.csv",
        sites = "integration_engine/subintegration/group_{groupNumber}/sites.csv"
    params:
        output_dir="integration_engine/subintegration/group_{groupNumber}/"
    benchmark:
        "benchmarks/integration_manager/subintegration_{groupNumber}.benchmark.txt"
    run:
        # List files which contain scan information
        scan_info_files = []
        percolator_files = []
        ascore_files = []
        sample_manifest = SampleManager(DATABASE_PATH).get_sample_manifest()
        for file_ind in range(sample_manifest.shape[0]):
            file_info = sample_manifest.iloc[file_ind, :]
            
            scan_info_files.extend(
              2 * ["samples/{parentDataset}/{sampleName}/{sampleName}.scan_info.tsv".format(**file_info)]
            )

            percolator_files.extend(
                expand(
                    "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.{psmLabel}.psms.txt",
                    psmLabel=["target", "decoy"],
                    **file_info
                )
            )

            ascore_files.extend(
                expand(
                    "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.{psmLabel}.txt",
                    psmLabel=["target", "decoy"],
                    **file_info
                )
           )
        
        # Feed files to subintegration     
        manager = SubintegrationManager(nworkers=8, 
                                        file_chunk_size=1e3, 
                                        record_chunk_size=1e4,
                                        fdr_filter=config["integration"]["fdr_filter"]
                                       )

        t0 = time.time()
        manager.map_psms(
            scan_info_files, percolator_files, ascore_files,
            input.group_file, int(wildcards.groupNumber)
        )
        print("PSMs took {} seconds".format(time.time() - t0))

        t0 = time.time()
        manager.map_peptides()
        print("Peptides took {} seconds".format(time.time() - t0))

        manager.read_in_fasta(input.ref)

        t0 = time.time()
        manager.infer_protein_coverage()
        print("Coverage estimation took {} seconds".format(time.time() - t0))

        t0 = time.time()
        manager.drop_low_coverage()
        print("High coverage selection took {} seconds".format(time.time() - t0))

        t0 = time.time()
        manager.map_modifications()
        print("Modifications took {} seconds".format(time.time() - t0))

        t0 = time.time()
        manager.dump(params.output_dir)
        print("Dump took {} seconds".format(time.time() - t0))


rule finalize_integration:
    input:
        expand(".pipeline_flags/subintegration/subintegration_{groupNumber}.complete",
               groupNumber = range(config["integration"]["ngroups"]))
    output:
        psms = "integration_engine/psms.csv",
        peptides = "integration_engine/peptides.csv",
        peptide_modifications = "integration_engine/peptide_modifications.csv",
        proteins = "integration_engine/proteins.csv",
        peptide_protein = "integration_engine/peptide_protein.csv",
        sites = "integration_engine/sites.csv"
    params:
        dirs = expand("integration_engine/subintegration/group_{groupNumber}",
                      groupNumber = range(config["integration"]["ngroups"]))
    benchmark:
        "benchmarks/integration_manager/finalize_integration.benchmark.txt"
    run:
        print("Merging data")
        merger = SubintegrationMerger("integration_engine")
        merger.process(params.dirs)

        print("Calculating FDR on full data")
        fdr_calculator = FDRCalculator()
        fdr_calculator.process_path("integration_engine/")

        print("Calculating the iRT on full data")
        rt_calculator = RTCalculator("isotonic",
                                     98115,
                                     method="descent",
                                     lr=1e-2,
                                     max_epochs=25)
        rt_calculator.process_path("integration_engine/")

##############################################
#                                            #
# Write output state to finalize integration #
#                                            #
##############################################

def wait_on_search(wildcards):
    checkpoints.search_checkpoint.get(**wildcards)
    return ["integration_engine/psms.csv",
            "integration_engine/peptides.csv",
            "integration_engine/peptide_modifications.csv",
            "integration_engine/proteins.csv",
            "integration_engine/peptide_protein.csv",
            "integration_engine/sites.csv"]

checkpoint integration_checkpoint:
    input:
        wait_on_search
    output:
        touch(".pipeline_flags/integration.complete")
