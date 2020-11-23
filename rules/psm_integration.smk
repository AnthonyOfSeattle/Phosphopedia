protein_group_filename = os.path.split(config["comet"]["ref"])[1]
protein_group_filename = ".".join(protein_group_filename.split(".")[:-1]) + ".protein_groups"
rule group_proteins:
    input:
        ref = "config/" + config["comet"]["ref"],
        flag = ancient("flags/search_flags/search.complete")
    output:
        group_file = "integration_engine/" + protein_group_filename
    run:
        # List percolator files
        sample_manifest = DatabaseInterface(DATABASE_PATH).get_sample_manifest()
        percolator_files = expand(
            expand(
                "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.{{psmLabel}}.psms.txt", zip, **sample_manifest
            ), psmLabel=["target", "decoy"]
        )
        percolator_files.sort()
        pg = PercolatorProteinGrouper(n_groups = config["integration"]["ngroups"], 
                                      n_workers = 8, 
                                      chunk_size=1e3)
        pg.group(percolator_files)
        pg.to_json(output.group_file)


rule subintegration:
    input:
        ref = "config/" + config["comet"]["ref"],
        group_file = "integration_engine/" + protein_group_filename,
        flag = ancient("flags/search_flags/search.complete")
    output:
        flag = touch("flags/integration_flags/subintegration_{groupNumber}.complete"),
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
        from glob import glob

        # List files which contain scan information
        sample_manifest = DatabaseInterface(DATABASE_PATH).get_sample_manifest()
        scan_info_files = 2 * expand(
            "samples/{parentDataset}/{sampleName}/{sampleName}.scan_info.tsv", zip, **sample_manifest
        )
        scan_info_files.sort()

        # List percolator files
        percolator_files = expand(
            expand(
                "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.{{psmLabel}}.psms.txt", zip, **sample_manifest
            ), psmLabel=["target", "decoy"]
        )
        percolator_files.sort()

        # List corresponding ascore files
        ascore_files = expand(
            expand(
                "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.{{psmLabel}}.txt", zip, **sample_manifest
            ), psmLabel=["target", "decoy"]
        )
        ascore_files.sort()

        manager = SubintegrationManager(
            nworkers=8, file_chunk_size=1e3, record_chunk_size=1e4
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
        expand("flags/integration_flags/subintegration_{groupNumber}.complete",
               groupNumber = range(config["integration"]["ngroups"]))
    output:
        touch("flags/integration_flags/integration.complete"),
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
                                     max_epochs=50)
        rt_calculator.process_path("integration_engine/")
