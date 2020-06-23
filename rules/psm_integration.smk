protein_group_filename = os.path.split(config["comet"]["ref"])[1]
protein_group_filename = ".".join(protein_group_filename.split(".")[:-1]) + ".protein_groups"
rule group_proteins:
    input:
        ref = "config/" + config["comet"]["ref"],
        flag = ancient("flags/search_flags/search.complete")
    output:
        group_file = "integration_engine/" + protein_group_filename
    run:
        pg = ProteinGrouper(
            PeptideExtractor(input.ref, "trypsin"),
            config["integration"]["ngroups"]
        )
        pg.group()
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

        sample_manifest = DatabaseInterface(DATABASE_PATH).get_sample_manifest()
        percolator_files = expand(
            expand(
                "percolator/{parentDataset}/{sampleName}/{sampleName}.percolator.{{psmLabel}}.psms.txt", zip, **sample_manifest
            ), psmLabel=["target", "decoy"]
        )
        ascore_files = expand(
            expand(
                "ascore/{parentDataset}/{sampleName}/{sampleName}.ascore.{{psmLabel}}.txt", zip, **sample_manifest
            ), psmLabel=["target", "decoy"]
        )
        #print(percolator_files[:10])
        #print(ascore_files[:10])
        #quit()

        #print("Finding files")
        #percolator_files = glob("percolator/LOC00000[1234]/*/*.percolator.*.psms.txt")
        #percolator_files += glob("percolator/PXD000293/*/*.percolator.*.psms.txt")
        #percolator_files += glob("percolator/PXD000612/*/*.percolator.*.psms.txt")
        #percolator_files.sort(
        #    key=lambda path: (os.path.split(path)[1].split(".")[0], 
        #                      os.path.split(path)[1].split(".")[2])
        #)

        #ascore_files = glob("ascore/LOC00000[1234]/*/*.ascore.*.txt")
        #ascore_files += glob("ascore/PXD000293/*/*.ascore.*.txt")
        #ascore_files += glob("ascore/PXD000612/*/*.ascore.*.txt")
        #ascore_files.sort(
        #    key=lambda path: (os.path.split(path)[1].split(".")[0], 
        #                      os.path.split(path)[1].split(".")[2])
        #)
        
        manager = SubintegrationManager(
            nworkers=8, file_chunk_size=1e3, record_chunk_size=1e4
        )

        t0 = time.time()
        manager.map_psms(
            percolator_files, ascore_files,
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
        manager.update_peptide_fdr()
        print("Peptide FDR update took {} seconds".format(time.time() - t0))
    
        t0 = time.time()
        manager.update_ptm_fdr()
        print("PTM FDR update took {} seconds".format(time.time() - t0))

        t0 = time.time()
        manager.dump(params.output_dir)
        print("Dump took {} seconds".format(time.time() - t0))


rule finalize_integration:
    input:
        expand("flags/integration_flags/subintegration_{groupNumber}.complete",
               groupNumber = range(config["integration"]["ngroups"]))
    output:
        touch("flags/integration_flags/integration.complete")
    params:
        dirs = expand("integration_engine/subintegration/group_{groupNumber}",
                      groupNumber = range(config["integration"]["ngroups"]))
    benchmark:
        "benchmarks/integration_manager/finalize_integration.benchmark.txt"
    run:
        merger = SubintegrationMerger("integration_engine")
        merger.process(params.dirs)
