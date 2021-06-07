rule upload_data:
    input:
        ".pipeline_flags/integration.complete",
        psms = "integration_engine/psms.csv",
        peptides = "integration_engine/peptides.csv",
        peptide_modifications = "integration_engine/peptide_modifications.csv",
        proteins = "integration_engine/proteins.csv",
        peptide_protein = "integration_engine/peptide_protein.csv",
        sites = "integration_engine/sites.csv"
    output:
        touch(".pipeline_flags/upload.complete")
    params:
        ref = "config/" + config["comet"]["ref"],
        ann = config.get("upload", {"annotations" : None})["annotations"]
    run:
        if params.ann is not None:
            params.ann = "config/" + params.ann

        uploader = BuildUploader(
                       DATABASE_PATH, 
                       "integration_engine", 
                       params.ref, 
                       params.ann
                   )
        uploader.convert()
        uploader.upload()
