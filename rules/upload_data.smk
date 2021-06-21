def wait_on_integration(wildcards):
    checkpoints.integration_checkpoint.get(**wildcards)
    return ["integration_engine/psms.csv",
            "integration_engine/peptides.csv",
            "integration_engine/peptide_modifications.csv",
            "integration_engine/proteins.csv",
            "integration_engine/peptide_protein.csv",
            "integration_engine/sites.csv"]

rule upload_data:
    input:
        wait_on_integration
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
