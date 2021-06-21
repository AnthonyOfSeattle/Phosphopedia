##############################
#                            #
# Download and convert files #
#                            #
##############################

rule obtain_data:
    output:
        "samples/{parentDataset}/{sampleName}/{sampleName}.raw"
    params:
        output_path="samples/{parentDataset}/{sampleName}/",
        download_error = "samples/{parentDataset}/{sampleName}/download.error"
    run:
        manager = SampleManager(DATABASE_PATH)
        sample_id = manager.lookup_id(wildcards.parentDataset,
                                      wildcards.sampleName)
        sample = manager.get_sample(sample_id)
        
        try:
            if sample.fileLocation == "remote":
                project = ppx.find_project(wildcards.parentDataset, 
                                           local=params.output_path)
                project.download(sample.fileName)
            else:
                os.symlink(os.path.join(sample.fileLacation, sample.fileName),
                           output[0])
        except Exception as e:
            print("Error occured while downloading file:", e)
            open(output[0], "w").close()
            open(params.download_error, "w").close()

rule raw_to_mzml:
    input:
        "samples/{parentDataset}/{sampleName}/{sampleName}.raw"
    output:
        "samples/{parentDataset}/{sampleName}/{sampleName}.mzML"
    conda:
        SNAKEMAKE_DIR + "/envs/raw_parser.yaml"
    benchmark:
        "benchmarks/raw_to_mzml/{parentDataset}/{sampleName}.benchmark.txt"
    params:
        conversion_error = "samples/{parentDataset}/{sampleName}/conversion.error"
    shell:
        """
        {{
            ThermoRawFileParser.sh -i {input} \
                                   -b {output} \
                                   -f 2
        }} || {{
            touch {output} {params.conversion_error}
        }}
        """

###################################
#                                 #
# Extract information about files #
#                                 #
###################################

def get_number(scan):
    return scan["index"] + 1

def get_level(scan):
    return scan["ms level"]

def get_rt(scan):
    return scan["scanList"]["scan"][0]["scan start time"]

def get_precursor(scan):
    return scan["precursorList"]\
               ["precursor"][0]\
               ["selectedIonList"]\
               ["selectedIon"][0]\
               .get("selected ion m/z", np.nan)

def get_charge(scan):
    return scan["precursorList"]\
               ["precursor"][0]\
               ["selectedIonList"]\
               ["selectedIon"][0]\
               .get("charge state", np.nan)

def get_activation(scan):
    filter_string = scan["scanList"]["scan"][0]["filter string"]
    return re.search("(?<=@)\w+\.\w+", filter_string).group()

def get_analyzer(scan):
    filter_string = scan["scanList"]["scan"][0]["filter string"]
    return re.search("^[A-Z]+", filter_string).group()

rule extract_scan_information:
    input:
        "samples/{parentDataset}/{sampleName}/{sampleName}.mzML"
    output:
        "samples/{parentDataset}/{sampleName}/{sampleName}.scan_info.tsv"
    run:
        from pyteomics.mzml import MzML

        global_data = {}
        scan_keys = ["scan", "ms_level", 
                     "rt", "activation",
                     "precursor", "charge"]
        scan_data = {k: [] for k in scan_keys}
        try:
            mzml = MzML(input[0])
            for scan in mzml:
                scan_data["scan"].append(get_number(scan))
                scan_data["ms_level"].append(get_level(scan))
                scan_data["rt"].append(get_rt(scan))
                if get_level(scan) == 2:
                    global_data.setdefault("ms2_analyzer", get_analyzer(scan))
                    scan_data["activation"].append(get_activation(scan))
                    scan_data["precursor"].append(get_precursor(scan))
                    scan_data["charge"].append(get_charge(scan))
                else:
                    global_data.setdefault("ms1_analyzer", get_analyzer(scan))
                    scan_data["activation"].append('')
                    scan_data["precursor"].append(np.nan)
                    scan_data["charge"].append(np.nan)
        except Exception as e:
            print("Error occured getting scan info:", e)

        
        sample_manager = SampleManager(DATABASE_PATH)
        sample_id = sample_manager.lookup_id(wildcards.parentDataset, 
                                             wildcards.sampleName)
        sample_manager.add_mass_analyzers(sample_id,
                                          global_data.get("ms1_analyzer", ""), 
                                          global_data.get("ms2_analyzer", ""))
        pd.DataFrame(scan_data).to_csv(output[0], sep="\t", index=False)

################################################
#                                              #
# Write output state to finalize preprocessing #
#                                              #
################################################

rule finalize_preprocessing:
    input:
        scan_info = ancient("samples/{parentDataset}/{sampleName}/{sampleName}.scan_info.tsv"),
        mzml = ancient("samples/{parentDataset}/{sampleName}/{sampleName}.mzML")
    output:
        touch(".pipeline_flags/{parentDataset}/{sampleName}/preprocess.complete")
    params:
        download_error = "samples/{parentDataset}/{sampleName}/download.error",
        conversion_error = "samples/{parentDataset}/{sampleName}/conversion.error"
    run:
        if os.path.isfile(params.download_error) or os.path.isfile(params.conversion_error):
            add_sample_error(wildcards.parentDataset,
                             wildcards.sampleName,
                             "preprocessError")
        else:
            add_sample_flag(wildcards.parentDataset,
                            wildcards.sampleName, 
                            "preprocessComplete")

def evaluate_samples_to_preprocess(wildcards):
    checkpoints.update_database.get(**wildcards)
    return expand(
        ".pipeline_flags/{parentDataset}/{sampleName}/preprocess.complete", zip, 
        **SampleManager(DATABASE_PATH).get_incomplete_samples("preprocessComplete")
    )

checkpoint preprocessing_checkpoint:
    input:
        evaluate_samples_to_preprocess
    output:
        touch(".pipeline_flags/preprocess.complete")
