rule unique_topgenomes_new_full:
    input:
        whatsgnu_output = [os.path.join(config["output"],"whatsgnu","whatsgnu_results",f"{genome}_WhatsGNU_topgenomes.txt") for genome in list(new_genome_path_dict.keys())]
    output:
        expected_download_topgenomes = os.path.join(config["output"],"datasets_topgenomes","expected_download_topgenomes.txt")
    params:
        output_dir = os.path.join(config["output"],"datasets_topgenomes"),
        input_dir = os.path.join(config["output"],"whatsgnu"),
        study_accession = new_genome_accession
    script:
        os.path.join(BASE_PATH,"scripts","unique_topgenomes.py")