rule unique_topgenomes:
    input:
        whatsgnu_output = [os.path.join(config["output"],"whatsgnu",f"{genome}",f"{genome}_WhatsGNU_topgenomes.txt") for genome in list(genome_path_dict.keys())]
    output:
        unique_topgenomes = os.path.join(config["output"],"datasets_topgenomes","topgenomes.txt")
    params:
        output_dir = os.path.join(config["output"],"datasets_topgenomes"),
        input_dir = os.path.join(config["output"],"whatsgnu"),
        study_accession = genome_accession
    script:
        os.path.join(BASE_PATH,"scripts","unique_topgenomes.py")