rule thresher_plateau_redo_endpoint:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        thresher_input = os.path.join(config["thresher_output"], "thresher", "input","thresher_input.RDS")
    params:
        plateau_length = config["plateau_length"],
        singleton_threshold = config["singleton_threshold"],
        output_dir = os.path.join(config["output"], "files")
    threads:
        1
    output:
        plateau_strains_rds = os.path.join(config["output"], "files", "plateau_strains.RDS"),
        plateau_strains_csv = os.path.join(config["output"], "files", "plateau_strains.csv"),
        group_plateau = os.path.join(config["output"], "files", "group_plateau.csv")
    script:
        os.path.join(BASE_PATH,"scripts","thresher_plateau.R")
