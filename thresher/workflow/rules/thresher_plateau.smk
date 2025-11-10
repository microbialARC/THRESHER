rule thresher_plateau:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        thresher_input = os.path.join(config["output"], "thresher", "input","thresher_input.RDS"),
        hc_groups = os.path.join(config["output"], "thresher", "input","hierarchical_clustering_groups.RDS")
    params:
        plateau_length = config["plateau_length"],
        output_dir = os.path.join(config["output"], "thresher", "output")
    threads:
        config["threads"]
    output:
        plateau_strains_rds = os.path.join(config["output"],"thresher", "output", "plateau_strains.RDS"),
        plateau_strains_csv = os.path.join(config["output"], "thresher", "output", "plateau_strains.csv"),
        group_plateau = os.path.join(config["output"], "thresher", "output", "group_plateau.csv")
    script: 
        os.path.join(BASE_PATH,"scripts","thresher_plateau.R")