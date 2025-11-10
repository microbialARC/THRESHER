rule thresher_discrepancy:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        thresher_input = os.path.join(config["output"], "thresher", "input", "thresher_input.RDS"),
        hc_groups = os.path.join(config["output"], "thresher", "input", "hierarchical_clustering_groups.RDS")
    params:
        output_dir = os.path.join(config["output"], "thresher", "output")
    threads:
        config["threads"]
    output:
        discrepancy_strains_rds = os.path.join(config["output"], "thresher", "output", "discrepancy_strains.RDS"),
        discrepancy_strains_csv = os.path.join(config["output"],"thresher", "output",  "discrepancy_strains.csv"),
        group_discrepancy = os.path.join(config["output"], "thresher", "output", "group_discrepancy.csv")
    script: 
        os.path.join(BASE_PATH,"scripts","thresher_discrepancy.R")