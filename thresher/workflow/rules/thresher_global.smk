rule thresher_global:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        thresher_input = os.path.join(config["output"], "thresher", "input","thresher_input.RDS"),
        hc_groups = os.path.join(config["output"], "thresher", "input","hierarchical_clustering_groups.RDS"),
        global_snp_matrix = os.path.join(config["output"], "mummer4_global", "global_snp_matrix.RDS")
    params:
        output_dir = os.path.join(config["output"], "thresher", "output")
    threads:
        config["threads"]
    output:
        global_strains_rds = os.path.join(config["output"], "thresher", "output", "global_strains.RDS"),
        global_strains_csv = os.path.join(config["output"], "thresher", "output", "global_strains.csv"),
        group_global = os.path.join(config["output"], "thresher", "output", "group_global.csv")
    script: 
        os.path.join(BASE_PATH,"scripts","thresher_global.R")