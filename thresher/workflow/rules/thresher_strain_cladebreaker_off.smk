rule thresher_plateau_cladebreaker_off:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        thresher_input = os.path.join(config["output"], "thresher", "input","thresher_input.RDS"),
        hc_groups = os.path.join(config["thresher_output"], "thresher", "input","hierarchical_clustering_groups.RDS")
    params:
        plateau_length = config["plateau_length"],
        singleton_threshold = config["singleton_threshold"],
        output_dir = os.path.join(config["output"], "thresher", "output")
    threads:
        config["threads"]
    output:
        plateau_strains_rds = os.path.join(config["output"],"thresher", "output", "plateau_strains.RDS"),
        plateau_strains_csv = os.path.join(config["output"], "thresher", "output", "plateau_strains.csv"),
        group_plateau = os.path.join(config["output"], "thresher", "output", "group_plateau.csv")
    script: 
        os.path.join(BASE_PATH,"scripts","thresher_plateau.R")

rule thresher_peak_cladebreaker_off:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        thresher_input = os.path.join(config["output"],"thresher", "input", "thresher_input.RDS"),
        hc_groups = os.path.join(config["thresher_output"], "thresher", "input","hierarchical_clustering_groups.RDS")
    params:
        output_dir = os.path.join(config["output"], "thresher", "output")
    threads:
        config["threads"]
    output:
        peak_strains_rds = os.path.join(config["output"],"thresher", "output",  "peak_strains.RDS"),
        peak_strains_csv = os.path.join(config["output"], "thresher", "output", "peak_strains.csv"),
        group_peak = os.path.join(config["output"], "thresher", "output", "group_peak.csv")
    script: 
        os.path.join(BASE_PATH,"scripts","thresher_peak.R")

rule thresher_discrepancy_cladebreaker_off:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        thresher_input = os.path.join(config["output"], "thresher", "input", "thresher_input.RDS"),
        hc_groups = os.path.join(config["thresher_output"], "thresher", "input", "hierarchical_clustering_groups.RDS")
    params:
        output_dir = os.path.join(config["output"], "thresher", "output"),
        singleton_threshold = config["singleton_threshold"]
    threads:
        config["threads"]
    output:
        discrepancy_strains_rds = os.path.join(config["output"], "thresher", "output", "discrepancy_strains.RDS"),
        discrepancy_strains_csv = os.path.join(config["output"],"thresher", "output",  "discrepancy_strains.csv"),
        group_discrepancy = os.path.join(config["output"], "thresher", "output", "group_discrepancy.csv")
    script: 
        os.path.join(BASE_PATH,"scripts","thresher_discrepancy.R")

rule thresher_global_cladebreaker_off:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        thresher_input = os.path.join(config["output"], "thresher", "input","thresher_input.RDS"),
        hc_groups = os.path.join(config["thresher_output"], "thresher", "input","hierarchical_clustering_groups.RDS"),
        global_snp_matrix = os.path.join(config["thresher_output"], "mummer4_global", "global_snp_matrix.RDS")
    params:
        output_dir = os.path.join(config["output"], "thresher", "output"),
        threshold_ceiling = config["threshold_ceiling"]
    threads:
        config["threads"]
    output:
        global_strains_rds = os.path.join(config["output"], "thresher", "output", "global_strains.RDS"),
        global_strains_csv = os.path.join(config["output"], "thresher", "output", "global_strains.csv"),
        group_global = os.path.join(config["output"], "thresher", "output", "group_global.csv")
    script: 
        os.path.join(BASE_PATH,"scripts","thresher_global.R")
        