rule thresher_peak:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        thresher_input = os.path.join(config["output"],"thresher", "input", "thresher_input.RDS"),
        hc_groups = os.path.join(config["output"], "thresher", "input","hierarchical_clustering_groups.RDS")
    params:
        output_dir = os.path.join(config["output"], "thresher", "output")
    output:
        peak_strains_rds = os.path.join(config["output"],"thresher", "output",  "peak_strains.RDS"),
        peak_strains_csv = os.path.join(config["output"], "thresher", "output", "peak_strains.csv"),
        group_peak = os.path.join(config["output"], "thresher", "output", "group_peak.csv")
    script: 
        os.path.join(BASE_PATH,"scripts","thresher_peak.R")