rule plot_cluster_plots:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        mrsa_results = [os.path.join(config["output"],"blastx","mrsa","output","summary","blastx_MRSA_strains.csv")] if config["species"] == "sau" else [],
        peak_strains_rds = os.path.join(config["output"], "thresher", "output",  "peak_strains.RDS"),
        plateau_strains_rds = os.path.join(config["output"], "thresher", "output",  "plateau_strains.RDS"),
        global_strains_rds = os.path.join(config["output"],  "thresher", "output", "global_strains.RDS"),
        discrepancy_strains_rds = os.path.join(config["output"],  "thresher", "output", "discrepancy_strains.RDS"),
        mlst_results = os.path.join(config["output"], "mlst","summary","mlst_results.csv")
    output:
        plot_rds = os.path.join(config["output"],"plots","ClusterPlots","ClusterPlots.RDS"),
        clusters_summary_rds = os.path.join(config["output"],"thresher","output","clusters_summary.RDS"),
        clusters_summary_csv = os.path.join(config["output"],"thresher","output","clusters_summary.csv")
    params:
        metadata = config["input"],
        species = config["species"],
        endpoint = config["endpoint"],
        output_dir = os.path.join(config["output"],"plots","ClusterPlots")
    script:
        os.path.join(BASE_PATH,"scripts","plot_cluster_plots.R")