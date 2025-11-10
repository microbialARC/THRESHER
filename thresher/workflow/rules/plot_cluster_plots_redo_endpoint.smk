rule plot_cluster_plots_redo_endpoint:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        peak_strains_rds = os.path.join(config["thresher_output"], "thresher", "output",  "peak_strains.RDS"),
        plateau_strains_rds = os.path.join(config["thresher_output"], "thresher", "output",  "plateau_strains.RDS"),
        global_strains_rds = os.path.join(config["thresher_output"],  "thresher", "output", "global_strains.RDS"),
        discrepancy_strains_rds = os.path.join(config["thresher_output"],  "thresher", "output", "discrepancy_strains.RDS"),
        mlst_results = os.path.join(config["thresher_output"], "mlst","summary","mlst_results.csv")
    output:
        plot_rds = [os.path.join(config["output"], "plots", "ClusterPlots_redo_endpoint.RDS")],
        clusters_summary_rds = [os.path.join(config["output"], "files", "clusters_summary_redo_endpoint.RDS")],
        clusters_summary_csv = [os.path.join(config["output"], "files", "clusters_summary_redo_endpoint.csv")]
    params:
        thresher_output = config["thresher_output"],
        metadata = config["original_metadata"],
        endpoint = config["endpoint"],
        output_dir = config["output"]
    script:
        os.path.join(BASE_PATH,"scripts","plot_cluster_plots_redo_endpoint.R")