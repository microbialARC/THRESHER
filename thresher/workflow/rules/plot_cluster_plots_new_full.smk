rule plot_cluster_plots_new_full:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        plateau_strains_rds = os.path.join(config["output"], "thresher", "output",  "plateau_strains.RDS"),
        peak_strain_rds = os.path.join(config["output"],  "thresher", "output", "peak_strains.RDS"),
        global_strains_rds = os.path.join(config["output"],  "thresher", "output", "global_strains.RDS"),
        discrepancy_strains_rds = os.path.join(config["output"],  "thresher", "output", "discrepancy_strains.RDS"),
        new_mlst_results = os.path.join(config["output"], "mlst","summary","mlst_results.csv")
    output:
        plot_rds = [os.path.join(config["output"],"plots", "ClusterPlots_new_full.RDS")],
        clusters_summary_rds = [os.path.join(config["output"], "thresher", "output","clusters_summary.RDS")],
        clusters_summary_csv = [os.path.join(config["output"], "thresher", "output","clusters_summary.csv")],
        genomes_summary_csv = [os.path.join(config["output"], "thresher", "output","genomes_summary_new_full.csv")]
    params:
        thresher_output = config["thresher_output"],
        output_dir = config["output"],
        original_mlst_results = os.path.join(config["thresher_output"], "mlst","summary","mlst_results.csv"),
        new_metadata = config["new_metadata"],
        original_metadata = config["original_metadata"]
    script:
        os.path.join(BASE_PATH,"scripts","plot_cluster_plots_new_full.R")