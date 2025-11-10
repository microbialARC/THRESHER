rule plot_cluster_plots_new_snps:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        new_plateau_strains_rds = os.path.join(config["output"], "thresher", "output",  "new_plateau.RDS"),
        new_peak_strain_rds = os.path.join(config["output"],  "thresher", "output", "new_peak.RDS"),
        new_global_strains_rds = os.path.join(config["output"],  "thresher", "output", "new_global.RDS"),
        new_discrepancy_strains_rds = os.path.join(config["output"],  "thresher", "output", "new_discrepancy.RDS"),
        new_mlst_results = os.path.join(config["output"], "mlst","summary","mlst_results.csv")
    output:
        plot_rds = [os.path.join(config["output"],"plots", "ClusterPlots_new_snps.RDS")],
        clusters_summary_rds = [os.path.join(config["output"], "thresher", "output","clusters_summary_new_snps.RDS")],
        clusters_summary_csv = [os.path.join(config["output"], "thresher", "output","clusters_summary_new_snps.csv")],
        genomes_summary_csv = [os.path.join(config["output"], "thresher", "output","genomes_summary_new_snps.csv")]
    params:
        original_plateau_strains_rds = os.path.join(config["thresher_output"], "thresher", "output",  "plateau_strains.RDS"),
        original_peak_strain_rds = os.path.join(config["thresher_output"], "thresher", "output", "peak_strain.RDS"),
        original_global_strains_rds = os.path.join(config["thresher_output"], "thresher", "output", "global_strains.RDS"),
        original_discrepancy_strains_rds = os.path.join(config["thresher_output"], "thresher", "output", "discrepancy_strains.RDS"),
        thresher_output = config["thresher_output"],
        output_dir = config["output"],
        original_mlst_results = os.path.join(config["thresher_output"], "mlst","summary","mlst_results.csv"),
        new_metadata = config["new_metadata"],
        original_metadata = config["original_metadata"]
    script:
        os.path.join(BASE_PATH,"scripts","plot_cluster_plots_new_snps.R")