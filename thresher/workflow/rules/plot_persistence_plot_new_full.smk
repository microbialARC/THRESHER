rule plot_persistence_plot_new_full:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        clusters_summary_rds = [os.path.join(config["output"], "thresher", "output", "clusters_summary.RDS")],
    output:
        persistence_plot_pdf = os.path.join(config["output"],"plots","PersistencePlot_new_full.pdf"),
        persistence_plot_rds = os.path.join(config["output"],"plots","PersistencePlot_new_full.RDS")
    params:
        output_dir = config["output"],
        original_metadata = config["original_metadata"],
        new_metadata = config["new_metadata"]
    script:
        os.path.join(BASE_PATH,"scripts","plot_persistence_plot_new.R")