rule plot_persistence_plot_redo_endpoint:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        clusters_summary_rds = [os.path.join(config['output'], "files", "clusters_summary_redo_endpoint.RDS")],
    output:
        persistence_plot_pdf = os.path.join(config["output"],"plots","PersistencePlot_redo_endpoint.pdf"),
        persistence_plot_rds = os.path.join(config["output"],"plots","PersistencePlot_redo_endpoint.RDS")
    params:
        output_dir = config["output"],
        metadata = config["original_metadata"]
    script:
        os.path.join(BASE_PATH,"scripts","plot_persistence_plot.R")