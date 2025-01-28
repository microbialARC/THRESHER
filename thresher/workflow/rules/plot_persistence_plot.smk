rule plot_persistence_plot:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        clusters_summary_rds = os.path.join(config["output"],"thresher","output","clusters_summary.RDS")
    output:
        persistence_plot_pdf = os.path.join(config["output"], "plots","PersistencePlot.pdf")
    params:
        output_dir = os.path.join(config["output"],"plots"),
        metadata = config["input"]
    script:
        os.path.join(BASE_PATH,"scripts","plot_persistence_plot.R")