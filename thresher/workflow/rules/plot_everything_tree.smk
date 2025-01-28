rule plot_everything_tree:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        contree = os.path.join(config["output"],"iqtree","everything_tree","everything_tree.contree"),
        hc_groups = os.path.join(config["output"], "thresher", "input", "hierarchical_clustering_groups_simplified.csv"),
        mlst_results = os.path.join(config["output"], "mlst","mlst_results.csv")
    output:
        os.path.join(config["output"],"plots","everything_tree.pdf")
    params:
        output_dir = os.path.join(config["output"],"plots")
    script:
        os.path.join(BASE_PATH,"scripts","plot_everything_tree.R")