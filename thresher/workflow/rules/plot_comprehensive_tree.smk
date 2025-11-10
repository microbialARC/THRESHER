rule plot_comprehensive_tree:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        contree = os.path.join(config["output"],"iqtree","comprehensive_tree","comprehensive_tree.contree"),
        hc_groups = os.path.join(config["output"], "thresher", "input", "hierarchical_clustering_groups_simplified.csv"),
        mlst_results = os.path.join(config["output"], "mlst","summary","mlst_results.csv")
    output:
        comprehensive_tree_mlst_rds = os.path.join(config["output"],"plots","comprehensive_tree_mlst.RDS"),
        comprehensive_tree_mlst_pdf = os.path.join(config["output"],"plots","comprehensive_tree_mlst.pdf"),
        comprehensive_tree_group_rds = os.path.join(config["output"],"plots","comprehensive_tree_group.RDS"),
        comprehensive_tree_group_pdf = os.path.join(config["output"],"plots","comprehensive_tree_group.pdf")
    params:
        output_dir = os.path.join(config["output"],"plots"),
        genome_species = config["species"],
    threads:
        config["threads"]
    script:
        os.path.join(BASE_PATH,"scripts","plot_comprehensive_tree.R")