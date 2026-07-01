rule plot_core_gene_tree:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        contree = os.path.join(config["output"],"iqtree","core_gene_tree","core_gene_tree.contree"),
        hc_groups = os.path.join(config["output"], "thresher", "input", "hierarchical_clustering_groups_simplified.csv"),
        mlst_results = os.path.join(config["output"], "mlst","summary","mlst_results.csv")
    output:
        core_gene_tree_mlst_rds = os.path.join(config["output"],"plots","core_gene_tree_mlst.RDS"),
        core_gene_tree_mlst_pdf = os.path.join(config["output"],"plots","core_gene_tree_mlst.pdf"),
        core_gene_tree_group_rds = os.path.join(config["output"],"plots","core_gene_tree_group.RDS"),
        core_gene_tree_group_pdf = os.path.join(config["output"],"plots","core_gene_tree_group.pdf")
    params:
        output_dir = os.path.join(config["output"],"plots"),
        genome_species = config["species"],
    threads:
        config["threads"]
    script:
        os.path.join(BASE_PATH,"scripts","plot_core_gene_tree.R")