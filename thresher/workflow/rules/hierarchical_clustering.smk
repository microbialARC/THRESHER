rule hierarchical_clustering:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        core_gene_tree_path = os.path.join(config["output"],"iqtree","core_gene_tree","core_gene_tree.treefile"),
        study_snp_matrix_path = os.path.join(config["output"],"mummer4_study","study_snp_matrix.RDS"),
        mlst_results_path = os.path.join(config["output"], "mlst","summary","mlst_results.csv")
    output:
        hc_groups = os.path.join(config["output"], "thresher", "input", "hierarchical_clustering_groups.RDS"),
        hc_groups_csv = os.path.join(config["output"], "thresher", "input", "hierarchical_clustering_groups_simplified.csv")
    params:
        thresher_input_dir = os.path.join(config["output"], "thresher", "input"),
        singleton_threshold = config["singleton_threshold"]
    script:
        os.path.join(BASE_PATH,"scripts","hierarchical_clustering.R")