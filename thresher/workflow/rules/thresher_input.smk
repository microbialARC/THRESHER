rule thresher_input:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        hc_groups = os.path.join(config["output"], "thresher", "input", "hierarchical_clustering_groups.RDS"),
        study_snp_matrix =  os.path.join(config["output"],"mummer4_study", "study_snp_matrix.RDS"),
        global_snp_matrix = os.path.join(config["output"],"mummer4_global", "global_snp_matrix.RDS"),
        iqtree_group_path = os.path.join(config["output"], "iqtree","group_tree","iqtree_group.txt")
    params:
        group_tree_dir = os.path.join(config["output"], "iqtree","group_tree"),
        thresher_input_dir = os.path.join(config["output"], "thresher", "input"),
        use_cladebreaker = config["use_cladebreaker"]
    threads:
        config["threads"]
    output:
       thresher_input = os.path.join(config["output"], "thresher", "input", "thresher_input.RDS")
    script:
        os.path.join(BASE_PATH,"scripts","thresher_input.R")