rule thresher_new_full:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        hc_groups = os.path.join(config["output"], "thresher", "input", "hierarchical_clustering_groups.RDS"),
        study_snp_matrix =  os.path.join(config["output"],"mummer4_study", "study_snp_matrix.RDS"),
        global_snp_matrix = os.path.join(config["output"],"mummer4_global", "global_snp_matrix.RDS"),
        iqtree_group_path = os.path.join(config["output"], "iqtree","group_tree","iqtree_group.txt")
    params:
        thresher_output = config["thresher_output"],
        new_metadata = config["new_metadata"],
        original_metadata = config["original_metadata"],
        output_dir = os.path.join(config["output"],"thresher","output"),
        group_tree_dir = os.path.join(config["output"], "iqtree","group_tree"),
        thresher_input_dir = os.path.join(config["output"], "thresher", "input")
    threads:
        config["threads"]
    output:
        updated_plateau_rds = os.path.join(config["output"], "thresher", "output", "plateau_strains.RDS"),
        updated_plateau_strains_csv = os.path.join(config["output"], "thresher", "output", "plateau_strains.csv"),
        updated_plateau_genomes_csv = os.path.join(config["output"], "thresher", "output", "updated_plateau_genomes.csv"),
        updated_peak_rds = os.path.join(config["output"], "thresher", "output", "peak_strains.RDS"),
        updated_peak_strains_csv = os.path.join(config["output"], "thresher", "output", "peak_strains.csv"),
        updated_peak_genomes_csv = os.path.join(config["output"], "thresher", "output", "updated_peak_genomes.csv"),
        updated_discrepancy_rds = os.path.join(config["output"], "thresher", "output", "discrepancy_strains.RDS"),
        updated_discrepancy_strains_csv = os.path.join(config["output"], "thresher", "output", "discrepancy_strains.csv"),
        updated_discrepancy_genomes_csv = os.path.join(config["output"], "thresher", "output", "updated_discrepancy_genomes.csv"),
        updated_global_rds = os.path.join(config["output"], "thresher", "output", "global_strains.RDS"),
        updated_global_strains_csv = os.path.join(config["output"], "thresher", "output", "global_strains.csv"),
        updated_global_genomes_csv = os.path.join(config["output"], "thresher", "output", "updated_global_genomes.csv")
    script:
        os.path.join(BASE_PATH,"scripts","thresher_new_full.R")