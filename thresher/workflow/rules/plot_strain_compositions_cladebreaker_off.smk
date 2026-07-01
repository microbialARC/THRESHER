rule plot_strain_compositions_cladebreaker_off:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        # snp matrices and group tree from previous thresher results
        study_snp_matrix = os.path.join(config["thresher_output"],"mummer4_study", "study_snp_matrix.RDS"),
        public_snp_matrix = os.path.join(config["thresher_output"],"mummer4_public", "public_snp_matrix.RDS"),
        iqtree_group_path = os.path.join(config["thresher_output"], "iqtree","group_tree","iqtree_group.txt"),
        # strain compositions from current thresher results with cladebreaker off
        plateau_strains_rds = os.path.join(config["output"],"thresher", "output", "plateau_strains.RDS"),
        peak_strains_rds = os.path.join(config["output"],"thresher", "output",  "peak_strains.RDS"),
        public_strains_rds = os.path.join(config["output"], "thresher", "output", "public_strains.RDS"),
        discrepancy_strains_rds = os.path.join(config["output"], "thresher", "output", "discrepancy_strains.RDS")
    params:
        output_dir = os.path.join(config["output"], "plots", "strain_compositions"),
        use_cladebreaker = config["use_cladebreaker"]
    output:
        plateau_strain_tree_snp_rds = os.path.join(config["output"], "plots", "strain_compositions","plateau","plateau_strain_tree_snp.RDS"),
        peak_strain_tree_snp_rds = os.path.join(config["output"], "plots", "strain_compositions","peak","peak_strain_tree_snp.RDS"),
        discrepancy_strain_tree_snp_rds = os.path.join(config["output"], "plots", "strain_compositions","discrepancy","discrepancy_strain_tree_snp.RDS"),
        public_strain_tree_snp_rds = os.path.join(config["output"], "plots", "strain_compositions","public","public_strain_tree_snp.RDS")
    script:
        os.path.join(BASE_PATH,"scripts","plot_strain_compositions.R")
