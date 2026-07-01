rule thresher_public:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        thresher_input = os.path.join(config["output"], "thresher", "input","thresher_input.RDS"),
        hc_groups = os.path.join(config["output"], "thresher", "input","hierarchical_clustering_groups.RDS"),
        public_snp_matrix = os.path.join(config["output"], "mummer4_public", "public_snp_matrix.RDS")
    params:
        output_dir = os.path.join(config["output"], "thresher", "output"),
        threshold_ceiling = config["threshold_ceiling"]
    threads:
        config["threads"]
    output:
        public_strains_rds = os.path.join(config["output"], "thresher", "output", "public_strains.RDS"),
        public_strains_csv = os.path.join(config["output"], "thresher", "output", "public_strains.csv"),
        group_public = os.path.join(config["output"], "thresher", "output", "group_public.csv")
    script: 
        os.path.join(BASE_PATH,"scripts","thresher_public.R")