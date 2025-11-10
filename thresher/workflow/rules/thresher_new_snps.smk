rule thresher_new_snps:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        snp_matrix_new_snps = os.path.join(config["output"], "mummer4", "snp_matrix_new_snps.RDS")
    params:
        thresher_output = config["thresher_output"],
        new_metadata = config["new_metadata"],
        original_metadata = config["original_metadata"],
        new_snp_matrix = os.path.join(config["output"], "mummer4", "snp_matrix_new_snps.RDS"),
        output_dir = os.path.join(config["output"],"thresher","output")
    threads:
        config["threads"]
    output:
        new_plateau_rds = os.path.join(config["output"], "thresher", "output", "new_plateau.RDS"),
        new_plateau_strains_csv = os.path.join(config["output"], "thresher", "output", "new_plateau_strains.csv"),
        new_plateau_genomes_csv = os.path.join(config["output"], "thresher", "output", "new_plateau_genomes.csv"),
        new_peak_rds = os.path.join(config["output"], "thresher", "output", "new_peak.RDS"),
        new_peak_strains_csv = os.path.join(config["output"], "thresher", "output", "new_peak_strains.csv"),
        new_peak_genomes_csv = os.path.join(config["output"], "thresher", "output", "new_peak_genomes.csv"),
        new_discrepancy_rds = os.path.join(config["output"], "thresher", "output", "new_discrepancy.RDS"),
        new_discrepancy_strains_csv = os.path.join(config["output"], "thresher", "output", "new_discrepancy_strains.csv"),
        new_discrepancy_genomes_csv = os.path.join(config["output"], "thresher", "output", "new_discrepancy_genomes.csv"),
        new_global_rds = os.path.join(config["output"], "thresher", "output", "new_global.RDS"),
        new_global_strains_csv = os.path.join(config["output"], "thresher", "output", "new_global_strains.csv"),
        new_global_genomes_csv = os.path.join(config["output"], "thresher", "output", "new_global_genomes.csv")
    script:
        os.path.join(BASE_PATH,"scripts","thresher_new_snps.R")