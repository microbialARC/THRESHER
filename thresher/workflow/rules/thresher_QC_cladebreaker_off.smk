rule thresher_QC_cladebreaker_off:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        # core tree and snp matrix from previous thresher results
        core_gene_tree = os.path.join(config["thresher_output"],"iqtree","core_gene_tree","core_gene_tree.contree"),
        study_snp_matrix = os.path.join(config["thresher_output"],"mummer4_study","study_snp_matrix.RDS"),
        # strain compositions from current thresher results with cladebreaker off
        peak_strains_rds = os.path.join(config["output"], "thresher", "output",  "peak_strains.RDS"),
        plateau_strains_rds = os.path.join(config["output"], "thresher", "output",  "plateau_strains.RDS"),
        public_strains_rds = os.path.join(config["output"],  "thresher", "output", "public_strains.RDS"),
        discrepancy_strains_rds = os.path.join(config["output"],  "thresher", "output", "discrepancy_strains.RDS")
    output:
        plateau_qc_plot = os.path.join(config["output"], "thresher", "output", "QC", "plateau_qc_plot.pdf"),
        peak_qc_plot = os.path.join(config["output"], "thresher", "output", "QC", "peak_qc_plot.pdf"),
        public_qc_plot = os.path.join(config["output"], "thresher", "output", "QC", "public_qc_plot.pdf"),
        discrepancy_qc_plot = os.path.join(config["output"], "thresher", "output", "QC", "discrepancy_qc_plot.pdf"),
        plateau_qc_table = os.path.join(config["output"], "thresher", "output", "QC", "plateau_qc_table.csv"),
        peak_qc_table = os.path.join(config["output"], "thresher", "output", "QC", "peak_qc_table.csv"),
        public_qc_table = os.path.join(config["output"], "thresher", "output", "QC", "public_qc_table.csv"),
        discrepancy_qc_table = os.path.join(config["output"], "thresher", "output", "QC", "discrepancy_qc_table.csv")
    params:
        QC_dir = os.path.join(config["output"],"thresher", "output", "QC")
    script:
        os.path.join(BASE_PATH,"scripts","thresher_QC.R")