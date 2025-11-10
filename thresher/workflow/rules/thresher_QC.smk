rule thresher_QC:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        comprehensive_tree = os.path.join(config["output"],"iqtree","comprehensive_tree","comprehensive_tree.contree"),
        study_snp_matrix = os.path.join(config["output"],"mummer4_study","study_snp_matrix.RDS"),
        peak_strains_rds = os.path.join(config["output"], "thresher", "output",  "peak_strains.RDS"),
        plateau_strains_rds = os.path.join(config["output"], "thresher", "output",  "plateau_strains.RDS"),
        global_strains_rds = os.path.join(config["output"],  "thresher", "output", "global_strains.RDS"),
        discrepancy_strains_rds = os.path.join(config["output"],  "thresher", "output", "discrepancy_strains.RDS")
    output:
        plateau_qc_plot = os.path.join(config["output"], "thresher", "output", "QC", "plateau_qc_plot.pdf"),
        peak_qc_plot = os.path.join(config["output"], "thresher", "output", "QC", "peak_qc_plot.pdf"),
        global_qc_plot = os.path.join(config["output"], "thresher", "output", "QC", "global_qc_plot.pdf"),
        discrepancy_qc_plot = os.path.join(config["output"], "thresher", "output", "QC", "discrepancy_qc_plot.pdf"),
        plateau_qc_table = os.path.join(config["output"], "thresher", "output", "QC", "plateau_qc_table.csv"),
        peak_qc_table = os.path.join(config["output"], "thresher", "output", "QC", "peak_qc_table.csv"),
        global_qc_table = os.path.join(config["output"], "thresher", "output", "QC", "global_qc_table.csv"),
        discrepancy_qc_table = os.path.join(config["output"], "thresher", "output", "QC", "discrepancy_qc_table.csv")
    params:
        QC_dir = os.path.join(config["output"],"thresher", "output", "QC")
    script:
        os.path.join(BASE_PATH,"scripts","thresher_QC.R")