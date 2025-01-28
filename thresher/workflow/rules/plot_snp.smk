rule plot_snp:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        global_snp_matrix = os.path.join(config["output"], "mummer4_global", "global_snp_matrix.RDS"),
        study_snp_matrix = os.path.join(config["output"], "mummer4_study", "study_snp_matrix.RDS")
    output:
        os.path.join(config["output"],"plots","SNP_Distance.pdf")
    params:
        output_dir = os.path.join(config["output"],"plots")
    script:
        os.path.join(BASE_PATH,"scripts","plot_snp.R")