rule study_snp_matrix_new_snps:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        expand(os.path.join(config["output"], "mummer4", "{genome_name}_concatenated.report"),genome_name=new_genome_path_dict.keys())
    output:
        study_snp_matrix_new = os.path.join(config["output"],"mummer4","study_snp_matrix_new.RDS")
    params:
        report_dir = os.path.join(config["output"],"mummer4"),
        original_snp_matrix = os.path.join(config["thresher_output"], "mummer4_study", "study_snp_matrix.RDS")
    threads:
        config["threads"]
    script:
        os.path.join(BASE_PATH,"scripts","study_snp_matrix_new.R")
