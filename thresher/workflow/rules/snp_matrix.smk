rule study_snp_matrix:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        expand(os.path.join(config["output"], "mummer4_study", "{genome_name}_concatenated.report"),genome_name=genome_path_dict.keys())
    output:
        os.path.join(config["output"],"mummer4_study","study_snp_matrix.RDS")
    params:
        report_dir = os.path.join(config["output"],"mummer4_study")
    script:
        os.path.join(BASE_PATH,"scripts","study_snp_matrix.R")


rule global_snp_matrix:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        expand(os.path.join(config["output"], "mummer4_global", "{genome_name}_concatenated.report"),genome_name=genome_path_dict.keys())
    output:
        os.path.join(config["output"], "mummer4_global", "global_snp_matrix.RDS")
    params:
        report_dir = os.path.join(config["output"],"mummer4_global"),
        whatsgnu_dir = os.path.join(config["output"],"whatsgnu")
    script:
        os.path.join(BASE_PATH,"scripts","global_snp_matrix.R")