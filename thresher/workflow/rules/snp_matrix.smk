rule study_snp_matrix:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        study_concatenated_reports = expand(os.path.join(config["output"], "mummer4_study", "{genome_name}_concatenated.report"),genome_name=genome_path_dict.keys())
    output:
        os.path.join(config["output"],"mummer4_study","study_snp_matrix.RDS")
    params:
        metadata = config["metadata"],
        snp_coverage_threshold = config["snp_coverage_threshold"],
        report_dir = os.path.join(config["output"],"mummer4_study")
    script:
        os.path.join(BASE_PATH,"scripts","study_snp_matrix.R")


rule global_snp_matrix:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        actual_download_topgenomes = os.path.join(config["output"],"datasets_topgenomes","actual_download_topgenomes.txt"),
        global_concatenated_reports = expand(os.path.join(config["output"], "mummer4_global", "{genome_name}_concatenated.report"),genome_name=genome_path_dict.keys())
    output:
        global_snp_matrix = os.path.join(config["output"], "mummer4_global", "global_snp_matrix.RDS")
    params:
        metadata = config["metadata"],
        report_dir = os.path.join(config["output"],"mummer4_global"),
        snp_coverage_threshold = config["snp_coverage_threshold"],
        whatsgnu_dir = os.path.join(config["output"],"whatsgnu")
    script:
        os.path.join(BASE_PATH,"scripts","global_snp_matrix.R")