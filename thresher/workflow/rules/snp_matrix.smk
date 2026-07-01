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


rule public_snp_matrix:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        actual_download_topgenomes = os.path.join(config["output"],"datasets_topgenomes","actual_download_topgenomes.txt"),
        public_concatenated_reports = expand(os.path.join(config["output"], "mummer4_public", "{genome_name}_concatenated.report"),genome_name=genome_path_dict.keys())
    output:
        public_snp_matrix = os.path.join(config["output"], "mummer4_public", "public_snp_matrix.RDS")
    params:
        metadata = config["metadata"],
        report_dir = os.path.join(config["output"],"mummer4_public"),
        snp_coverage_threshold = config["snp_coverage_threshold"],
        whatsgnu_dir = os.path.join(config["output"],"whatsgnu")
    script:
        os.path.join(BASE_PATH,"scripts","public_snp_matrix.R")