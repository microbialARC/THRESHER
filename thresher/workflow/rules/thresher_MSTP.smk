rule thresher_MSTP:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        thresher_input = os.path.join(config["output"],  "thresher", "input", "thresher_input.RDS"),
        peak_strains_rds = os.path.join(config["output"], "thresher", "output",  "peak_strains.RDS"),
        plateau_strains_rds = os.path.join(config["output"], "thresher", "output",  "plateau_strains.RDS"),
        global_strains_rds = os.path.join(config["output"],  "thresher", "output", "global_strains.RDS"),
        discrepancy_strains_rds = os.path.join(config["output"],  "thresher", "output", "discrepancy_strains.RDS")
    output:
        MSTP_rds = os.path.join(config["output"], "thresher", "output", "MSTP", "MSTP.RDS")
    params:
        MSTP_dir = os.path.join(config["output"],"thresher", "output", "MSTP")
    script:
        os.path.join(BASE_PATH,"scripts","thresher_MSTP.R")