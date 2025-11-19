rule simulator_visualization:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        # Substitution Simulation Outputs
        mutations = os.path.join(config["output"], "evo_simulator", "output", "mutations.csv"),
        weight_mutation = os.path.join(config["output"], "evo_simulator", "intermediate", "weighted_mutation.csv"),
        # MGEs Gain and Loss Simulation Outputs
        gain = os.path.join(config["output"], "evo_simulator", "output", "gain.csv"),
        loss = os.path.join(config["output"], "evo_simulator", "output", "loss.csv"),

        # SNP-dists Outputs
        snp_matrix = os.path.join(config["output"], "evo_simulator", "output", "snp_matrix.txt"),
    params:
        chromosomal_bins = config["bin"]
    output:
        snp_matrix_rds = os.path.join(config["output"], "evo_simulator", "output", "snp_matrix_visualization.RDS"),
        snp_matrix_pdf = os.path.join(config["output"], "evo_simulator", "output", "snp_matrix_visualization.pdf"),
        simulation_visualization = os.path.join(config["output"], "evo_simulator", "output", "simulation_visualization.png")
    script:
        os.path.join(BASE_PATH, "scripts", "simulator_visualization.R")