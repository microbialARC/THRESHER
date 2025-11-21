rule simulator_simulation:
    conda:
        os.path.join(BASE_PATH, "envs/evo_simulator.yaml")
    input:
        ancestor = config["ancestor"],
        names = os.path.join(config["output"], "evo_simulator", "intermediate", "names.txt"),
        renamed_txt = os.path.join(config["output"], "evo_simulator", "intermediate", "renamed.txt"),
        dichotomies = os.path.join(config["output"], "evo_simulator", "intermediate", "dichotomies.txt"),
        roots = os.path.join(config["output"], "evo_simulator", "intermediate", "roots.txt")
    output:
        stats = os.path.join(config["output"], "evo_simulator", "output", "stats.txt"),
        mutations = os.path.join(config["output"], "evo_simulator", "output", "mutations.csv"),
        weight_mutation = os.path.join(config["output"], "evo_simulator", "intermediate", "weight_mutation.csv") if config["use_weighted_mutation"] else [],
        gain = os.path.join(config["output"], "evo_simulator", "output", "gain.csv"),
        loss = os.path.join(config["output"], "evo_simulator", "output", "loss.csv"),
        recombination = os.path.join(config["output"], "evo_simulator", "output", "recombination.csv"),
        simulated_genomes = os.path.join(config["output"], "evo_simulator", "output", "simulated_genomes.fasta"),
        snp_matrix = os.path.join(config["output"], "evo_simulator", "output", "snp_matrix.txt"),
        chromosome_bins = os.path.join(config["output"], "evo_simulator", "intermediate", "chromosome_bins.csv"),
        nu = os.path.join(config["output"], "evo_simulator", "intermediate", "nu.txt"),
        rm = os.path.join(config["output"], "evo_simulator", "intermediate", "rm.txt")
    params:
        # General Parameters
        prefix = config["prefix"],
        output_dir = os.path.join(config["output"], "evo_simulator", "output"),
        intermediate_dir = os.path.join(config["output"], "evo_simulator", "intermediate"),
        seed = config["seed"],
        # Substitution Simulation Parameters
        mutation_rate = config["mutation_rate"],
        substitution_model = config["substitution_model"],
        model_parameters = config["model_parameters"],
        use_weighted_mutation = config["use_weighted_mutation"],
        weighted_mutation_data = config["weighted_mutation_file"],
        # Recombination Simulation Parameters
        use_recombination = config["use_recombination"],
        recombination_rate = config["recombination_rate"],
        min_recombination_size = config["min_recombination_size"],
        mean_recombination_size = config["mean_recombination_size"],
        nu = config["nu"],
        # MGEs Gain and Loss Simulation Parameters
        use_gain_loss = config["use_gain_loss"],
        gain_rate = config["gain_rate"],
        loss_rate = config["loss_rate"],
        chromosomal_bins = config["bin"],
        mge_fasta = config["mge_fasta"],
        mge_entropy = config["mge_entropy"]
    script:
        os.path.join(BASE_PATH, "scripts", "simulator_simulation.py")