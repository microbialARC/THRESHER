rule simulator_input:
    conda:
        os.path.join(BASE_PATH, "envs/evo_simulator.yaml")
    input:
        simu_tree = os.path.join(config["output"], "treesim", "simu_tree.nwk")
    output:
        names = os.path.join(config["output"], "evo_simulator", "intermediate", "names.txt"),
        renamed_txt = os.path.join(config["output"], "evo_simulator", "intermediate", "renamed.txt"),
        dichotomies = os.path.join(config["output"], "evo_simulator", "intermediate", "dichotomies.txt"),
        roots = os.path.join(config["output"], "evo_simulator", "intermediate", "roots.txt")
    params:
        intermediate_dir = os.path.join(config["output"], "evo_simulator", "intermediate")
    script:
        os.path.join(BASE_PATH,"scripts","simulator_input.py")
        