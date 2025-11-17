rule phylogeny_simulation:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        ancestor_path = config["ancestor"]
    output:
        sim_tree = os.path.join(config["output"], "treesim", "simu_tree.nwk")
    params:
        taxa = config["taxa"],
        years = config["years"],
        ancestor_name=ancestor_name,
        seed = config["seed"],
        output_dir=os.path.join(config["output"], "treesim")
    script:
        os.path.join(BASE_PATH,"scripts","treesim.R")