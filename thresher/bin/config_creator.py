import os
import yaml
# Create configuration files for different THRESHER functions and modes
# Each function has its own config creation function
# The config files are in YAML format
# The config files are saved in the output/config/ directory with a prefix in the format of config_<prefix>.yaml

# strain_identifier function, full-pipeline mode config creators
def strain_identifier_full_config(args):
    """Create config file in yaml format for strain_identifier function"""
    config = {
        'prefix': args.prefix,
        'metadata': args.metadata,
        'output': args.output,
        'species': args.species,
        'mode': args.analysis_mode,
        'whatsgnu_db_path': args.whatsgnu_db_path if args.whatsgnu_db_path else "",
        'bakta_db_type': args.bakta_db_type if args.bakta_db_type else "full",
        'bakta_db_path': args.bakta_db_path if args.bakta_db_path else os.path.join(args.output, "bakta", "db"),
        'core_bootstrap_method': args.core_bootstrap_method,
        'core_bootstrap_number': args.core_bootstrap_number,
        'group_bootstrap_method': args.group_bootstrap_method,
        'group_bootstrap_number': args.group_bootstrap_number,
        'endpoint': args.endpoint,
        'plateau_length': args.plateau_length,
        'threads': args.threads
    }

    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    
    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
# strain_identifier function, redo-endpoint mode config creators
def strain_identifier_redo_endpoint_config(args):
    """Create config file in yaml format for strain_identifier redo-endpoint function"""
    config = {
        'original_metadata': args.original_metadata,
        'thresher_output': args.thresher_output,
        'endpoint': args.endpoint,
        'output': args.output,
        'prefix': args.prefix,
    }

    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
# strain_identifier function, new-snps mode config creators
def strain_identifier_new_snps_config(args):
    """Create config file in yaml format for strain_identifier new-snps function"""
    config = {
        'prefix': args.prefix,
        'species': args.species,
        'new_metadata': args.new_metadata,
        "original_metadata": args.original_metadata,
        'thresher_output': args.thresher_output,
        'output': args.output,
        'threads': args.threads
    }

    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
# strain_identifier function, new-full mode config creators
def strain_identifier_new_full_config(args):
    """Create config file in yaml format for strain_identifier new-full function"""
    config = {
        'prefix': args.prefix,
        'original_metadata': args.original_metadata,
        'new_metadata': args.new_metadata,
        'output': args.output,
        'species': args.species,
        'whatsgnu_db_path': args.whatsgnu_db_path if args.whatsgnu_db_path else "",
        'bakta_db_type': args.bakta_db_type if args.bakta_db_type else "full",
        'bakta_db_path': args.bakta_db_path if args.bakta_db_path else os.path.join(args.output, "bakta", "db"),
        'thresher_output': args.thresher_output,
        'core_bootstrap_method': args.core_bootstrap_method,
        'core_bootstrap_number': args.core_bootstrap_number,
        'group_bootstrap_method': args.group_bootstrap_method,
        'group_bootstrap_number': args.group_bootstrap_number,
        'endpoint': args.endpoint,
        'plateau_length': args.plateau_length,
        'threads': args.threads
    }

    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
# genome_profiler function config creators
def genome_profiler_config(args):
    """Create config file in yaml format for genome_profiler function"""
    config = {
        'prefix': args.prefix,
        'input_genome': args.input_genome,
        'output': args.output,
        'top_genomes': args.top_genomes,
        'ani_threshold': args.ani_threshold,
        'bakta_db_type': args.bakta_db_type if args.bakta_db_type else "full",
        'bakta_db_path': args.bakta_db_path if args.bakta_db_path else os.path.join(args.output, "bakta", "db"),
        'whatsgnu_db_path': args.whatsgnu_db_path if args.whatsgnu_db_path else "",
        'species': args.species,
        'threads': args.threads
    }

    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
# evo_simulator function, preset mode config creators
def evo_simulator_preset_config(args):
    """Create config file in yaml format for evo_simulator preset function"""
    config = {
        'prefix': args.prefix,
        'seed': args.seed,
        'output': args.output,
        'years': args.years,
        'taxa': args.taxa,
        'weighted_mutation': args.weighted_mutation,
        'recombination': args.recombination,
        'gain_loss': args.gain_loss,
        'species': args.species,
        'st': args.st,
        'threads': args.threads
    }

    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
# evo_simulator function, custom mode config creators
def evo_simulator_custom_config(args):
    """Create config file in yaml format for evo_simulator custom function"""
    config = {
        'prefix': args.prefix,
        'seed': args.seed,
        'output': args.output,
        'ancestor_genome': args.ancestor_genome,
        'mutation_rate': args.mutation_rate,
        'recombination_rate': args.recombination_rate,
        'mean_recombination_size': args.mean_recombination_size,
        'gain_rate': args.gain_rate,
        'loss_rate': args.loss_rate,
        'weighted_mutation_file': args.weighted_mutation_file,
        'bin': args.bin,
        'mge_fasta': args.mge_fasta,
        'mge_entropy': args.mge_entropy,
        'threads': args.threads
   }
    
    if not os.path.exists(os.path.join(args.output, "config")):
       os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
       yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")