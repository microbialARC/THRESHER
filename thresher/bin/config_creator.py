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
        'whatsgnu_db_path': args.whatsgnu_db_path,
        'bakta_db_type': args.bakta_db_type,
        'bakta_db_path': args.bakta_db_path,
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
        'whatsgnu_db_path': args.whatsgnu_db_path,
        'bakta_db_type': args.bakta_db_type,
        'bakta_db_path': args.bakta_db_path,
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
        'bakta_db_type': args.bakta_db_type,
        'bakta_db_path': args.bakta_db_path,
        'whatsgnu_db_path': args.whatsgnu_db_path,
        'species': args.species,
        'threads': args.threads
    }

    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
# evo_simulator function config creators
def evo_simulator_config(args):
    """Create config file in yaml format for evo_simulator function"""
    config = {
        'prefix': args.prefix,
        'preset': args.preset,
        'ancestor': args.ancestor,
        'years': args.years,
        'taxa': args.taxa,
        'substitution_model': args.substitution_model,
        'model_parameters': args.model_parameters,
        'mutation_rate': args.mutation_rate,
        'use_weighted_mutation': args.use_weighted_mutation,
        'weighted_mutation_file': args.weighted_mutation_file,
        'use_gain_loss': args.use_gain_loss,
        'gain_rate': args.gain_rate,
        'loss_rate': args.loss_rate,
        'bin': args.bin,
        'position_coverage': args.position_coverage,
        'mge_data': args.mge_data,
        'mge_fasta': args.mge_fasta,
        'mge_entropy': args.mge_entropy,
        'use_recombination': args.use_recombination,
        'recombination_rate': args.recombination_rate,
        'min_recombination_size': args.min_recombination_size,
        'mean_recombination_size': args.mean_recombination_size,
        'nu': args.nu,
        'output': args.output,
        'seed': args.seed    
   }
    
    if not os.path.exists(os.path.join(args.output, "config")):
       os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
       yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")