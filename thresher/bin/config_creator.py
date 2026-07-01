import os
import yaml
# Create configuration files for different THRESHER modes
# Each mode has its own config creation function
# The config files are in YAML format
# The config files are saved in the output/config/ directory with a prefix in the format of config_<prefix>.yaml

# full-pipeline mode config creators
def full_config(args):
    """Create config file in yaml format for full-pipeline mode"""
    config = {
        'prefix': args.prefix,
        'metadata': args.metadata,
        'output': args.output,
        'species': args.species,
        'epi_mode': args.epi_mode,
        'whatsgnu_db_path': args.whatsgnu_db_path,
        'bakta_db_type': args.bakta_db_type,
        'bakta_db_path': args.bakta_db_path,
        'snp_coverage_threshold': args.snp_coverage_threshold,
        'core_threshold': args.core_threshold,
        'core_bootstrap_method': args.core_bootstrap_method,
        'core_bootstrap_number': args.core_bootstrap_number,
        'group_bootstrap_method': args.group_bootstrap_method,
        'group_bootstrap_number': args.group_bootstrap_number,
        'threshold_floor': args.threshold_floor,
        'threshold_ceiling': args.threshold_ceiling,
        'singleton_threshold': args.singleton_threshold,
        'correction_bootstrap': args.correction_bootstrap,
        'use_cladebreaker': args.use_cladebreaker,
        'endpoint': args.endpoint,
        'plateau_length': args.plateau_length,
        'threads': args.threads
    }

    # If thread is 1, give warning about long runtime for large datasets
    if args.threads == 1:
        print("Warning: Using a single thread may result in long runtimes. Consider using multiple threads for better performance.")
    
    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    
    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
# redo-endpoint mode config creators
def redo_endpoint_config(args):
    """Create config file in yaml format for redo-endpoint mode"""
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
# new-snps mode config creators
def new_snps_config(args):
    """Create config file in yaml format for new-snps mode"""
    config = {
        'prefix': args.prefix,
        'species': args.species,
        'new_metadata': args.new_metadata,
        "original_metadata": args.original_metadata,
        'snp_coverage_threshold': args.snp_coverage_threshold,
        'thresher_output': args.thresher_output,
        'output': args.output,
        'threads': args.threads
    }

    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
# new-full mode config creators
def new_full_config(args):
    """Create config file in yaml format for new-full mode"""
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
        'snp_coverage_threshold': args.snp_coverage_threshold,
        'core_threshold': args.core_threshold,
        'core_bootstrap_method': args.core_bootstrap_method,
        'core_bootstrap_number': args.core_bootstrap_number,
        'group_bootstrap_method': args.group_bootstrap_method,
        'group_bootstrap_number': args.group_bootstrap_number,
        'threshold_floor': args.threshold_floor,
        'threshold_ceiling': args.threshold_ceiling,
        'singleton_threshold': args.singleton_threshold,
        'correction_bootstrap': args.correction_bootstrap,
        'use_cladebreaker': args.use_cladebreaker,
        'endpoint': args.endpoint,
        'plateau_length': args.plateau_length,
        'threads': args.threads
    }
    # If thread is 1, give warning about long runtime for large datasets
    if args.threads == 1:
        print("Warning: Using a single thread may result in long runtimes. Consider using multiple threads for better performance.")
    
    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
# cladebreaker_off mode config creators
def cladebreaker_off_config(args):
    """Create config file in yaml format for cladebreaker_off mode"""
    config = {
        'prefix': args.prefix,
        'output': args.output,
        'thresher_output': args.thresher_output,
        'output': args.output,
        'threshold_floor': args.threshold_floor,
        'threshold_ceiling': args.threshold_ceiling,
        'singleton_threshold': args.singleton_threshold,
        'correction_bootstrap': args.correction_bootstrap,
        'plateau_length': args.plateau_length,
        'use_cladebreaker': False,
        'threads': args.threads
    }
    
    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
