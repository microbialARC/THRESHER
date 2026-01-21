#!/usr/bin/env python3
"""THRESHER Main Entry
This script creates the configuration files and executes the Snakemake workflow.
"""
VERSION = "0.2.3-beta"
# Import standard libraries and custom modules
import argparse
import os
import sys
from pathlib import Path
import traceback
# subprocess is needed to call snakemake
import subprocess
# import shutil to get terminal size for centering text
import shutil

# Import modules for parsers
from thresher.bin.parsers.strain_identifier_parser import add_strain_identifier_parser
from thresher.bin.parsers.genome_profiler_parser import add_genome_profiler_parser
from thresher.bin.parsers.evo_simulator_parser import add_evo_simulator_parser

# Import input validators
from thresher.bin.args_validator import (
    validate_function,
    validate_strain_identifier_full,
    validate_strain_identifier_redo_endpoint,
    validate_strain_identifier_new_snps,
    validate_strain_identifier_new_full,
    validate_genome_profiler,
    validate_evo_simulator,
    ValidationError
)

# Import config creator
from thresher.bin.config_creator import (
    strain_identifier_full_config,
    strain_identifier_redo_endpoint_config,
    strain_identifier_new_snps_config,
    strain_identifier_new_full_config,
    genome_profiler_config,
    evo_simulator_config
)

def build_parser():
    """Build the main argument parser with all subcommands"""
    description = f"""
Main functions:
1. Strain Identifier
2. Genome Profiler
3. Evolution Simulator

Usage:
    thresher <positional arguments>
"""

    parser = argparse.ArgumentParser(
        prog="thresher",
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS
    )

    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"THRESHER version {VERSION}",
        help="show version and exit"
    )
    
    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
        help="Available functions"
    )

    # Add command parsers
    add_strain_identifier_parser(subparsers)
    add_genome_profiler_parser(subparsers)
    add_evo_simulator_parser(subparsers)

    return parser

# Helper functions to check OS and RAM since these checks are used multiple times
def check_os(term_width):
    """Check the operating system"""
    print()
    print("=" * term_width)
    print("Checking operating system".center(term_width))
    print("=" * term_width)
    system_os = sys.platform
    if system_os != 'linux':
        print("THRESHER supports Linux only. Required dependencies are currently available only on Linux.")
        print("Please run THRESHER on a Linux system.")
        print("Thresher quitting...")
        return 1
    else:
        print(f"Current Operating system: {system_os}. Proceeding...")
        return 0

def check_ram(term_width,
              min_ram_gb):
    """Check if the system has at least min_ram_gb of RAM"""
    print()
    print("=" * term_width)
    print("Checking system RAM".center(term_width))
    print("=" * term_width)
    ram_gb =  os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / (1024 ** 3)
    if ram_gb < min_ram_gb:
        print("At least \033[91m40 GB\033[0m memory is required for THRESHER, as WhatsGNU requires loading large databases into memory.")
        print("Please run THRESHER on a system with sufficient RAM.")
        print("Thresher quitting...")
        return 1
    else:
        print(f"Detected RAM: \033[92m{int(ram_gb)} GB\033[0m. Proceeding...")
        return 0
    

def execute_snakemake(config_path, snakefile, conda_prefix, threads, output_dir, prefix):
    """
    Execute the Snakemake workflow
    
    Args:
        config_file: Path to the config YAML file
        snakefile: Name of the Snakefile to use
        conda_prefix: Path to conda environments directory
        threads: Number of threads to use
        output_dir: Output directory for resume command scripts
        prefix: Prefix for the output files
        
    Returns:
        int: return code from Snakemake
    """
    # Locate the workflow directory
    snakefile_path = Path(__file__).parent.parent / "workflow" / snakefile
    # Print execution info
    print(f"Config file: {config_path}")
    print(f"Snakefile: {snakefile_path}")
    print(f"Threads: {threads}")
    print(f"Conda environments: {conda_prefix}")
    
    # Build and execute Snakemake command
    snakemake_cmd = [
        "snakemake",
        "--cores", str(threads),
        "--configfile", str(config_path),
        "--snakefile", str(snakefile_path),
        "--use-conda",
        "--conda-prefix", str(conda_prefix)
    ]
    # If accidentally interrupted, allow resuming using resume_cmd
    # --rerun-incomplete is a flag that added to the snakemake_cmd to allow re-running incomplete jobs
    resume_cmd = [
        "snakemake",
        "--cores", str(threads),
        "--configfile", str(config_path),
        "--snakefile", str(snakefile_path),
        "--use-conda",
        "--conda-prefix", str(conda_prefix),
        "--rerun-incomplete"
    ]
    # First unlock the working directory
    unlock_cmd = [
        "snakemake",
        "--unlock",
        "--configfile", str(config_path),
        "--snakefile", str(snakefile_path),
    ]
    # Create resume_cmd directory
    resume_cmd_dir = Path(output_dir) / "resume"
    resume_cmd_dir.mkdir(parents=True, exist_ok=True)

    # Export unlock command script
    unlock_script_path = resume_cmd_dir / f"unlock_snakemake_{prefix}.sh"
    with open(unlock_script_path, 'w') as f:
        f.write("#!/usr/bin/env bash\n\n")
        f.write("# Unlock Snakemake working directory after interruption\n")
        f.write(' '.join(unlock_cmd) + '\n')
    unlock_script_path.chmod(0o755)
    
    # Export resume command script
    resume_script_path = resume_cmd_dir / f"resume_snakemake_{prefix}.sh"
    with open(resume_script_path, 'w') as f:
        f.write("#!/usr/bin/env bash\n\n")
        f.write("# Resume Snakemake workflow after interruption\n")
        f.write(' '.join(resume_cmd) + '\n')
    resume_script_path.chmod(0o755)

    print(f"Executing Snakemake command: \033[93m{' '.join(snakemake_cmd)}\033[0m")
    # note!! 
    # If the snakemake job is interrupted, to resume, run the same command again
    print("\033[33m" + """
NOTE: To resume an interrupted workflow, simply unlock and run the executed Snakemake command with "--rerun-incomplete" flag, highlighted below.
Snakemake will skip completed steps and continue where it left off.
""" + "\033[0m")
    print(f"To unlock the working directory, run: \033[93m{' '.join(unlock_cmd)}\033[0m")
    print()
    print(f"To resume the workflow, run: \033[93m{' '.join(resume_cmd)}\033[0m")
    print("\033[33mNOTE: This is a temporary workflow resumption method. A dedicated resume feature will be implemented in a future release.\033[0m")
    print()
    print(f"\033[33mUnlock and resume scripts saved to: {unlock_script_path}, {resume_script_path}\033[0m")
    print()
    print("If you are using the same function and mode after the first finished run,")
    print("you can reuse the conda environments by specifying the same --conda-prefix path.")
    print("This avoids re-creating conda environments and saves time and disk space.")
    print(f"Conda environments are stored in: \033[93m{conda_prefix}\033[0m")
    print()
    # Flush stdout to ensure all prints appear before Snakemake output
    sys.stdout.flush()
    sys.stderr.flush()
    # Execute with subprocess.run()
    try:
        result = subprocess.run(snakemake_cmd, check=False)
        return result.returncode
    except FileNotFoundError:
        print("Error: Snakemake command not found. Please ensure Snakemake is installed and accessible in your PATH.")
        return 1
    except Exception as e:
        print(f"Error executing Snakemake: {e}")
        return 1

def main():
    """Main entry point for THRESHER CLI"""

    # Display ASCII art banner
    ascii_art = """
███████╗██╗  ██╗██████╗ ███████╗███████╗██╗  ██╗███████╗██████╗
╚═██╔══╝██║  ██║██╔══██╗██╔════╝██╔════╝██║  ██║██╔════╝██╔══██╗
   ██║   ███████║██████╔╝█████╗  ███████╗███████║█████╗  ██████╔╝
   ██║   ██╔══██║██╔══██╗██╔══╝  ╚════██║██╔══██║██╔══╝  ██╔══██╗
   ██║   ██║  ██║██║  ██║███████╗███████║██║  ██║███████╗██║  ██║
   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚══════╝╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝
   Version {VERSION}
   Center for Microbial Medicine, Children's Hospital of Philadelphia; Philadelphia, PA, USA
    """
    term_width = shutil.get_terminal_size((80, 20)).columns
    ascii_art = ascii_art.format(VERSION=VERSION)
    centered_art = "\n".join(line.center(term_width) for line in ascii_art.splitlines())
    print(centered_art)
    # Build and parse arguments
    parser = build_parser()
    args = parser.parse_args()

    # Print the mode and function being executed
    if args.command == "strain_identifier":
        print(f"Function: {args.command} Mode: {args.mode}".center(term_width))
    elif args.command in ["genome_profiler","evo_simulator"]:
        print(f"Function: {args.command}".center(term_width))

    # Check the OS and RAM for strain_identifier function and full-pipeline and new-full modes
    # This tool can only run on Linux for now
    # Windows and MacOS are not supported due to dependency issues
    
    # If the strain_identifier mode, full-pipeline or new-full, check RAM and OS
    if args.command == "strain_identifier" and args.mode in ["full-pipeline","new-full"]:
        if not args.force:
            # Check OS and abort early on failure
            if (check_os_rc := check_os(term_width)):
                return check_os_rc
            # Check if RAM is enough (at least 40 GB for WhatsGNU database loading)
            # This is required for WhatsGNU database loading
            # the number is highlighted in red for visibility
            if (check_ram_rc := check_ram(term_width, 40)):
                return check_ram_rc
            print(f"RAM check passed for Function: {args.command} Mode: {args.mode}. Proceeding...")

        elif args.force:
            # This elif might seem redundant but it reminds me to handle the --force flag specifically for strain_identifier full-pipeline and new-full modes
            print()
            print("=" * term_width)
            print("Bypassing OS and RAM checks".center(term_width))
            print("=" * term_width)
            print(f"\033[93m--force flag set: skipping OS and RAM checks for function: {args.command} mode: {args.mode}.\nProceeding may lead to WhatsGNU-related and other dependency errors during execution.\033[0m")
    elif args.command == "genome_profiler":
        # For genome_profiler, check OS and RAM
        if not args.force:
            # Check OS and abort early on failure
            if (check_os_rc := check_os(term_width)):
                return check_os_rc
            # Check if RAM is enough (at least 40 GB for WhatsGNU database loading)
            if (check_ram_rc := check_ram(term_width, 40)):
                return check_ram_rc
            print(f"RAM check passed for Function: {args.command}. Proceeding...")
        elif args.force:
            print()
            print("=" * term_width)
            print("Bypassing OS and RAM checks".center(term_width))
            print("=" * term_width)
            print(f"\033[93m--force flag set: skipping OS and RAM checks for function: {args.command}.\nProceeding may lead to dependency errors during execution.\033[0m")
    elif args.command == "evo_simulator":
        # For evo_simulator, only check OS
        if not args.force:
            if (check_os_rc := check_os(term_width)):
                return check_os_rc
        elif args.force:
            print()
            print("=" * term_width)
            print("Bypassing OS check. Proceeding...".center(term_width))
            print("=" * term_width)
    # Validate function arguments
    print()
    print("=" * term_width)
    print("Validating input arguments".center(term_width))
    print("=" * term_width)
    try:
        validate_function(args)
    except Exception as e:
        print(f"Error during validation: {e}")
        return 1

    # Initialize variables
    config_path = None
    snakefile = None
    try:
        # Route to appropriate function parser
        if args.command == "strain_identifier":
            if args.mode == "full-pipeline":
                validate_strain_identifier_full(args)
                config_path = strain_identifier_full_config(args)
                snakefile = "Snakefile_strain_identifier_full"
            elif args.mode == "redo-endpoint":
                validate_strain_identifier_redo_endpoint(args)
                config_path = strain_identifier_redo_endpoint_config(args)
                snakefile = "Snakefile_strain_identifier_redo_endpoint"
                # only use 1 thread for redo-endpoint mode
                args.threads = 1
            elif args.mode == "new-snps":
                validate_strain_identifier_new_snps(args)
                config_path = strain_identifier_new_snps_config(args)
                snakefile = "Snakefile_strain_identifier_new_snps"
            elif args.mode == "new-full":
                validate_strain_identifier_new_full(args)
                config_path = strain_identifier_new_full_config(args)
                snakefile = "Snakefile_strain_identifier_new_full"
            else:
                raise ValidationError(f"Unknown Strain Identifier mode: {args.mode}")
            
        elif args.command == "genome_profiler":
            validate_genome_profiler(args)
            config_path = genome_profiler_config(args)
            snakefile = "Snakefile_genome_profiler"

        elif args.command == "evo_simulator":
            # if preset mode, validate preset-specific args
            validate_evo_simulator(args)
            config_path = evo_simulator_config(args)
            snakefile = "Snakefile_evo_simulator"
            # evo_simulator only uses 1 thread
            args.threads = 1  
        else:
            raise ValidationError(f"Unknown function: {args.command}")
            
        # Ensure config path and snakefile are set
        if config_path is None or snakefile is None:
            raise ValidationError("Configuration path or Snakefile not set properly.")
        # Print a full-width separator
        print()
        print("=" * term_width)
        print("Launching Snakemake workflow".center(term_width))
        print("=" * term_width)

        # Force flush all pending prints before Snakemake execution to ensure proper ordering
        sys.stdout.flush()
        sys.stderr.flush()

        # Execute Snakemake workflow
        # config_path, snakefile, conda_prefix, threads, output_dir, prefix
        returncode = execute_snakemake(
            config_path,
            snakefile,
            os.path.abspath(args.conda_prefix),
            args.threads,
            os.path.abspath(args.output),
            args.prefix
        )

        # Print completion message based on return code
        if returncode == 0:
            if args.command == "strain_identifier":
                print()
                print("=" * term_width)
                print(f"THRESHER function: {args.command} mode:{args.mode} completed successfully!".center(term_width))
                print("=" * term_width)
                print()
            else:
                print()
                print("=" * term_width)
                print(f"THRESHER function: {args.command} completed successfully!".center(term_width))
                print("=" * term_width)
                print()
        else:
            print()
            print("=" * term_width)
            print(f"THRESHER pipeline failed with return code {returncode}".center(term_width))
            print("=" * term_width)
            print()
        
        return returncode
    
    # Handle exceptions and print error messages for debugging purposes
    except ValidationError as e:
        print(f"Validation Error: {e}", file=sys.stderr)
        return 1
    except FileNotFoundError as e:
        print(f"File Not Found Error: {e}", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as e:
        print(f"Snakemake execution failed: {e}", file=sys.stderr)
        return e.returncode
    except Exception as e:
        print(f"Unexpected Error: {e}", file=sys.stderr)
        print("\nTraceback:", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())
    