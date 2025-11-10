"""Parser for simulation command"""
import argparse
import pandas as pd
import os
# The parser for the overall simulation command
def add_evo_simulator_parser(subparsers):
    """Add the simulation command parser to the subparsers"""

    evo_simulator_parser = subparsers.add_parser(
        "evo_simulator",
        help="Simulate bacterial genome evolution with substitution, gene gain/loss, and recombination events."
    )

    evo_simulator_subparsers = evo_simulator_parser.add_subparsers(
        dest="mode",
        required=True,
        help="Modes: 'preset' for preset parameters and 'custom' for custom configurations using your own genome and parameters."
    )

    # Add subparsers for preset and custom evo_simulator modes
    add_evo_simulator_preset_parser(evo_simulator_subparsers)
    add_evo_simulator_custom_parser(evo_simulator_subparsers)  

def add_evo_simulator_preset_parser(subparsers):
    """Add the evo_simulator preset command parser to the subparsers"""

    evo_simulator_preset_parser = subparsers.add_parser(
        "preset",
        help="Data-driven evolution simulation of bacterial genomes over time using the data inferred by THRESHER Profile, incorporating mutation, recombination, and gene gain/loss events."
    )

    evo_simulator_preset_parser.add_argument(
        "--output",
        required=False,
        help="Path to the output directory. If not provided, defaults to thresher_evo_simulator_output_<YYYY_MM_DD_HHMMSS> under the current working directory."
    )

    evo_simulator_preset_parser.add_argument(
        "--prefix",
        default=None,
        help = "Prefix for config file, output files, and analysis naming. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS"
    )

    evo_simulator_preset_parser.add_argument(
        "--conda_prefix",
        default=None,
        help = "Directory for conda environments needed for this analysis. If not provided, defaults to <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>"
    )

    evo_simulator_preset_parser.add_argument(
        "--seed",
        type=int,
        required=False,
        default=int(pd.Timestamp.now().strftime("%Y%m%d")),
        help="""Random seed for reproducibility.
        default is the current date in the format of YYYYMMDD"""
    )

    evo_simulator_preset_parser.add_argument(
        "--years",
        type=int,
        required=True,
        help="""Evolution time to simulate in years."""
    )

    evo_simulator_preset_parser.add_argument(
        "--taxa",
        type=int,
        required=True,
        help="""Number of taxa() in the phylogeny (offspring) for simulation."""
    )

    evo_simulator_preset_parser.add_argument(
        "--weighted_mutation",
        type=bool,
        required=False,
        default=True,
        help="""Whether to use weighted random mutation site selection.
        Default is True.""",
    )

    evo_simulator_preset_parser.add_argument(
        "--recombination",
        type=bool,
        required=False,
        default=True,
        help="""Whether to use recombination simulation.
        Default is True.""",
    )

    evo_simulator_preset_parser.add_argument(
        "--gain_loss",
        type=bool,
        required=False,
        default=True,
        help="""Whether to use gene gain/loss simulation.
        Default is True."""
    )

    evo_simulator_preset_parser.add_argument(
        "--species",
        required=True,
        choices=["sau", "kp", "cdiff"],
        help="""Preset name for simulation.
            Available options: [sau, kp, cdiff]
            sau: Staphylococcus aureus
            kp: Klebsiella pneumoniae
            cdiff: Clostridium difficile"""
    )

    evo_simulator_preset_parser.add_argument(
        "--st",
        required=True,
        type=int,
        help="""Sequence Type (ST) to use for the selected species.
        Specify the ST number for the simulation preset.""",
    )

    evo_simulator_preset_parser.add_argument(
        "--threads",
        required=False,
        type=int,
        default=os.cpu_count(),
        help="""Number of threads to use for the simulation. Default is maximum cores available.""",
    )

def add_evo_simulator_custom_parser(subparsers):
    """Add custom subparser"""

    evo_simulator_custom_parser = subparsers.add_parser(
        "custom",
        help="Evolution simulation of bacterial genomes over time with full customized configuration.",
    )

    evo_simulator_custom_parser.add_argument(
        "--seed",
        required=False,
        type=int,
        default=int(pd.Timestamp.now().strftime("%Y%m%d")),
        help="""Random seed for reproducibility.
        default is the current date in the format of YYYYMMDD"""
    )

    evo_simulator_custom_parser.add_argument(
        "--output",
        required=False,
        help="Path to the output directory. If not provided, defaults to thresher_evo_simulator_output_<YYYY_MM_DD_HHMMSS> under the current working directory."
    )


    evo_simulator_custom_parser.add_argument(
        "--ancestor_genome",
        required=True,
        help="Path to genome assembly as the ancestor for simulation in FASTA format"
    )

    evo_simulator_custom_parser.add_argument(
        "--prefix",
        default=None,
        help = "Prefix for config file, output files, and analysis naming. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS"
    )

    evo_simulator_custom_parser.add_argument(
        "--conda_prefix",
        default=None,
        help = "Directory for conda environments needed for this analysis. If not provided, defaults to <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>"
    )

    evo_simulator_custom_parser.add_argument(
        "--mutation_rate",
        type=float,
        required=True,
        help="Mutation rate (substitutions per site per year) for simulation"
    )

    evo_simulator_custom_parser.add_argument(
        "--recombination_rate",
        type=float,
        required=True,
        help="Recombination rate (event per genome per year) for simulation"
    )

    evo_simulator_custom_parser.add_argument(
        "--mean_recombination_size",
        type=int,
        required=True,
        help="Mean recombination size (in bp) for simulation"
    )

    evo_simulator_custom_parser.add_argument(
        "--gain_rate",
        type=float,
        required=True,
        help="Gene gain rate for simulation"
    )

    evo_simulator_custom_parser.add_argument(
        "--loss_rate",
        type=float,
        required=True,
        help="Gene loss rate for simulation"
    )

    evo_simulator_custom_parser.add_argument(
        "--weighted_mutation_file",
        type=str,
        required=True,
        help="""Path to the file containing weighted mutation sites.\n
        The file should be a csv file, separated with "," with two columns titled position and entropy."""
    )

    evo_simulator_custom_parser.add_argument(
        "--bin",
        type=str,
        required=True,
        help="""Path to the file containing the bins of the genome."""
    )

    evo_simulator_custom_parser.add_argument(
        "--mge_fasta",
        type=str,
        required=True,
        help=""""Path to the database containing mobile genetic elements to use for simulation in the fasta format."""
    )

    evo_simulator_custom_parser.add_argument(
        "--mge_entropy",
        type=str,
        required=True,
        help="""Path to the directory containing the entropy of all mobile genetic elements in the database provided with --mge_fasta."""
    )

    evo_simulator_custom_parser.add_argument(
        "--threads",
        required=False,
        type=int,
        default=os.cpu_count(),
        help="""Number of threads to use for the simulation. Default is maximum cores available.""",
    )