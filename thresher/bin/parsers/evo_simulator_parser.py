"""Parser for simulation command"""
import argparse
import pandas as pd
import os
# The parser for the overall simulation command
def add_evo_simulator_parser(subparsers):
    """Add the simulation command parser to the subparsers"""
    evo_simulator_parser = subparsers.add_parser(
        "evo_simulator",
        formatter_class=argparse.RawTextHelpFormatter,
        help="Models bacterial evolution through substitution, gene gain/loss, and recombination events." \
        "This tool is modified and adapted from CoreSimul (https://github.com/lbobay/CoreSimul ,DOI: 10.1186/s12859-020-03619-x)" \
        "Please cite the original publication when using this tool." \
        ""
    )

    evo_simulator_parser.add_argument(
        "--preset",
        required=False,
        type = int,
        help="The preset ID for the simulation. Please refer to the documentation for available preset IDs.\n" \
        "If the preset_id is provided, other parameters will be overwritten."
    )

    evo_simulator_parser.add_argument(
        "--ancestor",
        required=True,
        type=str,
        help="Path to the most common ancestor genome in FASTA format used in the simulation."
    )

    evo_simulator_parser.add_argument(
        "--years",
        type=int,
        required=True,
        help="Number of years to simulate evolution."
    )

    evo_simulator_parser.add_argument(
        "--taxa",
        type=int,
        required=True,
        help="Number of taxa(offspring) at the end of the simulation."
    )

    evo_simulator_parser.add_argument(
        "--substitution_model",
        type=str,
        required=True,
        help="Substitution model to use for simulation. Default is GTR."
        "Options are JC69, K2P, K3P, JC, GTR."
    )

    evo_simulator_parser.add_argument(
        "--model_parameters",
        type=str,
        required=True,
        help="""The parameters in the substitution model.
If the substitution model is JC69, no parameters are needed because all substitution rates are equal, and this argument can be left empty.

If the substitution model is K2P, the parameter indicates the transition/transversion ratio(kappa). The parameter should be provided as a single float value greater than 0. 

If the substitution model is K3P, the parameters indicate the transition(alpha) and the two transversion ratios(beta,gamma). The parameters should be provided as three comma-separated values without spaces. 
1st: transition (alpha)
2nd: beta (A↔C, G↔T rates)
3rd: gamma (A↔T, C↔G rates)
Note: Both values should be greater than 0.

If the substitution model is GTR, the paramaters indicate the rates of the six types of nucleotide substitutions.
The parameters should be provided as six comma-separated values without spaces.
The order of parameters should be:
1st: A↔G rate
2nd: A↔C rate
3rd: A↔T rate
4th: G↔C rate
5th: G↔T rate
6th: C↔T rate
        """ 
    )

    evo_simulator_parser.add_argument(
        "--mutation_rate",
        type=float,
        required=True,
        help="Mutation rate (substitutions per site per year) for simulation"
    )

    # Whether or not to use weighted mutation simulation
    evo_simulator_parser.add_argument(
        "--use_weighted_mutation",
        type=bool,
        required=False,
        default=True,
        help="""Whether to use weighted mutation. Default is True.
        To disable weighted mutation, set to False."""
    )

    evo_simulator_parser.add_argument(
        "--weighted_mutation_file",
        type=str,
        required=True,
        help="""Path to the file containing weighted mutation sites.\n
        The file should be a csv file, separated with "," with two columns titled position and entropy.
        Only required if --use_weighted_mutation is set to True."""
    )

    # Whether or not to simulate gene gain/loss
    evo_simulator_parser.add_argument(
        "--use_gain_loss",
        type=bool,
        required=False,
        default=True,
        help="Whether to simulate gene gain/loss events. Default is True." \
        "To disable gene gain/loss simulation, set to False."
    )

    evo_simulator_parser.add_argument(
        "--gain_rate",
        type=float,
        required=True,
        help="Gene gain rate for simulation. Only required if --use_gain_loss is set to True."
    )

    evo_simulator_parser.add_argument(
        "--loss_rate",
        type=float,
        required=True,
        help="Gene loss rate for simulation. Only required if --use_gain_loss is set to True."
    )

    evo_simulator_parser.add_argument(
        "--bin",
        type=str,
        required=True,
        help="Path to the file containing the bins of the genome. Only required if --use_gain_loss is set to True."
    )

    evo_simulator_parser.add_argument(
        "--mge_fasta",
        type=str,
        required=True,
        help="Path to the database containing mobile genetic elements to use for simulation in the fasta format. Only required if --use_gain_loss is set to True."
    )

    evo_simulator_parser.add_argument(
        "--mge_entropy",
        type=str,
        required=True,
        help="Path to the directory containing the entropy of all mobile genetic elements in the database provided with --mge_fasta. Only required if --use_gain_loss is set to True."
    )

    # Whether or not to simulate recombination
    evo_simulator_parser.add_argument(
        "--use_recombination",
        type=bool,
        required=False,
        default=True,
        help="Whether to simulate recombination events. Default is True." \
        "To disable recombination simulation, set to False."
    )

    evo_simulator_parser.add_argument(
        "--recombination_rate",
        type=float,
        required=True,
        help="Recombination rate (event per genome per year) for simulation. Only required if --use_recombination is set to True."
    )

    evo_simulator_parser.add_argument(
        "--min_recombination_size",
        type=int,
        required=True,
        help="Minimum recombination size (in bp) for simulation. Default is 1 bp. Only required if --use_recombination is set to True."
    )

    evo_simulator_parser.add_argument(
        "--mean_recombination_size",
        type=int,
        required=True,
        help="Mean recombination size (in bp) for simulation. Only required if --use_recombination is set to True."
    )

    evo_simulator_parser.add_argument(
        "--nu",
        type=float,
        required=True,
        help="Nu parameter (snps / total length of recombination) for recombination simulation. Only required if --use_recombination is set to True."
    )  

    evo_simulator_parser.add_argument(
        "--output",
        required=False,
        help="Output directory for the simulation results. If not provided, defaults to thresher_evo_simulator_output_<YYYY_MM_DD_HHMMSS> under the current working directory."
    )

    evo_simulator_parser.add_argument(
        "--threads",
        required=False,
        type=int,
        default=os.cpu_count(),
        help="""Number of threads to use for the simulation. Default is maximum cores available.""",
    )

    evo_simulator_parser.add_argument(
        "--seed",
        type=int,
        required=False,
        default=int(pd.Timestamp.now().strftime("%Y%m%d")),
        help="""Random seed for reproducibility.
        default is the current date in the format of YYYYMMDD"""
    )

    evo_simulator_parser.add_argument(
        "--prefix",
        default=None,
        help = "Prefix for config file, output files, and analysis naming. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS"
    )

    evo_simulator_parser.add_argument(
        "--conda_prefix",
        default=None,
        help = "Directory for conda environments needed for this analysis. If not provided, defaults to <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>"
    )