"""Parser for genome_profiler command"""
import argparse

# The parser for the overall genome_profiler command
def add_genome_profiler_parser(subparsers):
    genome_profiler_parser = subparsers.add_parser(
        "genome_profiler",
        help="Infer the probability of substitutions and mobile genetic elements."
    )

    genome_profiler_parser.add_argument(
        "--input_genome",
        required=True,
        help="Path to genome assembly in FASTA format"
    )

    genome_profiler_parser.add_argument(
        "--output",
        required=False,
        help="Output directory of the profiling results. If not provided, defaults to thresher_genome_profiler_output_<YYYY_MM_DD_HHMMSS> under the current working directory."
    )
    
    genome_profiler_parser.add_argument(
        "--species",
        type=str,
        required=True,
        choices=["sau", "sepi", "cdiff", "kp"],
        help="""Bacteria species.
        Available options: [sau, sepi, cdiff, kp]
        sau: Staphylococcus aureus
        sepi: Staphylococcus epidermidis
        cdiff: Clostridium difficile
        kp: Klebsiella pneumoniae"""
    )

    genome_profiler_parser.add_argument(
        "--top_genomes",
        type=int,
        default=1000,
        help="Number of top genomes to use for profiling (default: 1000)"
    )

    genome_profiler_parser.add_argument(
        "--ani_threshold",
        type=float,
        default=99.9,
        help="ANI threshold for profiling (default: 99.9)"
    )

    genome_profiler_parser.add_argument(
        "--bakta_db_type",
        required=False,
        default="full",
        help="""Bakta database.
        Available options: [full, light]
        Default is full"""
    )

    genome_profiler_parser.add_argument(
        "--bakta_db_path",
        required=False,
        type=str,
        help="""The path of the directory where the existing Bakta database locates.
        If provided, the Bakta database will not be downloaded.
        If not provided, defaults to <OUTPUT>/bakta/db""",
    )

    genome_profiler_parser.add_argument(
        "--whatsgnu_db_path",
        type=str,
        required=False,
        help="""The path to the existing WhatsGNU database.
        If provided, the WhatsGNU database will not be downloaded.
        If not provided, defaults to <OUTPUT>/whatsgnu/db.""",
    )
    genome_profiler_parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use (default: 4)"
    )

    genome_profiler_parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help="Prefix for config file, output files, and analysis naming. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS"
    )

    genome_profiler_parser.add_argument(
        "--conda_prefix",
        default=None,
        help="Directory for conda environments needed for this analysis. If not provided, defaults to OUTPUT/conda_envs_<YYYY_MM_DD_HHMMSS>"
    )
