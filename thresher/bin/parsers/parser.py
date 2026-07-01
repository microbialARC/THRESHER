"""Parser for THRESHER command"""
import argparse
import os
# The parser for the overall thresher command
def add_thresher_parser(subparsers):
    """Attach the thresher modes directly to the top-level subparsers."""
    _add_full_pipeline_parser(subparsers)
    _add_cladebreaker_off_parser(subparsers)
    _add_redo_endpoint_parser(subparsers)
    _add_new_snps_parser(subparsers)
    _add_new_full_parser(subparsers)

# The parser for the full-pipeline mode
def _add_full_pipeline_parser(subparsers):
    """
    Add full-pipeline parser to the  subparsers
    """

    full_parser = subparsers.add_parser(
        "full",
        formatter_class=argparse.RawTextHelpFormatter,
        help="Run the complete analysis from genome assemblies with no prior THRESHER output.",
    )


    full_parser.add_argument(
        "--metadata",
        type=str,
        required=True,
        help="""Path to the input metadata file.
The input metadata file should be a tab-delimited file with 3 or 5 columns.
3 columns: genome_name, genome_accession, genome_path
5 columns: genome_name, genome_accession, genome_path, patient_id, collection_date
At least 4 genomes should be provided to perform the analysis."""
    )

    full_parser.add_argument(
        "-o",
        "--output",
        required=False,
        help = "Path to output directory. If not provided, defaults to thresher_output_<YYYY_MM_DD_HHMMSS> under the current working directory."
    )

    full_parser.add_argument(
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

    full_parser.add_argument(
        "--epi_mode",
        required=False,
        type = bool,
        default=True,
        choices=[True, False],
        help="""Whether to perform epidemiological analysis.
Available Options: [True, False]
True: Determine strains, clusters and make plots.
False: Determine strains without determing clusters and making plots.
Default is True"""
    )

    full_parser.add_argument(
        "--whatsgnu_db_path",
        type=str,
        required=False,
        help="""The path to the existing WhatsGNU database.
If provided, the WhatsGNU database will not be downloaded.
If not provided, defaults to <OUTPUT>/whatsgnu/db""",
    )

    full_parser.add_argument(
        "--bakta_db_type",
        required=False,
        default="full",
        help="""The type of Bakta database.
Available options: [full, light]
Default is full"""
    )

    full_parser.add_argument(
        "--bakta_db_path",
        required=False,
        type=str,
        help="""The path of the directory where the existing Bakta database locates.
If provided, the Bakta database will not be downloaded.
If not provided, defaults to <OUTPUT>/bakta/db""",
    )

    full_parser.add_argument(
        "--snp_coverage_threshold",
        required=False,
        type=int,
        default=80,
        help="""Minimum alignment coverage (0-100) required for pairwise SNP distances to be included in analysis. 
Low-coverage alignments can yield unreliable SNP count.
SNP distance below this threshold are excluded. Default: 80."""
    )

    full_parser.add_argument(
        "--core_threshold",
        required=False,
        type=float,
        default=0.95,
        help="""Panaroo Core-genome sample threshold.The frequency of a gene in your sample required to classify it as 'core'.
Range is 0.0 to 1.0. Default is 0.95."""
    )

    full_parser.add_argument(
        "--core_bootstrap_method",
        required=False,
        type=str,
        default="ultrafast",
        choices=["ultrafast", "nonparametric"],
        help="""The bootstrap method for core genome phylogeny used for hierarchical clustering.
Available options: [ultrafast, nonparametric]
ultrafast: Use ultrafast bootstrap method.
nonparametric: Use standard nonparametric bootstrap method.
Default is ultrafast."""
    )

    full_parser.add_argument(
        "--core_bootstrap_number",
        required=False,
        type=int,
        default=0,
        help="""The number of bootstrap replicates for core genome phylogeny.
If method is ultrafast, default is 1000.
If method is nonparametric, default is 100."""
    )

    full_parser.add_argument(
        "--group_bootstrap_method",
        required=False,
        default="ultrafast",
        choices=["ultrafast", "nonparametric"],
        help="""The bootstrap method for phylogeny of each hierarchical group.
Available options: [ultrafast, nonparametric]
ultrafast: Use ultrafast bootstrap method.
nonparametric: Use standard nonparametric bootstrap method.
Default is ultrafast."""
    )

    full_parser.add_argument(
        "--group_bootstrap_number",
        required=False,
        type=int,
        default=0,
        help="""The number of bootstrap replicates for phylogeny of each hierarchical group.
If method is ultrafast, default is 1000.
If method is nonparametric, default is 100."""
    )

    full_parser.add_argument(
        "--threshold_floor",
        required=False,
        type=int,
        default=5,
        help="""The floor of the range of SNP distances to search for the optimal phylothreshold.
Default is 5."""
    )

    full_parser.add_argument(
        "--threshold_ceiling",
        required=False,
        type=int,
        default=500,
        help="""The ceiling of the range of SNP distances to search for the optimal phylothreshold.
Default is 500."""
    )

    full_parser.add_argument(
        "--singleton_threshold",
        required=False,
        type=int,
        default=100,
        help="""The SNP distance threshold above which every genome in the group is considered a singleton.
Default is 100."""
    )

    full_parser.add_argument(
        "--correction_bootstrap",
        required=False,
        type=int,
        default=0,
        help="""Minimum bootstrap support threshold for applying phylogenetic corrections to strain composition (default: 0)."""
    )

    full_parser.add_argument(
        "--use_cladebreaker",
        required=False,
        type=bool,
        default=True,
        choices=[True, False],
        help="""Use CladeBreaker to restrain the strain composition.
Options are [True, False].
Default is True.""")

    full_parser.add_argument(
        "--endpoint",
        required=False,
        default= "plateau",
        choices=["plateau", "peak", "discrepancy", "public"],
        help="""The endpoint method to use for determing clusters and making plots.
Available Options: [plateau, peak, discrepancy, public]
plateau : Phylothreshold set at a plateau where further increases no longer change the number or composition of strains within the group
peak: Phylothreshold set at the peak number of clones defined within the group.
discrepancy: Phylothreshold set at the point where the discrepancy is minimized within the group.
public: Phylothreshold set at the first time a public genome is included in any strain within the group.
Default is plateau."""
    )

    full_parser.add_argument(
        "--plateau_length",
        type=int,
        required=False,
        default=15,
        help="""The plateau length for the plateau endpoint method. Default is 15.
Only used when endpoint method is 'plateau'."""
    )

    full_parser.add_argument(
        "-t",
        "--threads",
        default=1,
        type=int,
        help = "Thread number. Default is 1."
    )

    full_parser.add_argument(
        "--prefix",
        default=None,
        help = "Prefix for config file, output files, and analysis naming. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS"
    )

    full_parser.add_argument(
        "--conda_prefix",
        default=None,
        help = "Directory for conda environments needed for this analysis. If not provided, defaults to <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>"
    )

    # Override system compatibility checks (OS and minimum RAM) and run pipeline regardless
    full_parser.add_argument(
        "--force",
        action="store_true",
        help="""Bypass system compatibility checks (operating system and available RAM) and force execution of the pipeline.
This may cause instability or failures."""
    )
# The parser for the redo-endpoint mode
def _add_redo_endpoint_parser(subparsers):
    """Add redo-endpoint subparser"""
    redo_parser = subparsers.add_parser(
        "redo-endpoint",
        formatter_class=argparse.RawTextHelpFormatter,
        help="Revisit the endpoint choice without rerunning the core algorithm.",
    )

    redo_parser.add_argument(
        "--original_metadata",
        type=str,
        required=True,
        help="""Path to the original input file used for the THRESHER full-pipeline.
The file must be tab-delimited and contain 5 columns because the full-pipeline mode requires patient ID and collection date information.
5 columns: genome_name, genome_accession, genome_path, patient_id, collection_date."""
    )

    redo_parser.add_argument(
        "--thresher_output",
        type=str,
        required=True,
        help="""Path to the existing THRESHER directory.
The existing analysis directory should contain the previous analysis results."""
    )

    redo_parser.add_argument(
        "--output",
        type=str,
        required=False,
        help = """Path to output directory.
If not provided, defaults to thresher_redo_endpoint_<YYYY_MM_DD_HHMMSS> under the current working directory."""
        )
        
    redo_parser.add_argument(
        "--endpoint",
        required=False,
        default= "plateau",
        choices=["plateau", "peak", "discrepancy", "public"],
        help="""The endpoint method to use for determing clusters and making plots.
Available Options: [plateau, peak, discrepancy, public]
plateau : Phylothreshold set at a plateau where further increases no longer change the number or composition of strains within the group
peak: Phylothreshold set at the peak number of clones defined within the group.
discrepancy: Phylothreshold set at the point where the discrepancy is minimized within the group.
public: Phylothreshold set at the first time a public genome is included in any strain within the group.
Default is plateau."""
    )

    redo_parser.add_argument(
        "--prefix",
        default=None,
        help = "Prefix for config files. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS"
    )

    redo_parser.add_argument(
        "--conda_prefix",
        default=None,
        help = "Directory for conda environments needed for this analysis. If not provided, defaults to OUTPUT/conda_envs_<YYYY_MM_DD_HHMMSS>"
    )
# The parser for the new-snps mode
def _add_new_snps_parser(subparsers):
    """Add new-snps subparser"""
    new_snps_parser = subparsers.add_parser(
        "new-snps",
        formatter_class=argparse.RawTextHelpFormatter,
        help=(
            "Add new genomes by reusing the previously inferred phylothresholds."
        ),
    )

    new_snps_parser.add_argument(
        "--original_metadata",
        type=str,
        required=True,
        help="""Path to the input metadata file containing original genomes.
The input file should be a tab-delimited file with 3 or 5 columns.
3 columns: genome_name, genome_accession, genome_path
5 columns: genome_name, genome_accession, genome_path, patient_id, collection_date"""
    )

    new_snps_parser.add_argument(
        "--new_metadata",
        type=str,
        required=True,
        help="""Path to the input metadata file containing new genomes.
At least 1 new genome should be provided to perform new-snps analysis.
The description of the file format is the same as that for --original_metadata."""
    )

    new_snps_parser.add_argument(
        "--thresher_output",
        type=str,
        required=True,
        help="""Path to the existing THRESHER directory.
The existing analysis directory should contain the previous transmission cluster results.
(epi_mode should be 'True' in the prior run)"""
    )

    new_snps_parser.add_argument(
        "--output",
        type=str,
        required=False,
        help = """Path to output directory of new-snps analysis.
If not provided, defaults to thresher_new_snps_<YYYY_MM_DD_HHMMSS> under the current working directory."""
        )
    
    new_snps_parser.add_argument(
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

    new_snps_parser.add_argument(
        "--snp_coverage_threshold",
        required=False,
        type=int,
        default=80,
        help="""Minimum alignment coverage (0-100) required for pairwise SNP distances to be included in analysis. 
Low-coverage alignments can yield unreliable SNP count.
SNP distance below this threshold are excluded. Default: 80."""
    )
    
    new_snps_parser.add_argument(
        "-t",
        "--threads",
        default=1,
        type=int,
        help = "Thread number. Default is 1."
    )

    new_snps_parser.add_argument(
        "--prefix",
        default=None,
        help = "Prefix for config file. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS"
    )

    new_snps_parser.add_argument(
        "--conda_prefix",
        default=None,
        help = "Directory for conda environments needed for this analysis. If not provided, defaults to <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>"
    )
# The parser for the new-full mode
def _add_new_full_parser(subparsers):
    """Add new-full subparser"""
    
    new_full_parser = subparsers.add_parser(
        "new-full",
        formatter_class=argparse.RawTextHelpFormatter,
        help="Add new genomes and re-derive the phylothresholds from the combined dataset."
    )

    new_full_parser.add_argument(
        "--original_metadata",
        type=str,
        required=True,
        help="""Path to the input metadata file containing original genomes.
The input file should be a tab-delimited file with 3 or 5 columns.
3 columns: genome_name, genome_accession, genome_path
5 columns: genome_name, genome_accession, genome_path, patient_id, collection_date"""
    )

    new_full_parser.add_argument(
        "--new_metadata",
        type=str,
        required=True,
        help="""Path to the input metadata file containing new genomes.
At least 1 new genome should be provided to perform new-snps analysis.
The description of the file format is the same as that for --original_metadata."""
    )

    new_full_parser.add_argument(
        "--thresher_output",
        type=str,
        required=True,
        help="""Path to the existing THRESHER directory.
The existing analysis directory should contain the previous transmission cluster results.
(epi_mode should be 'True' in the prior run)"""
    )

    new_full_parser.add_argument(
        "--output",
        type=str,
        required=False,
        help = """Path to output directory of new-full analysis.
If not provided, defaults to thresher_new_full_<YYYY_MM_DD_HHMMSS> under the current working directory."""
        )
    
    new_full_parser.add_argument(
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

    new_full_parser.add_argument(
        "--whatsgnu_db_path",
        type=str,
        required=False,
        help="""The path to the existing WhatsGNU database.
If provided, the WhatsGNU database will not be downloaded.
If not provided, defaults to <OUTPUT>/whatsgnu/db""",
    )

    new_full_parser.add_argument(
        "--bakta_db_type",
        required=False,
        default="full",
        help="""The type of Bakta database.
Available options: [full, light]
Default is full"""
    )

    new_full_parser.add_argument(
        "--bakta_db_path",
        required=False,
        type=str,
        help="""The path of the directory where the existing Bakta database locates.
If provided, the Bakta database will not be downloaded.
If not provided, defaults to <OUTPUT>/bakta/db""",
    )

    new_full_parser.add_argument(
        "--snp_coverage_threshold",
        required=False,
        type=int,
        default=80,
        help="""Minimum alignment coverage (0-100) required for pairwise SNP distances to be included in analysis. 
Low-coverage alignments can yield unreliable SNP count.
SNP distance below this threshold are excluded. Default: 80."""
    )

    new_full_parser.add_argument(
        "--core_threshold",
        required=False,
        type=float,
        default=0.95,
        help="""Panaroo Core-genome sample threshold. The frequency of a gene in your sample required to classify it as 'core'.
Range is 0.0 to 1.0. Default is 0.95."""
    )

    new_full_parser.add_argument(
        "--core_bootstrap_method",
        required=False,
        default="ultrafast",
        choices=["ultrafast", "nonparametric"],
        help="""The bootstrap method for core genome phylogeny used for hierarchical clustering.
Available options: [ultrafast, nonparametric]
ultrafast: Use ultrafast bootstrap method.
nonparametric: Use standard nonparametric bootstrap method.
Default is ultrafast"""
    )

    new_full_parser.add_argument(
        "--core_bootstrap_number",
        required=False,
        type=int,
        default=0,
        help="""The number of bootstrap number for core genome phylogeny.
If method is ultrafast, default is 1000.
If method is nonparametric, default is 100."""
    )

    new_full_parser.add_argument(
        "--group_bootstrap_method",
        required=False,
        default="ultrafast",
        choices=["ultrafast", "nonparametric"],
        help="""The bootstrap method for phylogeny of each hierarchical group.
Available options: [ultrafast, nonparametric]
ultrafast: Use ultrafast bootstrap method.
nonparametric: Use standard nonparametric bootstrap method.
Default is ultrafast"""
    )

    new_full_parser.add_argument(
        "--group_bootstrap_number",
        required=False,
        type=int,
        default=0,
        help="""The number of bootstrap replicates for phylogeny of each hierarchical group.
If method is ultrafast, default is 1000.
If method is nonparametric, default is 100."""
    )
    
    new_full_parser.add_argument(
        "--threshold_floor",
        required=False,
        type=int,
        default=5,
        help="""The floor of the range of SNP distances to search for the optimal phylothreshold.
Default is 5."""
    )

    new_full_parser.add_argument(
        "--threshold_ceiling",
        required=False,
        type=int,
        default=500,
        help="""The ceiling of the range of SNP distances to search for the optimal phylothreshold.
Default is 500."""
    )

    new_full_parser.add_argument(
        "--singleton_threshold",
        required=False,
        type=int,
        default=100,
        help="""The SNP distance threshold above which every genome in the group is considered a singleton.
Default is 100."""
    )

    new_full_parser.add_argument(
        "--correction_bootstrap",
        required=False,
        type=int,
        default=0,
        help="Minimum bootstrap support threshold for applying phylogenetic corrections to strain composition (default: 0)."
    )

    new_full_parser.add_argument(
        "--use_cladebreaker",
        required=False,
        type=bool,
        default=True,
        choices=[True, False],
        help="""Use CladeBreaker to restrain the strain composition.
Options are [True, False].
Default is True.""")

    new_full_parser.add_argument(
        "--endpoint",
        required=False,
        default= "plateau",
        choices=["plateau", "peak", "discrepancy", "public"],
        help="""The endpoint method to use for determing clusters and making plots.
Available Options: [plateau, peak, discrepancy, public]
plateau : Phylothreshold set at a plateau where further increases no longer change the number or composition of strains within the group
peak: Phylothreshold set at the peak number of clones defined within the group.
discrepancy: Phylothreshold set at the point where the discrepancy is minimized within the group.
public: Phylothreshold set at first time a public genome is included in any strain within the group.
Default is plateau."""
    )

    # Always add plateau_length, but only use it when endpoint is "plateau"
    new_full_parser.add_argument(
        "--plateau_length",
        type=int,
        default=15,
        required=False,
        help="""The plateau length for the plateau endpoint method.
Only used when endpoint method is 'plateau'. Default is 15"""
    )

    new_full_parser.add_argument(
        "-t",
        "--threads",
        default=1,
        type=int,
        help = "Thread number. Default is 1."
    )

    new_full_parser.add_argument(
        "--prefix",
        default=None,
        help = "Prefix for config file. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS"
    )

    new_full_parser.add_argument(
        "--conda_prefix",
        default=None,
        help = f"Directory for conda environments needed for this analysis. If not provided, defaults to <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>"
    )

    # Override system compatibility checks (OS and minimum RAM) and run pipeline regardless
    new_full_parser.add_argument(
        "--force",
        action="store_true",
        help="""Bypass system compatibility checks (operating system and available RAM) and force execution of the pipeline.
This may cause instability or failures."""
    )
# The parser for cladebreaker-off mode
def _add_cladebreaker_off_parser(subparsers):
    """Add cladebreaker-off subparser"""
    cladebreaker_off_parser = subparsers.add_parser(
        "cladebreaker-off",
        formatter_class=argparse.RawTextHelpFormatter,
        help="Re-run identification with the CladeBreaker step disabled, allowing public genomes to cluster with study genomes.",
    )

    cladebreaker_off_parser.add_argument(
        "--thresher_output",
        type=str,
        required=True,
        help="""Path to the existing THRESHER directory.
The existing analysis directory should contain the previous analysis results.""")
    
    cladebreaker_off_parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="""Path to output directory.
If not provided, defaults to thresher_cladebreaker_off_<YYYY_MM_DD_HHMMSS> under the current working directory.""")
    
    cladebreaker_off_parser.add_argument(
        "--threshold_floor",
        required=False,
        type=int,
        default=5,
        help="""The floor of the range of SNP distances to search for the optimal phylothreshold.
Default is 5."""
    )

    cladebreaker_off_parser.add_argument(
        "--threshold_ceiling",
        required=False,
        type=int,
        default=500,
        help="""The ceiling of the range of SNP distances to search for the optimal phylothreshold.
Default is 500."""
    )

    cladebreaker_off_parser.add_argument(
        "--singleton_threshold",
        required=False,
        type=int,
        default=100,
        help="""The SNP distance threshold above which every genome in the group is considered a singleton.
Default is 100."""
    )

    cladebreaker_off_parser.add_argument(
        "--plateau_length",
        required=False,
        type=int,
        default=15,
        help="The plateau length for the plateau endpoint method. Default is 15."
    )

    cladebreaker_off_parser.add_argument(
        "--correction_bootstrap",
        required=False,
        type=int,
        default=0,
        help="Minimum bootstrap support threshold for applying phylogenetic corrections to strain composition (default: 0)."
    )

    cladebreaker_off_parser.add_argument(
        "-t",
        "--threads",
        default=1,
        type=int,
        help = "Thread number. Default is 1."
    )

    cladebreaker_off_parser.add_argument(
        "--prefix",
        default=None,
        help = "Prefix for config files. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS"
    )

    cladebreaker_off_parser.add_argument(
        "--conda_prefix",
        default=None,
        help = "Directory for conda environments needed for this analysis. If not provided, defaults to OUTPUT/conda_envs_<YYYY_MM_DD_HHMMSS>"
    )