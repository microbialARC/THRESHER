# This script serves as the wrapper for the workflow
# It creates the config for snakemake and executes the workflow
import argparse
import os
import sys
import pandas as pd

# Create the parser
def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter
    )
    # Add the arguments
    # Argument for the input file
    parser.add_argument(
    "-i",
    "--input",
    required=True,
    help="""Path to the input file, which is a CSV file with at least 3 columns(lite mode), separated by tab, no header.
The first column is the genome name(lite & full).
The second column is the GenBank accession number(lite & full). If the accession number is not available, put "new"(all lowercase) in the column.
The third column is the path to the genome(lite & full).
The forth column is the patient ID (full).
The fifth column is the collection date. Format: yyyy-mm-dd (full).
"""
    )
    # Argument for the output directory
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="""Path to the output directory."""
    )

    # Argument for bacteria species
    # This will be used to determine the WhatsGNU database to use
    parser.add_argument(
        "--species",
        required=True,
        help="""Bacteria species.
Available options: [sau, sepi, cdiff, kp]
sau: Staphylococcus aureus
sepi: Staphylococcus epidermidis
cdiff: Clostridium difficile
kp: Klebsiella pneumoniae"""
    )

    # Argument for the WhatsGNU database if already downloaded
    parser.add_argument(
        "--whatsgnu_db_path",
        required=False,
        help="""The path to the existing WhatsGNU database.
If provided, the WhatsGNU database will not be downloaded"""
    )
    # Argument for full or light batka database
    parser.add_argument(
        "-db",
        "--bakta_db",
        required=False,
        default="full",
        help="""Bakta database.
Available options: [full, light]
Default is full"""
    )

    # Argument if the bakta database is already downloaded and the path to the database
    parser.add_argument(
        "--bakta_db_path",
        required=False,
        help="""The path of the directory where the existing Bakta database locates.
If provided, the Bakta database will not be downloaded"""
    )

    # Argument for the number of threads
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=16,
        help="""Number of threads to use.
Default is 16"""
    )

    # Argument for the RAM(mb)
    parser.add_argument(
        "--memory",
        type=int,
        default=64,
        help="""RAM(Gb) to use.
WhatsGNU database requires loading the database into memory.
All WhatsGNU databases require at least 64Gb to run.
Default is 64Gb"""
    )

    # Argument for whether to make cluster plots and persistence plot
    parser.add_argument(
        "-m",
        "--mode",
        required=False,
        default="lite",
        help="""Whether to make cluster plots and persistence plot.
Availabel Options: [full, lite] 
full: Determine strains, clusters and make plots.
lite: Determine strains without determing clusters and making plots.
Default is lite"""
    )

    # Argument for the endpoint method
    parser.add_argument(
        "-e",
        "--endpoint",
        required=False,
        default= "plateau",
        help="""The endpoint method to use for determing clusters and making plots.
Availabel Options: [plateau, peak, discrepancy, global]
plateau : Threshold set at the plateau where strain compositions are stable in the group.
peak: Threshold set at the peak number of clones defined in the group.
discrepancy: Threshold set at the point where the discrepancy is minimized in the group. 
global: Threshold set at the first time a global genome is included in any strain in the group.
Default is plateau"""
    )

    # The plateau length for the plateau endpoint method
    parser.add_argument(
        "--plateau_length",
        type=int,
        default=15,
        help="""The plateau length for the plateau endpoint method.
Default is 15"""
    )

    # Parse the arguments
    args = parser.parse_args()
    
    if args.input:
        input_df = pd.read_csv(args.input, header=None, sep="\t")
        if input_df.shape[1] == 3:
            input_df[4] = None
            input_df[5] = None
        elif input_df.shape[1] not in [3, 5]:
            raise ValueError("Invalid number of columns in the input dataframe, provide 3(lite) or 5 columns(full)")
        
        input_df.columns = ["genome_name", "genome_accession", "genome_path", "patient_id", "collection_date"]
        if input_df["patient_id"].isnull().any():
            # If the patient ID is not provided, skip the rule of making cluster plots
            args.mode = "lite"
            print("Patient ID is not provided or not complete, skipping the rule of making cluster plots")
        elif args.mode == "full" and input_df["patient_id"].notnull().all():
        # If the patient ID is provided, make the cluster plots
            print("Patient ID is complete, will make cluster plots")

        # If bakta database path is provided, use the provided path, and skip the rule of downloading the database
        # Otherwise download the database
        if args.bakta_db_path:
            print("Using the provided Bakta database")
        elif not args.bakta_db_path and args.bakta_db == "full":
            print("Full bakta database will be downloaded")
        elif not args.bakta_db_path and args.bakta_db == "light":
            print("Light bakta database will be downloaded")
        
        if args.whatsgnu_db_path:
            print("Using the provided WhatsGNU database")
        elif not args.whatsgnu_db_path and args.species == "sau":
            print("Staphylococcus aureus WhatsGNU database will be downloaded")
        elif not args.whatsgnu_db_path and args.species == "cdiff":
            print("Clostridium difficile WhatsGNU database will be downloaded")
        elif not args.whatsgnu_db_path and args.species == "kp":
            print("Klebsiella pneumoniae WhatsGNU database will be downloaded")
        elif not args.whatsgnu_db_path and args.species == "sepi":
            print("Staphylococcus epidermidis WhatsGNU database will be downloaded")

        # If the species is not provided, raise an error
        if not args.species and not args.whatsgnu_db_path:
            raise SystemExit(
                "Please provide the bacteria species"
            )
        # If the endpoint method is not provided, raise an error
        if not args.endpoint:
            raise SystemExit(
                "Please provide the endpoint method"
            ) 
        # If the endpoint method is not one of the available options, raise an error
        if args.endpoint not in ["plateau", "peak", "discrepancy", "global"]:
            raise SystemExit(
                "Invalid endpoint method"
            )
    elif not args.input:
    # If the input file is not provided, raise an error
        raise SystemExit(
            "Please provide the input file"
        )
    return args

def create_config(args):

    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    # Convert input and output path to absolute path
    input_path = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)
    os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","workflow","config"), exist_ok=True)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","workflow","config", "config.yaml"), "w") as f:
        f.write("input: {}\n".format(input_path))
        f.write("output: {}\n".format(output_path))
        f.write("species: {}\n".format(args.species))
        if args.whatsgnu_db_path:
            f.write("whatsgnu_db_path: {}\n".format(args.whatsgnu_db_path))
        else:
            f.write("whatsgnu_db_path: None\n")
        if args.bakta_db_path:
            f.write("bakta_db_path: {}\n".format(args.bakta_db_path))
            f.write("bakta_db: None\n")
        else:
            f.write("bakta_db_path: None\n")
            f.write("bakta_db: {}\n".format(args.bakta_db))
        f.write("threads: {}\n".format(args.threads))
        f.write("endpoint: {}\n".format(args.endpoint))
        f.write("plateau_length: {}\n".format(args.plateau_length))
        f.write("mode: {}\n".format(args.mode))
        f.write("memory: {}\n".format(args.memory))


# Main function
def main():
    args = parse_args()
    create_config(args)
    snakefile_abs_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","workflow", "Snakefile")
    # Execute the workflow using the Snakefile
    # This needs to be fixed in the future
    # Because I haven't figured out how to use the Snakemake API for the newer version of Snakemake
    # https://snakemake-api.readthedocs.io/en/stable/api_reference/snakemake_api.html
    # So for now I will just use the os.system to execute the workflow
    os.system(f"snakemake --snakefile {snakefile_abs_path} --use-conda -f")

if __name__ == "__main__":
    sys.exit(main())
