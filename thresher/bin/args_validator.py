"""Input validation functions for THRESHER"""
import os
import pandas as pd
import time
import re
from thresher.bin.parse_genome_name import parse_genome_name

class ValidationError(Exception):
    """Validation error"""
    pass

# Validate the function argument
def validate_function(args):
    """Validate the function argument"""
    valid_functions = {"strain_identifier", "genome_profiler", "evo_simulator"}
    if args.command not in valid_functions:
        raise ValidationError(f"Function must be one of: {', '.join(valid_functions)}")
# Validate arguments for strain_identifier full-pipeline mode
def validate_strain_identifier_full(args):
    """Validate the arguments used in full-pipeline of strain_identifier function"""
    
    # Check input file
    if not os.path.exists(args.metadata):
        raise ValidationError(f"Input file not found: {args.metadata}")

    input_df = pd.read_csv(args.metadata, header=None, sep="\t")
    
    if input_df.shape[1] not in [3, 5]:
        raise ValidationError("Input file must have 3 or 5 columns")
    if input_df.shape[1] == 3:
        input_df[4] = None
        input_df[5] = None
    input_df.columns = ["genome_name", "accession", "genome_path", "patient_id", "collection_date"]
    # Check if there is space in the genome_name column
    # if so, parse the genome_name
    # Parse genome names with invalid characters
    for row_index, name_entry in enumerate(input_df["genome_name"]):
        sanitized = parse_genome_name(name_entry)
        if sanitized != name_entry:
            print(f"Warning: Genome name '{name_entry}' contains invalid characters, renamed to '{sanitized}'")
            input_df.at[row_index, "genome_name"] = sanitized
    # use the absolute path for genome_path
    input_df["genome_path"] = input_df["genome_path"].apply(lambda x: os.path.abspath(x))
    if input_df.shape[0] < 4:
        raise ValidationError("At least 4 genomes required for analysis")
    if input_df["accession"].isnull().any():
        raise ValidationError("Accession column contains missing values. If you don't have GenBank accession, please use 'new'(all lower case)")
    
    # Recursively check if genome files exist otherwise snakemake will fail right away
    for genome_path in input_df["genome_path"]:
        if not os.path.exists(genome_path):
            raise ValidationError(f"Genome file not found: {genome_path}")
        
    # Ensure no duplicate genome names
    if input_df["genome_name"].duplicated().any():
        raise ValidationError("Duplicate genome names found in the input file. Please ensure all genome names are unique.")
    # Ensure no duplicate genome paths
    if input_df["genome_path"].duplicated().any():
        raise ValidationError("Duplicate genome paths found in the input file. Please ensure all genome paths are unique.")
    
    # Print info about number of genomes and patient IDs if applicable
    print(f"Input file read successfully. Number of genomes: {input_df.shape[0]}")
    if args.analysis_mode == "lite":
        print("Lite mode selected. Clusters will not be determined.")
    elif args.analysis_mode == "full":
        patient_missing = input_df["patient_id"].isnull().any()
        date_missing = input_df["collection_date"].isnull().any()
        # Check if patient_id column is provided if full analysis mode is selected
        # If any patient_id is missing, switch to lite mode
        if patient_missing:
            print("Patient ID is not complete. Switch to lite mode and no clusters will be determined.")
            args.analysis_mode = "lite"
        else:
            unique_patient_count = input_df['patient_id'].nunique(dropna=True)
            print(f"Number of unique patient IDs: {unique_patient_count}")
        
        # Check collection_date column, if any collection_date is missing, switch to lite mode
        if date_missing:
            print("Collection date column contains missing values. Switch to lite mode and no clusters will be determined.")
            args.analysis_mode = "lite"
        else:
            unique_collection_date_count = input_df['collection_date'].nunique(dropna=True)
            print(f"Number of unique collection dates: {unique_collection_date_count}")
        
        if not patient_missing and not date_missing:
            print("Both patient ID and collection date columns are complete. Transmission clusters will be determined.")
    # Check prefix
    if not args.prefix:
        args.prefix = time.strftime("%Y_%m_%d_%H%M%S")
        print(f"Prefix not provided, using default prefix: current timestamp in the format of YYYY_MM_DD_HHMMSS, which is {args.prefix}")
    
    # Check conda prefix
    if not args.conda_prefix:
        print(f"Conda environment path not provided, using default path: {args.output}/conda_envs_{args.prefix}")
        args.conda_prefix = os.path.abspath(f"{args.output}/conda_envs_{args.prefix}")

    # Check output directory
    if not args.output:
        print(f"Output directory not provided, creating a directory named 'thresher_strain_identifier_output_{args.prefix} under current working directory({os.getcwd()}) as output directory")
        args.output = os.path.abspath(os.path.join(os.getcwd(), f"thresher_strain_identifier_output_{args.prefix}"))

    if not os.path.exists(args.output):
        print(f"Output directory {args.output} does not exist, creating it.")
        os.makedirs(args.output)
        args.output = os.path.abspath(args.output)

    # Check species
    if not args.species:
        raise ValidationError("Species must be provided")
    elif args.species not in {"sau", "sepi", "cdiff", "kp"}:
        raise ValidationError(f"Species must be one of: sau, sepi, cdiff, kp")
    
    # Check full or lite mode
    if args.analysis_mode not in {"full", "lite"}:
        print("Unsupported mode, using default mode 'full'")
        args.analysis_mode = "full"

    # Check whatsgnu_db_path
    if not args.whatsgnu_db_path:
        print(f"WhatsGNU database path not provided. WhatsGNU database will be downloaded to {args.output}/whatsgnu/db")
        args.whatsgnu_db_path = "None"
    else:
        # Check if the file exists in the provided whatsgnu_db_path)
        if not os.path.exists(args.whatsgnu_db_path):
            print(f"WhatsGNU database file not found in the provided whatsgnu_db_path: {args.whatsgnu_db_path}")
            print(f"WhatsGNU database will be downloaded to {args.output}/whatsgnu/db")
            args.whatsgnu_db_path = "None"
        else:
            print(f"WhatsGNU database file found at: {args.whatsgnu_db_path}")

    # Check bakta_db_path and bakta_db_type
    if not args.bakta_db_path:
        print(f"Bakta database path not provided. Bakta database will be downloaded to {args.output}/bakta/db")
        args.bakta_db_path = "None"
    elif args.bakta_db_path:
        # Check if bakta.db exists in the provided bakta_db_path
        bakta_db_file = os.path.join(args.bakta_db_path, "bakta.db")
        if not os.path.exists(bakta_db_file):
            print(f"Bakta database file not found in the provided bakta_db_path: {bakta_db_file}")
            print(f"Bakta database will be downloaded to {args.output}/bakta/db")
            args.bakta_db_path = "None"
        else:
            print(f"Bakta database file found at: {bakta_db_file}")

    if not args.bakta_db_type:
        print("Bakta database type not provided, using default 'full'")
        args.bakta_db_type = "full"
    elif args.bakta_db_type not in {"full", "light"}:
        print("Unsupported Bakta database type, using default 'full'")
        args.bakta_db_type = "full"

    # SNP coverage threshold
    if not args.snp_coverage_threshold:
        print("SNP coverage threshold not provided, using default threshold of 80")
        args.snp_coverage_threshold = 80
    elif args.snp_coverage_threshold < 0 or args.snp_coverage_threshold > 100:
        print("SNP coverage threshold must be between 0 and 100, using default threshold of 80")
        args.snp_coverage_threshold = 80
    # Check core gene threshold
    if not args.core_threshold:
        print("Core gene threshold not provided, using default threshold of 0.95")
        args.core_threshold = 0.95
    elif args.core_threshold <= 0.0 or args.core_threshold > 1.0:
        print("Core gene threshold must be between 0.0 and 1.0, using default threshold of 0.95")
        args.core_threshold = 0.95

    # Check core bootstrap method
    if not args.core_bootstrap_method or args.core_bootstrap_method not in {"ultrafast", "nonparametric"}:
        print("Unsupported core bootstrap method, using default method 'ultrafast'")
        args.core_bootstrap_method = "ultrafast"

    # Check core bootstrap number
    if not args.core_bootstrap_number:
        if args.core_bootstrap_method == "ultrafast":
            print("Core bootstrap number not provided, using default 1000 replicates for ultrafast method")
            args.core_bootstrap_number = 1000
        elif args.core_bootstrap_method == "nonparametric":
            print("Core bootstrap number not provided, using default 100 replicates for standard nonparametric method")
            args.core_bootstrap_number = 100
        
    # Check group bootstrap method
    if not args.group_bootstrap_method or args.group_bootstrap_method not in {"ultrafast", "nonparametric"}:
        print("Unsupported group bootstrap method, using default method 'ultrafast'")
        args.group_bootstrap_method = "ultrafast"
    
    # Check group bootstrap number
    if not args.group_bootstrap_number:
        if args.group_bootstrap_method == "ultrafast":
            print("Group bootstrap number not provided, using default 1000 replicates for ultrafast method")
            args.group_bootstrap_number = 1000
        elif args.group_bootstrap_method == "nonparametric":
            print("Group bootstrap number not provided, using default 100 replicates for standard nonparametric method")
            args.group_bootstrap_number = 100
    
    # Check the value of threshold_ceiling
    if type(args.threshold_ceiling) is not int:
        print("Threshold ceiling(--threshold_ceiling) must be a positive integer, using default 500")
        args.threshold_ceiling = 500
    if not args.threshold_ceiling:
        print("Threshold ceiling(--threshold_ceiling) not provided, using default 500")
        args.threshold_ceiling = 500
    elif args.threshold_ceiling < 0:
        print("Threshold ceiling(--threshold_ceiling) must be a positive integer, using default 500")
        args.threshold_ceiling = 500

    # Check the value of singleton_threshold
    if type(args.singleton_threshold) is not int:
        print("Singleton SNP distance threshold(--singleton_threshold) must be a positive integer, using default 100")
        args.singleton_threshold = 100
    if not args.singleton_threshold:
        print("Singleton SNP distance threshold(--singleton_threshold) not provided, using default 100")
        args.singleton_threshold = 100
    elif args.singleton_threshold < 0:
        print("Singleton SNP distance threshold(--singleton_threshold) must be a positive integer, using default 100")
        args.singleton_threshold = 100

    # Check the value of correction_bootstrap
    if type(args.correction_bootstrap) is not int:
        print("Minimum bootstrap support threshold for applying phylogenetic corrections(--correction_bootstrap) must be an integer between 0 and 100, using default 0")
        args.correction_bootstrap = 0
    if not args.correction_bootstrap:
        print("Minimum bootstrap support threshold for applying phylogenetic corrections(--correction_bootstrap) not provided, using default 0")
        args.correction_bootstrap = 0
    elif args.correction_bootstrap < 0 or args.correction_bootstrap > 100:
        print("Minimum bootstrap support threshold for applying phylogenetic corrections(--correction_bootstrap) must be between 0 and 100, using default 0")
        args.correction_bootstrap = 0

    # Check the bool of use_cladebreaker
    if not args.use_cladebreaker:
        print("use_cladebreaker not provided, using default 'True'")
        args.use_cladebreaker = True
    elif args.use_cladebreaker not in {True, False}:
        print("Unsupported option for use_cladebreaker, using default 'True'")
        args.use_cladebreaker = True

    # Check endpoint method
    if args.endpoint not in {"plateau", "peak", "discrepancy", "global"}:
        print("Unsupported endpoint method, using default method 'plateau'")
        args.endpoint = "plateau"
    
    if args.endpoint == "plateau" and not args.plateau_length:
        print("Plateau length not provided, using default length of 15")
        args.plateau_length = 15

    # Check threads, if not provided, use 1 thread
    if not args.threads:
        print("Thread number not provided or invalid, using 1 thread")
        args.threads = 1
    elif args.threads < 1:
        print("Thread number must be a positive integer, using 1 thread")
        args.threads = 1
    elif args.threads > os.cpu_count():
        print(f"Thread number exceeds available threads ({os.cpu_count()}), using 1 thread")
        args.threads = 1
# Validate arguments for strain_identifier redo-endpoint mode
def validate_strain_identifier_redo_endpoint(args):
    """Validate the arguments used in redo-endpoint of strain_identifier function"""
    # Check Original metadata file (Original input file used in full-pipeline mode)
    if not os.path.exists(args.original_metadata):
        raise ValidationError(f"Input file not found: {args.original_metadata}")
    
    original_metadata_df = pd.read_csv(args.original_metadata, header=None, sep="\t")
    original_metadata_df.columns = ["genome_name", "accession", "genome_path", "patient_id", "collection_date"]
    original_metadata_df["genome_path"] = original_metadata_df["genome_path"].apply(lambda x: os.path.abspath(x))
    if original_metadata_df.shape[1] != 5:
        raise ValidationError("Input file must have 5 columns")
    # Check if there is space in the genome_name column
    # if so, parse the genome_name
    # Parse genome names with invalid characters
    for row_index, name_entry in enumerate(original_metadata_df["genome_name"]):
        sanitized = parse_genome_name(name_entry)
        if sanitized != name_entry:
            print(f"Warning: Genome name '{name_entry}' contains invalid characters, renamed to '{sanitized}'")
            original_metadata_df.at[row_index, "genome_name"] = sanitized
    # Check the exisiting THRESHER directory
    if not args.thresher_output or not os.path.exists(args.thresher_output):
        raise ValidationError("Existing THRESHER output directory not found")
    else:
        # Get absolute path
        args.thresher_output =  os.path.abspath(args.thresher_output)
        # Look for peak_strains.RDS, plateau_strains.RDS, global_strains.RDS, and discrepancy_strains.RDS within the existing THRESHER directory
        thresher_output = {
            "peak": os.path.join(args.thresher_output, "thresher", "output", "peak_strains.RDS"),
            "plateau": os.path.join(args.thresher_output, "thresher", "output", "plateau_strains.RDS"),
            "global": os.path.join(args.thresher_output, "thresher", "output", "global_strains.RDS"),
            "discrepancy": os.path.join(args.thresher_output, "thresher", "output", "discrepancy_strains.RDS")
        }
        # Only 4 files all exist can proceed
        if not all(os.path.exists(f) for f in thresher_output.values()):
            missing_files = [f for f in thresher_output.values() if not os.path.exists(f)]
            raise ValidationError(f"Missing THRESHER output files: {', '.join(missing_files)}")
        else:
            print(f"Existing THRESHER output directory found: {args.thresher_output}")
    # Check endpoint method
    # Available options: {plateau, peak, discrepancy, global}
    if args.endpoint not in {"plateau", "peak", "discrepancy", "global"}:
        print("Unsupported endpoint method, using default method 'plateau'")
        args.endpoint = "plateau"
    
    # Check prefix
    if not args.prefix:
        print("Prefix not provided, using default prefix: current timestamp in the format of YYYY_MM_DD_HHMMSS")
        args.prefix = time.strftime("%Y_%m_%d_%H%M%S")
    # Check conda prefix
    if not args.conda_prefix:
        print(f"Conda environment path not provided, using default path: {args.thresher_output}/conda_envs_{args.prefix}")
        args.conda_prefix = os.path.abspath(f"{args.thresher_output}/conda_envs_{args.prefix}")   

    # Check output directory
    if not args.output:
        print(f"Output directory not provided, creating a directory named 'thresher_strain_identifier_redo_endpoint_{args.prefix} under current working directory({os.getcwd()}) as output directory")
        args.output = os.path.abspath(os.path.join(os.getcwd(), f"thresher_strain_identifier_redo_endpoint_{args.prefix}"))
    if not os.path.exists(args.output):
        print(f"Output directory {args.output} does not exist, creating it.")
        os.makedirs(args.output)
        args.output = os.path.abspath(args.output)
# Validate arguments for strain_identifier new-snps mode
def validate_strain_identifier_new_snps(args):
    """Validate the arguments used in new-snps of strain_identifier function"""

    # Check new metadata file
    if not os.path.exists(args.new_metadata):
        raise ValidationError(f"New metadata file not found: {args.new_metadata}")

    new_metadata_df = pd.read_csv(args.new_metadata, header=None, sep="\t")

    if new_metadata_df.shape[1] not in [3, 5]:
        raise ValidationError("New metadata file must have 3 or 5 columns to use new-snps mode")
    elif new_metadata_df.shape[1] == 3:
        new_metadata_df[4] = None
        new_metadata_df[5] = None
    new_metadata_df.columns = ["genome_name", "accession", "genome_path", "patient_id", "collection_date"]
    # Check if there is space in the genome_name column
    # if so, parse the genome_name
    # Parse genome names with invalid characters
    for row_index, name_entry in enumerate(new_metadata_df["genome_name"]):
        sanitized = parse_genome_name(name_entry)
        if sanitized != name_entry:
            print(f"Warning: Genome name '{name_entry}' contains invalid characters, renamed to '{sanitized}'")
            new_metadata_df.at[row_index, "genome_name"] = sanitized
    # Use the absolute path for genome_path
    new_metadata_df["genome_path"] = new_metadata_df["genome_path"].apply(lambda x: os.path.abspath(x))
    if new_metadata_df["accession"].isnull().any():
        raise ValidationError("Accession column contains missing values. If you don't have GenBank accession, please use 'new'(all lower case)")
    
    # Recursively check if genome files exist otherwise snakemake will fail right away
    for genome_path in new_metadata_df["genome_path"]:
        if not os.path.exists(genome_path):
            raise ValidationError(f"Genome file not found: {genome_path}")
        
    # Ensure no duplicate genome names
    if new_metadata_df["genome_name"].duplicated().any():
        raise ValidationError("Duplicate genome names found in the new metadata file. Please ensure all genome names are unique.")
    # Ensure no duplicate genome paths
    if new_metadata_df["genome_path"].duplicated().any():
        raise ValidationError("Duplicate genome paths found in the new metadata file. Please ensure all genome paths are unique.")
    # Print info about number of genomes in the new metadata file
    print(f"New metadata file read successfully. Number of genomes: {new_metadata_df.shape[0]}")

    # Check if transmission clusters can be determined based on if previous THRESHER output contains cluster summary file
    if os.path.exists(os.path.join(args.thresher_output,  "thresher", "output", "clusters_summary.RDS")):
        print("Previous THRESHER output contains cluster summary file.")
        # Check if patient_id is complete
        if not new_metadata_df["patient_id"].isnull().any():
            raise ValidationError("Patient ID column contains missing values.")
        else:
            unique_new_patient_count = new_metadata_df['patient_id'].nunique(dropna=True)
            print(f"Number of unique patient IDs in new metadata file: {unique_new_patient_count}")
        # Check if collection_date is complete
        if not new_metadata_df["collection_date"].isnull().any():
            raise ValidationError("Collection date column contains missing values.")
    else:
        print("Previous THRESHER output does not contain cluster summary file. Transmission clusters will not be determined.")


    # Check original metadata file (Origninal input file used in full-pipeline mode)
    if not os.path.exists(args.original_metadata):
        raise ValidationError(f"Original metadata file not found: {args.original_metadata}")
    original_metadata_df = pd.read_csv(args.original_metadata, header=None, sep="\t")
    if original_metadata_df.shape[1] not in [3, 5]:
        raise ValidationError("Original metadata file must have 3 or 5 columns")
    if original_metadata_df.shape[1] == 3:
        original_metadata_df[4] = None
        original_metadata_df[5] = None
    original_metadata_df.columns = ["genome_name", "accession", "genome_path", "patient_id", "collection_date"]

    # Check if there is overlap between new and original genomes
    overlapping_genome = set(new_metadata_df["genome_name"]).intersection(set(original_metadata_df["genome_name"]))
    if overlapping_genome:
        raise ValidationError(f"Overlapping genome names found between new and original metadata files: {', '.join(overlapping_genome)}")
    
    # Ensure no duplicate genome names
    if original_metadata_df["genome_name"].duplicated().any():
        raise ValidationError("Duplicate genome names found in the original metadata file. Please ensure all genome names are unique.")
    # Ensure no duplicate genome paths
    if original_metadata_df["genome_path"].duplicated().any():
        raise ValidationError("Duplicate genome paths found in the original metadata file. Please ensure all genome paths are unique.")
    
    # Check the existing THRESHER directory
    if not args.thresher_output or not os.path.exists(args.thresher_output):
        raise ValidationError("Existing THRESHER output directory not found")
    else:
        # Get absolute path
        args.thresher_output =  os.path.abspath(args.thresher_output)
        # Look for output files within the existing THRESHER directory
        thresher_mandatory_output = {
            # Prior strain compositions determined by THRESHER
            "peak_strains": os.path.join(args.thresher_output, "thresher", "output", "peak_strains.RDS"),
            "plateau_strains": os.path.join(args.thresher_output, "thresher", "output", "plateau_strains.RDS"),
            "global_strains": os.path.join(args.thresher_output, "thresher", "output", "global_strains.RDS"),
            "discrepancy_strains": os.path.join(args.thresher_output, "thresher", "output", "discrepancy_strains.RDS"),
            # Prior SNP matrix from the study genomes
            "study_snp_matrix": os.path.join(args.thresher_output, "mummer4_study", "study_snp_matrix.RDS")
        }
        
        thresher_optional_output = {
            # Prior transmission cluster compositions 
            "cluster_summary": os.path.join(args.thresher_output, "thresher", "output", "clusters_summary.csv")
        }
        # Only all files all exist can proceed
        if not all(os.path.exists(f) for f in thresher_mandatory_output.values()):
            missing_files = [f for f in thresher_mandatory_output.values() if not os.path.exists(f)]
            raise ValidationError(f"Missing THRESHER output files: {', '.join(missing_files)}")
        else:
            print(f"Existing THRESHER output directory found: {args.thresher_output}")
        
        if not all(os.path.exists(f) for f in thresher_optional_output.values()):
            print("Warning: Missing prior transmission cluster compositions. " \
            "Clusters will not be analyzed or reported for the new genomes.")
        else:
            print("Prior transmission cluster compositions found.")

    # SNP coverage threshold
    if not args.snp_coverage_threshold:
        print("SNP coverage threshold not provided, using default threshold of 80")
        args.snp_coverage_threshold = 80
    elif args.snp_coverage_threshold < 0 or args.snp_coverage_threshold > 100:
        print("SNP coverage threshold must be between 0 and 100, using default threshold of 80")
        args.snp_coverage_threshold = 80
    
    # Check species
    if not args.species:
        raise ValidationError("Species must be provided")
    elif args.species not in {"sau", "sepi", "cdiff", "kp"}:
        raise ValidationError(f"Species must be one of: sau, sepi, cdiff, kp")
    
    # Check threads, if not provided, use 1 thread
    if not args.threads:
        print("Thread number not provided or invalid, using 1 thread")
        args.threads = 1
    elif args.threads < 1:
        print("Thread number must be a positive integer, using 1 thread")
        args.threads = 1
    elif args.threads > os.cpu_count():
        print(f"Thread number exceeds available threads ({os.cpu_count()}), using 1 thread")
        args.threads = 1
        
    # Check prefix
    if not args.prefix:
        print("Prefix not provided, using default prefix: current timestamp in the format of YYYY_MM_DD_HHMMSS")
        args.prefix = time.strftime("%Y_%m_%d_%H%M%S")
    
    # Check conda prefix
    if not args.conda_prefix:
        print("Conda environment path not provided, using default path: <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>")
        args.conda_prefix = os.path.abspath(f"{args.output}/conda_envs_{args.prefix}")

    # Check output directory
    if not args.output:
        print(f"Output directory not provided, creating a directory named 'thresher_strain_identifier_new_snps_{args.prefix} under current working directory({os.getcwd()}) as output directory")
        args.output = os.path.abspath(os.path.join(os.getcwd(), f"thresher_strain_identifier_new_snps_{args.prefix}"))
    if not os.path.exists(args.output):
        print(f"Output directory {args.output} does not exist, creating it.")
        os.makedirs(args.output)
        args.output = os.path.abspath(args.output)
# Validate arguments for clinality new-full mode
def validate_strain_identifier_new_full(args):
    """Validate the arguments used in new-full of strain_identifier function"""
    
    # Check new metadata file
    if not os.path.exists(args.new_metadata):
        raise ValidationError(f"New metadata file not found: {args.new_metadata}")

    new_metadata_df = pd.read_csv(args.new_metadata, header=None, sep="\t")

    if new_metadata_df.shape[1] not in [3, 5]:
        raise ValidationError("New metadata file must have 3 or 5 columns to use new-full mode")
    elif new_metadata_df.shape[1] == 3:
        new_metadata_df[4] = None
        new_metadata_df[5] = None
    new_metadata_df.columns = ["genome_name", "accession", "genome_path", "patient_id", "collection_date"]
    # Check if there is space in the genome_name column
    # if so, parse the genome_name
    # Parse genome names with invalid characters
    for row_index, name_entry in enumerate(new_metadata_df["genome_name"]):
        sanitized = parse_genome_name(name_entry)
        if sanitized != name_entry:
            print(f"Warning: Genome name '{name_entry}' contains invalid characters, renamed to '{sanitized}'")
            new_metadata_df.at[row_index, "genome_name"] = sanitized
    # Use the absolute path for genome_path
    new_metadata_df["genome_path"] = new_metadata_df["genome_path"].apply(lambda x: os.path.abspath(x))
    if new_metadata_df["accession"].isnull().any():
        raise ValidationError("Accession column contains missing values. If you don't have GenBank accession, please use 'new'(all lower case)")
    
    # Recursively check if genome files exist otherwise snakemake will fail right away
    for genome_path in new_metadata_df["genome_path"]:
        if not os.path.exists(genome_path):
            raise ValidationError(f"Genome file not found: {genome_path}")
    # Print info about number of genomes in the new metadata file
    print(f"New metadata file read successfully. Number of genomes: {new_metadata_df.shape[0]}")
    
    # Ensure no duplicate genome names
    if new_metadata_df["genome_name"].duplicated().any():
        raise ValidationError("Duplicate genome names found in the new metadata file. Please ensure all genome names are unique.")
    # Ensure no duplicate genome paths
    if new_metadata_df["genome_path"].duplicated().any():
        raise ValidationError("Duplicate genome paths found in the new metadata file. Please ensure all genome paths are unique.")
    # Check if transmission clusters can be determined based on if previous THRESHER output contains cluster summary file
    if os.path.exists(os.path.join(args.thresher_output,  "thresher", "output", "clusters_summary.RDS")):
        print("Previous THRESHER output contains cluster summary file.")
        # Check if patient_id is complete
        if not new_metadata_df["patient_id"].isnull().any():
            raise ValidationError("Patient ID column contains missing values.")
        else:
            unique_new_patient_count = new_metadata_df['patient_id'].nunique(dropna=True)
            print(f"Number of unique patient IDs in new metadata file: {unique_new_patient_count}")
        # Check if collection_date is complete
        if not new_metadata_df["collection_date"].isnull().any():
            raise ValidationError("Collection date column contains missing values.")
    else:
        print("Previous THRESHER output does not contain cluster summary file. Transmission clusters will not be determined.")

    # Check original metadata file (Origninal input file used in full-pipeline mode)
    if not os.path.exists(args.original_metadata):
        raise ValidationError(f"Original metadata file not found: {args.original_metadata}")
    original_metadata_df = pd.read_csv(args.original_metadata, header=None, sep="\t")
    if original_metadata_df.shape[1] not in [3, 5]:
        raise ValidationError("Original metadata file must have 3 or 5 columns")
    if original_metadata_df.shape[1] == 3:
        original_metadata_df[4] = None
        original_metadata_df[5] = None
    original_metadata_df.columns = ["genome_name", "accession", "genome_path", "patient_id", "collection_date"]
    #  Check if there is space in the genome_name column
    # if so, parse the genome_name
    # Parse genome names with invalid characters
    for row_index, name_entry in enumerate(original_metadata_df["genome_name"]):
        sanitized = parse_genome_name(name_entry)
        if sanitized != name_entry:
            print(f"Warning: Genome name '{name_entry}' contains invalid characters, renamed to '{sanitized}'")
            original_metadata_df.at[row_index, "genome_name"] = sanitized
    # Create another dataframe to keep those genomes without accession number, which are those not publicly available in GenBank
    filtered_original_metadata_df = original_metadata_df[original_metadata_df["accession"] == "new"]
    
    # Check if there is overlap between new and original genomes
    overlapping_genome = set(new_metadata_df["genome_name"]).intersection(set(original_metadata_df["genome_name"]))
    if overlapping_genome:
        raise ValidationError(f"Overlapping genome names found between new and original metadata files: {', '.join(overlapping_genome)}")
    
    # Check the existing THRESHER directory
    if not args.thresher_output or not os.path.exists(args.thresher_output):
        raise ValidationError("Existing THRESHER output directory not found")
    else:
        # Get absolute path
        args.thresher_output =  os.path.abspath(args.thresher_output)
        # Look for output files within the existing THRESHER directory
        thresher_mandatory_output = {
            # Prior whatsGNU top genome output
            "whatsgnu_results": {genome: os.path.join(args.thresher_output, "whatsgnu", genome, f"{genome}_WhatsGNU_topgenomes.txt") for genome in filtered_original_metadata_df["genome_name"]},
            # Prior strain compositions determined by THRESHER
            "peak_strains": os.path.join(args.thresher_output, "thresher", "output", "peak_strains.RDS"),
            "plateau_strains": os.path.join(args.thresher_output, "thresher", "output", "plateau_strains.RDS"),
            "global_strains": os.path.join(args.thresher_output, "thresher", "output", "global_strains.RDS"),
            "discrepancy_strains": os.path.join(args.thresher_output, "thresher", "output", "discrepancy_strains.RDS"),
            # Prior SNP matrix from the study genomes
            "study_snp_matrix": os.path.join(args.thresher_output, "mummer4_study", "study_snp_matrix.RDS"),
            # Prior SNP matrix from the global genomes
            "global_snp_matrix": os.path.join(args.thresher_output, "mummer4_global", "global_snp_matrix.RDS"),
            # Prior panaroo result directory.
            # Verify panaroo output by checking for core_gene_alignment_filtered.aln;
            # if this file exists, other panaroo outputs are expected to be present.
            "panaroo_results": os.path.join(args.thresher_output, "panaroo","core_gene_alignment_filtered.aln"),
            # Prior bakta annotation results
            "bakta_annotations": {genome: os.path.join(args.thresher_output, "bakta_annotation", genome, f"{genome}.gff3") for genome in original_metadata_df["genome_name"]}
        }
        
        thresher_optional_output = {
            # Prior transmission cluster compositions 
            "cluster_summary": os.path.join(args.thresher_output, "thresher", "output", "clusters_summary.csv")
        }
        # Only all files exist can proceed
        # Because whatsGNU_results and bakta_annotations are dictionaries, we need to flatten them first
        
        thresher_optional_output_flat = [file for value_entry in thresher_optional_output.values() for file in (value_entry.values() if isinstance(value_entry, dict) else [value_entry])]
        thresher_mandatory_output_flat = [file for value_entry in thresher_mandatory_output.values() for file in (value_entry.values() if isinstance(value_entry, dict) else [value_entry])]

        if not all(os.path.exists(file) for file in thresher_mandatory_output_flat):
            missing_files = [file for file in thresher_mandatory_output_flat if not os.path.exists(file)]
            raise ValidationError(f"Missing THRESHER output files: {', '.join(missing_files)}")
        else:
            print(f"Existing THRESHER output directory found: {args.thresher_output}")

        if not all(os.path.exists(file) for file in thresher_optional_output_flat):
            print("Warning: Missing prior transmission cluster compositions. " \
            "Clusters will not be analyzed or reported for the new genomes.")
        else:
            print("Prior transmission cluster compositions found.")
    # Check prefix
    if not args.prefix:
        args.prefix = time.strftime("%Y_%m_%d_%H%M%S")
        print(f"Prefix not provided, using default prefix: current timestamp in the format of YYYY_MM_DD_HHMMSS, which is {args.prefix}")
    
    # Check conda prefix
    if not args.conda_prefix:
        print(f"Conda environment path not provided, using default path: {args.output}/conda_envs_{args.prefix}")
        args.conda_prefix = os.path.abspath(f"{args.output}/conda_envs_{args.prefix}")
    
    # Check output directory
    if not args.output:
        print(f"Output directory not provided, creating a directory named 'thresher_strain_identifier_new_full_{args.prefix}' under current working directory({os.getcwd()}) as output directory")
        args.output = os.path.abspath(os.path.join(os.getcwd(), f"thresher_strain_identifier_new_full_{args.prefix}"))

    if not os.path.exists(args.output):
        print(f"Output directory {args.output} does not exist, creating it.")
        os.makedirs(args.output)
        args.output = os.path.abspath(args.output)

    # Check species
    if not args.species:
        raise ValidationError("Species must be provided")
    elif args.species not in {"sau", "sepi", "cdiff", "kp"}:
        raise ValidationError(f"Species must be one of: sau, sepi, cdiff, kp")
    
    # Check whatsgnu_db_path
    if not args.whatsgnu_db_path:
        print(f"WhatsGNU database path not provided. WhatsGNU database will be downloaded to {args.output}/whatsgnu/db")
    else:
        # Check if the file exists in the provided whatsgnu_db_path)
        if not os.path.exists(args.whatsgnu_db_path):
            print(f"WhatsGNU database file not found in the provided whatsgnu_db_path: {args.whatsgnu_db_path}")
            print(f"WhatsGNU database will be downloaded to {args.output}/whatsgnu/db")
            args.whatsgnu_db_path = "None"
        else:
            print(f"WhatsGNU database file found at: {args.whatsgnu_db_path}")

    # Check bakta_db_path and bakta_db_type
    if not args.bakta_db_path:
        print(f"Bakta database path not provided. Bakta database will be downloaded to {args.output}/bakta/db")
        if not args.bakta_db_type:
            print("Bakta database type not provided, using default 'full'")
            args.bakta_db_type = "full"
    elif args.bakta_db_path:
        # Check if bakta.db exists in the provided bakta_db_path
        bakta_db_file = os.path.join(args.bakta_db_path, "bakta.db")
        if not os.path.exists(bakta_db_file):
            print(f"Bakta database file not found in the provided bakta_db_path: {bakta_db_file}")
            print(f"Bakta database will be downloaded to Bakta database will be downloaded to {args.output}/bakta/db")
            args.bakta_db_path = "None"
        else:
            print(f"Bakta database file found at: {bakta_db_file}")

    if args.bakta_db_type not in {"full", "light"}:
        print("Unsupported Bakta database type, using default 'full'")
        args.bakta_db_type = "full"
    
    # Check the existing THRESHER directory
    if not args.thresher_output or not os.path.exists(args.thresher_output):
        raise ValidationError("Existing THRESHER output directory not found")
    else:
        # Get absolute path
        args.thresher_output =  os.path.abspath(args.thresher_output)
        # Look for peak_strains.RDS, plateau_strains.RDS, global_strains.RDS, and discrepancy_strains.RDS within the existing THRESHER directory
        thresher_output = {
            "peak": os.path.join(args.thresher_output, "thresher", "output", "peak_strains.RDS"),
            "plateau": os.path.join(args.thresher_output, "thresher", "output", "plateau_strains.RDS"),
            "global": os.path.join(args.thresher_output, "thresher", "output", "global_strains.RDS"),
            "discrepancy": os.path.join(args.thresher_output, "thresher", "output", "discrepancy_strains.RDS"),
        }
        # Only 4 files all exist can proceed
        if not all(os.path.exists(f) for f in thresher_output.values()):
            missing_files = [f for f in thresher_output.values() if not os.path.exists(f)]
            raise ValidationError(f"Missing THRESHER output files: {', '.join(missing_files)}")
        else:
            print(f"Existing THRESHER output directory found: {args.thresher_output}")
    
    # Check threads, if not provided, use 1 thread
    if not args.threads:
        print("Thread number not provided or invalid, using 1 thread")
        args.threads = 1
    elif args.threads < 1:
        print("Thread number must be a positive integer, using 1 thread")
        args.threads = 1
    elif args.threads > os.cpu_count():
        print(f"Thread number exceeds available threads ({os.cpu_count()}), using 1 thread")
        args.threads = 1

    # SNP coverage threshold
    if not args.snp_coverage_threshold:
        print("SNP coverage threshold not provided, using default threshold of 80")
        args.snp_coverage_threshold = 80
    elif args.snp_coverage_threshold < 0 or args.snp_coverage_threshold > 100:
        print("SNP coverage threshold must be between 0 and 100, using default threshold of 80")
        args.snp_coverage_threshold = 80

    # Check core gene threshold
    if not args.core_threshold:
        print("Core gene threshold not provided, using default threshold of 0.95")
        args.core_threshold = 0.95
    elif args.core_threshold <= 0.0 or args.core_threshold > 1.0:
        print("Core gene threshold must be between 0.0 and 1.0, using default threshold of 0.95")
        args.core_threshold = 0.95
    
    # Check core bootstrap method
    if not args.core_bootstrap_method or args.core_bootstrap_method not in {"ultrafast", "nonparametric"}:
        print("Unsupported core bootstrap method, using default method 'ultrafast'")
        args.core_bootstrap_method = "ultrafast"

    # Check core bootstrap number
    if not args.core_bootstrap_number:
        if args.core_bootstrap_method == "ultrafast":
            print("Core bootstrap number not provided, using default 1000 replicates for ultrafast method")
            args.core_bootstrap_number = 1000
        elif args.core_bootstrap_method == "nonparametric":
            print("Core bootstrap number not provided, using default 100 replicates for standard nonparametric method")
            args.core_bootstrap_number = 100
    
    # Check group bootstrap method
    if not args.group_bootstrap_method or args.group_bootstrap_method not in {"ultrafast", "nonparametric"}:
        print("Unsupported group bootstrap method, using default method 'ultrafast'")
        args.group_bootstrap_method = "ultrafast"
    
    # Check group bootstrap number
    if not args.group_bootstrap_number:
        if args.group_bootstrap_method == "ultrafast":
            print("Group bootstrap number not provided, using default 1000 replicates for ultrafast method")
            args.group_bootstrap_number = 1000
        elif args.group_bootstrap_method == "nonparametric":
            print("Group bootstrap number not provided, using default 100 replicates for standard nonparametric method")
            args.group_bootstrap_number = 100
    
    # Check the value of threshold_ceiling
    if type(args.threshold_ceiling) is not int:
        print("Threshold ceiling(--threshold_ceiling) must be a positive integer, using default 500")
        args.threshold_ceiling = 500
    if not args.threshold_ceiling:
        print("Threshold ceiling(--threshold_ceiling) not provided, using default 500")
        args.threshold_ceiling = 500
    elif args.threshold_ceiling < 0:
        print("Threshold ceiling(--threshold_ceiling) must be a positive integer, using default 500")
        args.threshold_ceiling = 500
    
    # Check the value of singleton_threshold
    if type(args.singleton_threshold) is not int:
        print("Singleton SNP distance threshold(--singleton_threshold) must be a positive integer, using default 100")
        args.singleton_threshold = 100
    if not args.singleton_threshold:
        print("Singleton SNP distance threshold(--singleton_threshold) not provided, using default 100")
        args.singleton_threshold = 100
    elif args.singleton_threshold < 0:
        print("Singleton SNP distance threshold(--singleton_threshold) must be a positive integer, using default 100")
        args.singleton_threshold = 100
    
    if type(args.correction_bootstrap) is not int:
        print("Minimum bootstrap support threshold for applying phylogenetic corrections(--correction_bootstrap) must be an integer between 0 and 100, using default 0")
        args.correction_bootstrap = 0
    if not args.correction_bootstrap:
        print("Minimum bootstrap support threshold for applying phylogenetic corrections(--correction_bootstrap) not provided, using default 0")
        args.correction_bootstrap = 0
    elif args.correction_bootstrap < 0 or args.correction_bootstrap > 100:
        print("Minimum bootstrap support threshold for applying phylogenetic corrections(--correction_bootstrap) must be between 0 and 100, using default 0")
        args.correction_bootstrap = 0
    
    # Check the bool of use_cladebreaker
    if not args.use_cladebreaker:
        print("use_cladebreaker not provided, using default 'True'")
        args.use_cladebreaker = True
    elif args.use_cladebreaker not in {True, False}:
        print("Unsupported option for use_cladebreaker, using default 'True'")
        args.use_cladebreaker = True

    # check endpoint method
    # available options: {plateau, peak, discrepancy, global}
    if args.endpoint not in {"plateau", "peak", "discrepancy", "global"}:
        print("Unsupported endpoint method, using default method 'plateau'")
        args.endpoint = "plateau"

    if args.endpoint == "plateau" and not args.plateau_length:
        print("Plateau length not provided, using default length of 15")
        args.plateau_length = 15
# Validate arguments for genome_profiler function
def validate_genome_profiler(args):
    """Validate the arguments used in genome_profiler function"""
    # Check if input_genome exists
    # For now I don't check if the input_genome is a valid fasta file for simplicity
    # But this will be checked during the actual analysis
    if not os.path.exists(args.input_genome):
        raise SystemExit(
            f"Input genome {args.input_genome} does not exist"
        )
    
    # Check prefix
    if not args.prefix:
        print("Prefix not provided, using default prefix: current timestamp in the format of YYYY_MM_DD_HHMMSS")
        args.prefix = time.strftime("%Y_%m_%d_%H%M%S")

    # Check output directory
    if not args.output:
        print(f"Output directory not provided, creating a directory named thresher_genome_profiler_{args.prefix} under current working directory({os.getcwd()}) as output directory")
        args.output = os.path.join(os.getcwd(), f"thresher_genome_profiler_{args.prefix}")
    if not os.path.exists(args.output):
        print(f"Output directory {args.output} does not exist, creating it.")
        os.makedirs(args.output)
        args.output = os.path.abspath(args.output)
    else:
        args.output = os.path.abspath(args.output)
    # Check top_genomes
    if args.top_genomes < 100:
        raise ValidationError("top_genomes must be at least 100 for reliable profiling")
    
    # Check ani_threshold
    if not args.ani_threshold:
        print("ani_threshold not provided, using default 99.5")
        args.ani_threshold = 99.5
    elif args.ani_threshold < 95 or args.ani_threshold > 100:
        print("Warning: ani_threshold must be between 95 and 100 for reliable profiling.\n"
        "Using default threshold 99.5.")
        args.ani_threshold = 99.5
    
    # Check conda prefix
    if not args.conda_prefix:
        print("Conda environment path not provided, using default path: <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>")
        args.conda_prefix = os.path.abspath(f"{args.output}/conda_envs_{args.prefix}")
    
    # Check whatsgnu_db_path
    if not args.whatsgnu_db_path:
        print(f"WhatsGNU database path not provided. WhatsGNU database will be downloaded to {args.output}/whatsgnu/db")
        args.whatsgnu_db_path = "None"
    else:
        # Check if the file exists in the provided whatsgnu_db_path)
        if not os.path.exists(args.whatsgnu_db_path):
            print(f"WhatsGNU database file not found in the provided whatsgnu_db_path: {args.whatsgnu_db_path}")
            print(f"WhatsGNU database will be downloaded to {args.output}/whatsgnu/db")
            args.whatsgnu_db_path = "None"
        else:
            print(f"WhatsGNU database file found at: {args.whatsgnu_db_path}")

    # Check bakta_db_path and bakta_db_type
    if not args.bakta_db_path:
        print(f"Bakta database path not provided. Bakta database will be downloaded to {args.output}/bakta/db")
        args.bakta_db_path = "None"
    elif args.bakta_db_path:
        # Check if bakta.db exists in the provided bakta_db_path
        bakta_db_file = os.path.join(args.bakta_db_path, "bakta.db")
        if not os.path.exists(bakta_db_file):
            print(f"Bakta database file not found in the provided bakta_db_path: {bakta_db_file}")
            print(f"Bakta database will be downloaded to {args.output}/bakta/db")
            args.bakta_db_path = "None"
        else:
            print(f"Bakta database file found at: {bakta_db_file}")

    if not args.bakta_db_type:
        print("Bakta database type not provided, using default 'full'")
        args.bakta_db_type = "full"
    elif args.bakta_db_type not in {"full", "light"}:
        print("Unsupported Bakta database type, using default 'full'")
        args.bakta_db_type = "full"
# Validate the shared arguments for evo_simulator function
def validate_evo_simulator(args):
    """Validate the shared arguments used in evo_simulator function"""
    # Check prefix
    if not args.prefix:
        print("Prefix not provided, using default prefix: current timestamp in the format of YYYY_MM_DD_HHMMSS")
        args.prefix = time.strftime("%Y_%m_%d_%H%M%S")
    
    # Check conda prefix
    if not args.conda_prefix:
        print("Conda environment path not provided, using default path: <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>")
        args.conda_prefix = os.path.abspath(f"{args.output}/conda_envs_{args.prefix}")
    
    # Check seed
    if not args.seed:
        print("Seed not provided, using current date in the format of YYYYMMDD as seed")
        args.seed = int(pd.Timestamp.now().strftime("%Y%m%d"))
    elif args.seed < 0 or not isinstance(args.seed, int):
        raise ValidationError("Seed must be a non-negative integer")
    
    # Check output directory
    if not args.output:
        print(f"Output directory not provided, creating a directory named 'thresher_evo_simulator_output_{args.prefix} under current working directory({os.getcwd()}) as output directory")
        args.output = os.path.join(os.getcwd(), f"thresher_evo_simulator_output_{args.prefix}")
    if not os.path.exists(args.output):
        print(f"Output directory {args.output} does not exist, creating it.")
        os.makedirs(args.output)
        args.output = os.path.abspath(args.output)
    
    # Check ancestor_genome
    if not os.path.exists(args.ancestor):
        raise ValidationError(f"Ancestor genome file not found: {args.ancestor}")
    else:
        # Use absolute path
        args.ancestor = os.path.abspath(args.ancestor)
    
    # Check years
    if not args.years :
        raise ValidationError("Years must be provided")
    elif args.years <= 0 or not isinstance(args.years, int):
        raise ValidationError("Years must be a positive integer")

    if not args.taxa:
        raise ValidationError("Taxa (offspring) must be provided")
    elif args.taxa < 3 or not isinstance(args.taxa, int):
        raise ValidationError("Taxa (offspring) must be an integer no less than 3")

    # Check substitution model 
    if not args.substitution_model:
        print("Substitution model not provided, using default 'GTR'")
        args.substitution_model = "GTR"
    elif args.substitution_model not in {"JC69","K2P","K3P","GTR"}:
        print("Unsupported substitution model, using default 'GTR'")
        args.substitution_model = "GTR"

    # Check substitution model parameters for different models
    if args.substitution_model == "JC69":
        print("JC69 model selected, no additional parameters required.")
        args.model_parameters = "None"
    elif args.substitution_model == "K2P":
        print("K2P model selected, validating parameter...")
        if not args.model_parameters:
            raise ValidationError("K2P model requires transition/transversion ratio(kappa) as parameter. The parameter should be provided as a single float value greater than 0.")
        elif isinstance(float(args.model_parameters), float) or isinstance(int(args.model_parameters), int):
            if args.model_parameters <= 0:
                raise ValidationError("Transition/transversion ratio(kappa) must be a float value greater than 0 for K2P model")
            else:
                print(f"K2P model parameter (transition/transversion ratio, kappa) set to {args.model_parameters}")
        else:
            raise ValidationError("Transition/transversion ratio must be a float value greater than 0 for K2P model")
    elif args.substitution_model == "K3P":
        print("K3P model selected, validating parameters...")
        if not args.model_parameters:
            raise ValidationError("K3P model requires three parameters (alpha, gamma, beta) as parameters. The parameters should be provided as three float values greater than 0, separated by a comma (e.g., 0.4,0.6,0.8).")
        else:
            splited_params = args.model_parameters.split(",")
            if len(splited_params) != 3:
                raise ValidationError("K3P model requires three parameters (alpha, beta, gamma) as parameters. The parameters should be provided as three float values greater than 0, separated by commas (e.g., 0.4,0.6,0.8).")
            try:
                alpha = float(splited_params[0])
                beta = float(splited_params[1])
                gamma = float(splited_params[2])
                if alpha <= 0 or gamma <= 0 or beta <= 0:
                    raise ValidationError("All alpha, beta, and gamma parameters must be float values greater than 0 for K3P model")
                else:
                    print(f"K3P model parameters set to alpha: {alpha}, beta: {beta}, gamma: {gamma}")
            except ValueError:
                raise ValidationError("All alpha, beta, and gamma parameters must be float values greater than 0 for K3P model")
    elif args.substitution_model == "GTR":
        print("GTR model selected, validating parameters...")
        if not args.model_parameters:
            raise ValidationError("GTR model requires six rate parameters (a, b, c, d, e, f) as parameters. The parameters should be provided as six float values greater than 0, separated by commas (e.g., 0.1,0.2,0.3,0.15,0.15,0.1).")
        else:
            splited_params = args.model_parameters.split(",")
            if len(splited_params) != 6:
                raise ValidationError("GTR model requires six rate parameters (a, b, c, d, e, f) as parameters. The parameters should be provided as six float values greater than 0, separated by commas (e.g., 0.1,0.2,0.3,0.15,0.15,0.1).")
            try:
                rates = [float(param) for param in splited_params]
                if any(rate <= 0 for rate in rates):
                    raise ValidationError("All six rate parameters must be float values greater than 0 for GTR model")
                else:
                    print(f"GTR model parameters set to: a={rates[0]}, b={rates[1]}, c={rates[2]}, d={rates[3]}, e={rates[4]}, f={rates[5]}")
            except ValueError:
                raise ValidationError("All six rate parameters must be float values greater than 0 for GTR model")
        
    # Check mutation rate
    if not args.mutation_rate or not isinstance(args.mutation_rate, (int, float)) or args.mutation_rate <= 0:
        raise ValidationError("Mutation rate must be a positive float value")
    
    # Check preset
    # If preset is provided, all rates and mge profiles will be overwritten by the preset values
    # This check is temporary unavailable until all presets are well established
    # Will be implemented in future versions before Q1 2026

    # If preset is provided, the minimal required arguments will be taxa and years
    # Whereas other parameters will be predefined
    if args.preset:
        print("Preset option is currently unavailable. Please provide all required parameters manually.")
        args.preset = "None"

    # Check if weighted mutation is on
    if not args.use_weighted_mutation:
        print("Whether using weighted mutation not provided, using default True")
        args.use_weighted_mutation = True
    elif args.use_weighted_mutation not in {True, False}:
        print("use_weighted_mutation must be either True or False. Using default True")
        args.use_weighted_mutation = True
    elif args.use_weighted_mutation:
        print("Weighted mutation enabled, validating weighted mutation input...")
        # Check weighted mutation file
        if not os.path.exists(args.weighted_mutation_file):
            raise ValidationError(f"Weighted mutation file not found: {args.weighted_mutation_file}")
    elif args.use_weighted_mutation == False:
        print("Weighted mutation disabled.")
        args.weighted_mutation_file = "None"
    
    # Check if recombination simulation is on
    if not args.use_recombination:
        print("Whether using recombination not provided, using default True")
        args.use_recombination = True
    elif args.use_recombination not in {True, False}:
        print("use_recombination must be either True or False. Using default True")
        args.use_recombination = True
    elif args.use_recombination:
        print("Recombination simulation enabled, validating recombination input...")
        # Check recombination rate 
        if not args.recombination_rate or not isinstance(args.recombination_rate, (int, float)) or args.recombination_rate <= 0:
            raise ValidationError("Recombination rate must be a positive float value. If no recombination simulation is desired, set --use_recombination to False")
        # Check mean recombination size
        if not args.mean_recombination_size or not isinstance(args.mean_recombination_size, float) or args.mean_recombination_size <= 0:
            raise ValidationError("Mean recombination size must be a positive float value")
        # Check minimal recombination size
        if not args.min_recombination_size or not isinstance(args.min_recombination_size, float) or args.min_recombination_size <= 0:
            raise ValidationError("Minimal recombination size must be a positive float value")
        elif args.min_recombination_size >= args.mean_recombination_size:
            raise ValidationError("Minimal recombination size must be smaller than mean recombination size")
        # Check nu
        if not args.nu or not isinstance(args.nu, (int, float)) or args.nu <= 0 or args.nu >=1:
            raise ValidationError("Nu must be a positive float value between 0 and 1")
        else: 
            print(f"Recombination parameters set to: recombination_rate={args.recombination_rate}, mean_recombination_size={args.mean_recombination_size}, minimal_recombination_size={args.min_recombination_size}, nu={args.nu}")
    elif args.use_recombination == False:
        print("Recombination simulation disabled.")
        args.recombination_rate = "None"
        args.mean_recombination_size = "None"
        args.min_recombination_size = "None"
        args.nu = "None"
    
    # Check if gene gain/loss simulation is on
    if not args.use_gain_loss:
        print("Whether using gene gain/loss not provided, using default True")
        args.use_gain_loss = True
    elif args.use_gain_loss not in {True, False}:
        print("use_gain_loss must be either True or False. Using default True")
        args.use_gain_loss = True
    elif args.use_gain_loss == False:
        print("Gene gain/loss simulation disabled.")
        args.gain_rate = "None"
        args.loss_rate = "None"
        args.bin = "None"
        args.mge_data = "None"
        args.position_coverage = "None"
        args.mge_fasta = "None"
        args.mge_entropy = "None"

    if args.use_gain_loss:
        print("Gene gain/loss simulation enabled, validating gene gain/loss input...")
        # Check gain rate
        if not args.gain_rate or not isinstance(args.gain_rate, (int, float)) or args.gain_rate < 0:
            raise ValidationError("Gain rate must be a non-negative float value")
        # Check loss rate
        if not args.loss_rate or not isinstance(args.loss_rate, (int, float)) or args.loss_rate < 0:
            raise ValidationError("Loss rate must be a non-negative float value")
        # Check chromosome bin file
        if not os.path.exists(args.bin):
            raise ValidationError(f"Chromosome bin file not found: {args.bin}")
        # Check MGE data file
        if not os.path.exists(args.mge_data):
            raise ValidationError(f"MGE data file not found: {args.mge_data}")
        # Check MGE fasta file
        if not os.path.exists(args.mge_fasta):
            raise ValidationError(f"MGE fasta file not found: {args.mge_fasta}")
        # Check MGE entropy file
        if not os.path.exists(args.mge_entropy):
            raise ValidationError(f"The directory of MGE entropy file not found: {args.mge_entropy}")
        # Check position coverage file
        if not os.path.exists(args.position_coverage):
            raise ValidationError(f"Position coverage file not found: {args.position_coverage}")
        else:
            print(f"Gene gain/loss parameters set to: gain_rate={args.gain_rate}, loss_rate={args.loss_rate}, bin={args.bin}, mge_data={args.mge_data}, mge_fasta={args.mge_fasta}, mge_entropy={args.mge_entropy}, position_coverage={args.position_coverage}")
        
    elif args.use_gain_loss == False:
        print("Gene gain/loss simulation disabled.")
        args.gain_rate = "None"
        args.loss_rate = "None"
        args.bin = "None"
        args.mge_fasta = "None"
        args.mge_entropy = "None"
        args.position_coverage = "None"