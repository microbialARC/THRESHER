"""Input validation functions for THRESHER"""
import os
import pandas as pd
import time

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
    if input_df.shape[0] < 4:
        raise ValidationError("At least 4 genomes required for analysis")
    if input_df["accession"].isnull().any():
        raise ValidationError("Accession column contains missing values. If you don't have GenBank accession, please use 'new'(all lower case)")


    # Print info about number of genomes and patient IDs if applicable
    print(f"Input file read successfully. Number of genomes: {input_df.shape[0]}")
    print(f"Number of unique patient IDs: {input_df['patient_id'].unique().shape[0]}")
    print("Patient ID provided. The cluster plots will be generated based on patient IDs.")
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

    # Check endpoint method
    if args.endpoint not in {"plateau", "peak", "discrepancy", "global"}:
        print("Unsupported endpoint method, using default method 'plateau'")
        args.endpoint = "plateau"
    
    if args.endpoint == "plateau" and not args.plateau_length:
        print("Plateau length not provided, using default length of 15")
        args.plateau_length = 15

    # Check threads, if not provided, use all available threads
    if not args.threads:
        print("Thread number not provided or invalid, using all available threads")
        args.threads = os.cpu_count()
    elif args.threads < 1:
        print("Thread number must be a positive integer, using all available threads")
        args.threads = os.cpu_count()
    elif args.threads > os.cpu_count():
        print(f"Thread number exceeds available threads ({os.cpu_count()}), using all available threads")
        args.threads = os.cpu_count()
# Validate arguments for strain_identifier redo-endpoint mode
def validate_strain_identifier_redo_endpoint(args):
    """Validate the arguments used in redo-endpoint of strain_identifier function"""
    # Check Origninal metadata file (Origninal input file used in full-pipeline mode)
    if not os.path.exists(args.original_metadata):
        raise ValidationError(f"Input file not found: {args.original_metadata}")

    original_metadata_df = pd.read_csv(args.original_metadata, header=None, sep="\t")

    if original_metadata_df.shape[1] != 5:
        raise ValidationError("Input file must have 5 columns")

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
        raise ValidationError("New metadata file must have 3 or 5 columns")
    if new_metadata_df.shape[1] == 3:
        new_metadata_df[4] = None
        new_metadata_df[5] = None
    new_metadata_df.columns = ["genome_name", "accession", "genome_path", "patient_id", "collection_date"]
    if new_metadata_df["accession"].isnull().any():
        raise ValidationError("Accession column contains missing values. If you don't have GenBank accession, please use 'new'(all lower case)")
    # Print info about number of genomes in the new metadata file
    print(f"New metadata file read successfully. Number of genomes: {new_metadata_df.shape[0]}")
    print(f"Number of unique patient IDs in new metadata file: {new_metadata_df['patient_id'].unique().shape[0]}")

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

    
    # Check species
    if not args.species:
        raise ValidationError("Species must be provided")
    elif args.species not in {"sau", "sepi", "cdiff", "kp"}:
        raise ValidationError(f"Species must be one of: sau, sepi, cdiff, kp")
    
    # Check threads, if not provided, use all available threads
    if not args.threads:
        print("Thread number not provided or invalid, using all available threads")
        args.threads = os.cpu_count()
    elif args.threads < 1:
        print("Thread number must be a positive integer, using all available threads")
        args.threads = os.cpu_count()
    elif args.threads > os.cpu_count():
        print(f"Thread number exceeds available threads ({os.cpu_count()}), using all available threads")
        args.threads = os.cpu_count()
        
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
        raise ValidationError("New metadata file must have 3 or 5 columns")
    if new_metadata_df.shape[1] == 3:
        new_metadata_df[4] = None
        new_metadata_df[5] = None
    new_metadata_df.columns = ["genome_name", "accession", "genome_path", "patient_id", "collection_date"]
    if new_metadata_df["accession"].isnull().any():
        raise ValidationError("Accession column contains missing values. If you don't have GenBank accession, please use 'new'(all lower case)")
    # Print info about number of genomes in the new metadata file
    print(f"New metadata file read successfully. Number of genomes: {new_metadata_df.shape[0]}")
    print(f"Number of unique patient IDs in new metadata file: {new_metadata_df['patient_id'].unique().shape[0]}")

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
    
    # Check threads, if not provided, use all available threads
    if not args.threads:
        print("Thread number not provided or invalid, using all available threads")
        args.threads = os.cpu_count()
    elif args.threads < 1:
        print("Thread number must be a positive integer, using all available threads")
        args.threads = os.cpu_count()
    elif args.threads > os.cpu_count():
        print(f"Thread number exceeds available threads ({os.cpu_count()}), using all available threads")
        args.threads = os.cpu_count()

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
    if not (95 <= args.ani_threshold <= 100):
        raise ValidationError("ani_threshold must be between 95 and 100 for reliable profiling")
    
    # Check conda prefix
    if not args.conda_prefix:
        print("Conda environment path not provided, using default path: <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>")
        args.conda_prefix = os.path.abspath(f"{args.output}/conda_envs_{args.prefix}")
    
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
            print(f"Bakta database will be downloaded to {args.output}/bakta/db")
            args.bakta_db_path = "None"
        else:
            print(f"Bakta database file found at: {bakta_db_file}")

    if args.bakta_db_type not in {"full", "light"}:
        print("Unsupported Bakta database type, using default 'full'")
        args.bakta_db_type = "full"
# Validate the shared arguments for evo_simulator function by preset mode
def validate_evo_simulator_preset(args):
    """Validate the shared arguments used in evo_simulator function by custom mode"""
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
    elif args.seed < 0:
        raise ValidationError("Seed must be a non-negative integer")
    
    # Check output directory
    if not args.output:
        print(f"Output directory not provided, creating a directory named 'thresher_evo_simulator_output_{args.prefix} under current working directory({os.getcwd()}) as output directory")
        args.output = os.path.join(os.getcwd(), f"thresher_evo_simulator_output_{args.prefix}")
    if not os.path.exists(args.output):
        print(f"Output directory {args.output} does not exist, creating it.")
        os.makedirs(args.output)
        args.output = os.path.abspath(args.output)
        
    # Check years
    if not args.years:
        raise ValidationError("Years must be provided")
    elif args.years <= 0:
        raise ValidationError("Years must be a positive integer")

    if not args.taxa:
        raise ValidationError("Taxa (Offspring) must be provided")
    elif args.taxa < 1:
        raise ValidationError("Taxa (Offspring) must be an integer no less than 1")
    
    # Check weighted_mutation
    if not args.weighted_mutation:
        print("Whether using weighted_mutation not provided, using default True")
        args.weighted_mutation = True
    elif args.weighted_mutation not in {True, False}:
        print("weighted_mutation must be either True or False. Using default True")
        args.weighted_mutation = True

    # Check recombination
    if not args.recombination:
        print("Whether using recombination not provided, using default True")
        args.recombination = True
    elif args.recombination not in {True, False}:
        print("recombination must be either True or False. Using default True")
        args.recombination = True

    # Check gain_loss
    if not args.gain_loss:
        print("Whether using gain_loss not provided, using default True")
        args.gain_loss = True
    elif args.gain_loss not in {True, False}:
        print("gain_loss must be either True or False. Using default True")
        args.gain_loss = True
    
    # Check species
    if not args.species:
        raise ValidationError("Species must be provided")
    elif args.species not in {"sau", "sepi", "cdiff", "kp"}:
        raise ValidationError(f"Species must be one of: sau, sepi, cdiff, kp")
    
    # Check ST
    if not args.st:
        raise ValidationError("ST must be provided")
    
    # Check threads, if not provided, use all available threads
    if not args.threads:
        print("Thread number not provided or invalid, using all available threads")
        args.threads = os.cpu_count()
    elif args.threads < 1:
        print("Thread number must be a positive integer, using all available threads")
        args.threads = os.cpu_count()
    elif args.threads > os.cpu_count():
        print(f"Thread number exceeds available threads ({os.cpu_count()}), using all available threads")
        args.threads = os.cpu_count()
    # Check the genome_profiler database path (TBD)
# Validate the arguments for evo_simulator function by custom mode
def validate_evo_simulator_custom(args):
    """Validate the prefix argument used in evo_simulator function by preset mode"""
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
        print(f"Output directory not provided, creating a directory named 'thresher_evo_simulator_output_{args.prefix} under current working directory({os.getcwd()}) as output directory")
        args.output = os.path.join(os.getcwd(), f"thresher_evo_simulator_output_{args.prefix}")
    if not os.path.exists(args.output):
        print(f"Output directory {args.output} does not exist, creating it.")
        os.makedirs(args.output)
        args.output = os.path.abspath(args.output)
    
    # Check ancestor_genome
    if not os.path.exists(args.ancestor_genome):
        raise ValidationError(f"Ancestor genome file not found: {args.ancestor_genome}")
    
    # Check mutation rate
    if not args.mutation_rate or args.mutation_rate <= 0:
        raise ValidationError("Mutation rate must be a positive float value")
    
    # Check recombination rate
    if not args.recombination_rate or args.recombination_rate < 0:
        raise ValidationError("Recombination rate must be a non-negative float value")
    
    # Check mean recombination size
    if not args.mean_recombination_size or args.mean_recombination_size <= 0:
        raise ValidationError("Mean recombination size must be a positive integer")
    
    # Check gain rate
    if not args.gain_rate or args.gain_rate < 0:
        raise ValidationError("Gain rate must be a non-negative float value")
    # Check loss rate
    if not args.loss_rate or args.loss_rate < 0:
        raise ValidationError("Loss rate must be a non-negative float value")
    # Check weighted_mutation_file
    if not os.path.exists(args.weighted_mutation_file):
        raise ValidationError(f"Weighted mutation file not found: {args.weighted_mutation_file}")
    # Check choromosome bin file
    if not os.path.exists(args.bin):
        raise ValidationError(f"Chromosome bin file not found: {args.bin}")
    # Check MGE fasta file
    if not os.path.exists(args.mge_fasta):
        raise ValidationError(f"MGE fasta file not found: {args.mge_fasta}")
    # Check MGE entropy file
    if not os.path.exists(args.mge_entropy):
        raise ValidationError(f"The directory of MGE entropy file not found: {args.mge_entropy}")
    # Check threads, if not provided, use all available threads
    if not args.threads:
        print("Thread number not provided or invalid, using all available threads")
        args.threads = os.cpu_count()
    elif args.threads < 1:
        print("Thread number must be a positive integer, using all available threads")
        args.threads = os.cpu_count()
    elif args.threads > os.cpu_count():
        print(f"Thread number exceeds available threads ({os.cpu_count()}), using all available threads")
        args.threads = os.cpu_count()
    