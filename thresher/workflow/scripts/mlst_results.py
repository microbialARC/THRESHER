import os
import pandas as pd
import glob as glob

#Import from snakemake
species = snakemake.params.species
metadata = pd.read_csv(snakemake.params.metadata, sep='\t', header=None)
analysis_mode = snakemake.params.analysis_mode
from thresher.bin.parse_genome_name import parse_genome_name

# Parse MLST output manually to handle variable column counts
# Because some mlst output lines can have only 3 columns like /path/to/genome - -
# which causes pandas to misalign the columns when reading the file directly with read_csv
mlst_raw_line_path = snakemake.input.mlst_output
mlst_raw_rows = []
with open(mlst_raw_line_path, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) >= 3:
            mlst_raw_rows.append((fields[0], fields[2]))
        elif len(fields) == 2:
            # Edge case: only path and scheme, no ST
            mlst_raw_rows.append((fields[0], '-'))
        else:
            mlst_raw_rows.append((fields[0], '-'))

mlst_raw = pd.DataFrame(mlst_raw_rows, columns=['path', 'ST'])
#Create output dataframe
# Actually we just need genome_name and genome_path from metadata
# But I leave the other columns in for now in case we need them later

if analysis_mode == "lite":
    metadata = metadata[[0,1,2]]
    metadata.columns = ["genome_name", "accession", "genome_path"]
elif analysis_mode == "full":
    metadata.columns = ["genome_name", "accession", "genome_path", "patient_id", "collection_date"]

mlst_results = pd.DataFrame()
mlst_results['genome'] = metadata['genome_name']
mlst_results['genome'] = mlst_results['genome'].apply(lambda x: parse_genome_name(x))

# Create a mapping from genome_path to ST using the raw MLST output
path_to_st = dict(zip(mlst_raw['path'], mlst_raw['ST']))
# Map ST values using genome_path from metadata
mlst_results['ST'] = metadata['genome_path'].map(path_to_st)

#Define functions to get clonal complex or clade
def get_sau_cc(st, db):
    if st == "-":
        return 'Unassigned'
    try:
        st = int(st)
        cc = db.loc[db['ST'] == st, 'clonal_complex'].values[0]
        return 'Unassigned' if pd.isna(cc) else cc
    except (ValueError, IndexError):
        return 'Unassigned'
    
def get_cdiff_clade(st, db):
    if st == "-":
        return 'Unassigned'
    try:
        st = int(st)
        clade = db.loc[db['ST'] == st, 'mlst_clade'].values[0]
        # Return Unassigned if clade is NaN before stringifying (avoids "Cladenan")
        if pd.isna(clade):
            return 'Unassigned'
        # Make the clade character but not number (e.g. Don't show 1.0 but just 1)
        clade = str(clade).split('.')[0]
        # Because for C. difficile, there is no prefix for the clades in the database
        # I add "Clade " here for clarity in the output
        # Return 'Unassigned' if clade is NaN or "", otherwise return "CladeX" where X is the clade number
        return 'Unassigned' if clade == "" or clade.lower() == "nan" else f"Clade{clade}"
    except (ValueError, IndexError):
        return 'Unassigned'
#Add clonal complex or clade to output dataframe
if species == "sau":
    db = pd.read_csv(snakemake.params.mlst_sau_db, sep='\t')
    mlst_results['MLST'] = mlst_results['ST'].apply(lambda x: get_sau_cc(x, db))
elif species == "cdiff":
    db = pd.read_csv(snakemake.params.mlst_cdiff_db, sep='\t')
    mlst_results['MLST'] = mlst_results['ST'].apply(lambda x: get_cdiff_clade(x, db))
elif species == "kp":
    # No clonal complex / clade equivalent for K. pneumoniae for now. Prefix with "ST" for clarity
    mlst_results['MLST'] = mlst_results['ST'].apply(lambda x: 'Unassigned' if x == "-" else f"ST{x}")
elif species == "sepi":
    # No clonal complex / clade equivalent for S. epidermidis for now. Prefix with "ST" for clarity
    mlst_results['MLST'] = mlst_results['ST'].apply(lambda x: 'Unassigned' if x == "-" else f"ST{x}")
#Write output
mlst_results.to_csv(snakemake.output["mlst_results"], sep='\t', index=False)
