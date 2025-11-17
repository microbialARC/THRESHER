import os
import pandas as pd
import glob as glob

#Import from snakemake
species = snakemake.params.species
mlst_raw = pd.read_csv(snakemake.input.mlst_output, sep=',', header=None)
metadata = pd.read_csv(snakemake.params.metadata, sep='\t', header=None)
metadata.columns = ["genome_name", "accession", "genome_path", "patient_id", "collection_date"]

#Create output dataframe
mlst_results = pd.DataFrame()
mlst_results['genome'] = metadata['genome_name']

# Create a mapping from genome_path to ST
# column 0 is path, column 2 is ST
path_to_st = dict(zip(mlst_raw[0], mlst_raw[2])) 
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
        return 'Unassigned' if pd.isna(clade) else clade
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
    # This is a placeholder because I haven't found something similar to clonal complex or clade for K. pneumoniae
    mlst_results['MLST'] = mlst_results['ST'].apply(lambda x: 'Unassigned' if x == "-" else x)
elif species == "sepi":
    # This is a placeholder because I haven't found something similar to clonal complex or clade for S. epidermidis
    mlst_results['MLST'] = mlst_results['ST'].apply(lambda x: 'Unassigned' if x == "-" else x)
#Write output
mlst_results.to_csv(snakemake.output["mlst_results"], sep='\t', index=False)
