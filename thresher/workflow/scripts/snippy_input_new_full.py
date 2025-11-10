import os
import pandas as pd
import glob

# Import from snakemake
# Compared to the full-pipeline snippy_input.py, this script handles both new genomes and original genomes
# Merge new and original variable and then proceed with the same logic as snippy_input.py

# HC group is updated to include both new and original genomes
hc_group_path = snakemake.input.hc_groups_csv
# Merge new and original study genome paths
new_genome_path = snakemake.params.new_genome_path
original_genome_path = snakemake.params.original_genome_path
genome_path = original_genome_path | new_genome_path
# Merge new and original study accessions
new_study_accession = snakemake.params.new_study_accession
new_study_accession.discard("new")
original_study_accession = snakemake.params.original_study_accession
original_study_accession.discard("new")
study_accession = new_study_accession | original_study_accession
# the global genome path is the path to the directory containing all global genomes used in the analysis
# so no need to merge the two paths here
new_global_genome_path = snakemake.params.new_global_genome_path
original_global_genome_path = snakemake.params.original_global_genome_path
# Merge new and original whatsgnu results
new_whatsgnu_path = snakemake.input.new_whatsgnu_results
original_whatsgnu_path = snakemake.params.original_whatsgnu_results
whatsgnu_path = new_whatsgnu_path + original_whatsgnu_path
# Merge new and original assembly scan results
new_assembly_scan_path = snakemake.input.new_assembly_scan_results
original_assembly_scan_path = snakemake.params.original_assembly_scan_results
assembly_scan_path = new_assembly_scan_path + original_assembly_scan_path
# tab_dir is the directory to store the tab files for new hierachical groups using snippy-multi
tab_dir = snakemake.params.tab_dir
script_dir = snakemake.params.script_dir
threads = snakemake.threads
# Use the maximum memory in the local system for snippy-multi in GB
memory = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') // (1024 ** 3)

os.makedirs(tab_dir, exist_ok=True)
os.makedirs(script_dir, exist_ok=True)
# Create dict of both whatsgnu_path and assembly_scan_path
whatsgnu_path_dict = {os.path.basename(path_entry).replace('_WhatsGNU_topgenomes.txt', ''): path_entry for path_entry in whatsgnu_path}
assembly_scan_path_dict = {os.path.basename(path_entry).replace('_assembly_scan.txt', ''): path_entry for path_entry in assembly_scan_path}

# Snippy groups with no less than 4 genomes
# IQtree: It makes no sense to perform bootstrap with less than 4 sequences
# The sequences include both study genomes and global genomes

hc_group = pd.read_csv(hc_group_path)
snippy_groups = sorted(hc_group['group'].value_counts().index.tolist())

# This iteration will find the genome with largest N50 as reference 
# Create the tab files for each group
# Used for creating the snippy-multi command
# Only groups with total genome count (study + global) >=4 will have tab files created and be included in reference_dict
reference_dict = {}

for group in snippy_groups:
    # Exclude the overlimit genomes from snippy and later IQtree
    group_genomes = hc_group[(hc_group['group'] == group) & (hc_group['overlimit'] != True)]['genome'].tolist()
    group_n50_dict = {}
    group_global_genomes = set()
    for genome in group_genomes:
        # N50 for each group genome
        assembly_scan = pd.read_csv(assembly_scan_path_dict[genome], sep='\t', header=None)
        group_n50_dict[genome] = int(assembly_scan[assembly_scan[1] == 'n50_contig_length'].iloc[0, 2])
        # Introduce top3 global genomes from whatsgnu results for each group genome
        whatsgnu = pd.read_csv(whatsgnu_path_dict[genome], sep='\t', skiprows=1, header=None)
        # Only keep the global genomes with GCA
        whatsgnu = whatsgnu[whatsgnu[0].str.contains('GCA_')]
        top3_genomes = [whatsgnu.iloc[i, 0] for i in range(3)]
        group_global_genomes.update(top3_genomes)
        if study_accession: 
            group_global_genomes = group_global_genomes - study_accession
    # If the count of global genomes and study genomes together is less than 4, skip
    # Normally this should not happen
    # Because by default each study genome brings in 3 global genomes from whatsgnu results
    group_total_genomes_count = len(group_genomes) + len(group_global_genomes)
    if group_total_genomes_count < 4:
        continue
    else:
        reference = max(group_n50_dict, key=group_n50_dict.get)
        reference_dict[group] = reference
        # Create the tab file
        with open(f'{tab_dir}/Group{group}.tab', 'w') as f:
            for genome in group_genomes:
                if genome != reference:
                    f.write(f"{genome}\t{genome_path[genome]}\n")
            for global_genome in group_global_genomes:
                # Search the global genomes in both new and original whatsgnu paths
                # new first, then original
                # If not found in either, raise error

                # I guess there must be a more elegant way to do this, which I will implement in the future if I think of one
                # but for now I am trying to get stuff working asap, and this works
                # so...
                putative_new_path = f"{new_global_genome_path}/{global_genome}.fna"
                putative_original_path = f"{original_global_genome_path}/{global_genome}.fna"
                if os.path.exists(putative_new_path):
                    f.write(f"{global_genome}\t{new_global_genome_path}/{global_genome}.fna\n")
                elif os.path.exists(putative_original_path):
                    f.write(f"{global_genome}\t{original_global_genome_path}/{global_genome}.fna\n")
                else:
                    raise FileNotFoundError(f"Global genome {global_genome} not found in either {new_global_genome_path} or {original_global_genome_path}")
                
# Create the reference.txt file summarizing the reference genome for each group in analysis_groups
with open(f'{tab_dir}/snippy_reference.txt', 'w') as f:
        f.write(f"Group\tReference\n")
        for group, reference in reference_dict.items():
            f.write(f"Group{group}\t{reference}\n")
# Create the snippy-multi command
with open(f'{script_dir}/snippy_multi.sh', 'w') as f:
    for group in reference_dict.keys():
        f.write(f"snippy-multi {tab_dir}/Group{group}.tab --ref {genome_path[reference_dict[group]]} --cpus {threads} --ram {memory} > Group{group}.sh\n")