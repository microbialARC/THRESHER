import os
import pandas as pd
import glob

# Import from snakemake
hc_group_path = snakemake.input.hc_groups_csv
genome_path = snakemake.params.genome_path
study_accession = snakemake.params.study_accession
study_accession.discard("new")
global_genome_path = snakemake.params.global_genome_path
whatsgnu_path = snakemake.input.whatsgnu_results
assembly_scan_path = snakemake.input.assembly_scan_results
tab_dir = snakemake.params.tab_dir
script_dir = snakemake.params.script_dir
threads = snakemake.threads
memory = snakemake.params.memory

os.makedirs(tab_dir, exist_ok=True)
os.makedirs(script_dir, exist_ok=True)
# Create dict of both whatsgnu_path and assembly_scan_path
whatsgnu_path_dict = {os.path.basename(path_entry).replace('_WhatsGNU_topgenomes.txt', ''): path_entry for path_entry in whatsgnu_path}
assembly_scan_path_dict = {os.path.basename(path_entry).replace('_assembly_scan.txt', ''): path_entry for path_entry in assembly_scan_path}

# Snippy groups with no less than 2 genomes
hc_group = pd.read_csv(hc_group_path)
snippy_groups = sorted(hc_group['group'].value_counts()[hc_group['group'].value_counts() >= 2].index.tolist())

# This iteration will find the genome with largest N50 as reference 
# Create the tab files for each group
# Used for creating the snippy-multi command
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
        
    reference = max(group_n50_dict, key=group_n50_dict.get)
    reference_dict[group] = reference

    # Create the tab file
    with open(f'{tab_dir}/Group{group}.tab', 'w') as f:
        for genome in group_genomes:
            if genome != reference:
                f.write(f"{genome}\t{genome_path[genome]}\n")
        for global_genome in group_global_genomes:
            f.write(f"{global_genome}\t{global_genome_path}/{global_genome}.fna\n")
# Create the reference.txt file summarizing the reference genome for each group
with open(f'{tab_dir}/snippy_reference.txt', 'w') as f:
        f.write(f"Group\tReference\n")
        for group, reference in reference_dict.items():
            f.write(f"Group{group}\t{reference}\n")
# Create the snippy-multi command
with open(f'{script_dir}/snippy_multi.sh', 'w') as f:
    for group in snippy_groups:
        f.write(f"snippy-multi {tab_dir}/Group{group}.tab --ref {genome_path[reference_dict[group]]} --cpus {threads} --ram {int(memory)} > Group{group}.sh\n")