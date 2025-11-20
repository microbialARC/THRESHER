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
# Use the maximum memory in the system for snippy-multi in GB
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
# Only groups with total genome count (study + global) >=4 or group study genome more than 1 
# will have tab files created and be included in reference_dict
reference_dict = {}
for group in snippy_groups:
    # Exclude the overlimit genomes from snippy and later IQtree
    group_genomes = hc_group[(hc_group['group'] == group) & (hc_group['overlimit'] != True)]['genome'].tolist()
    # Only consider the groups with more than 1 study genome
    if len(group_genomes) <= 1:
        continue
    else:
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
            whatsgnu_genomes = set(whatsgnu[0].unique())
            # Because for some whatsgnu database, the genome name may contain different naming styles
            # First split by "_", then find the position of "GCA".
            # Only take the "GCA" and the next part (9-digit accession)
            # Normalize genome names: extract the GCA accession and its following part if present
            extracted = set()
            for genome_entry in whatsgnu_genomes:
                genome_name_parts = genome_entry.split("_")
                # find first part starting with GCA
                gca_index = [i for i, part in enumerate(genome_name_parts) if part.startswith("GCA")]
                gca_index_value = gca_index[0]
                extracted.add("_".join(genome_name_parts[gca_index_value:gca_index_value+2]))

            top3_genomes = [list(extracted)[i] for i in range(3)]
            group_global_genomes.update(top3_genomes)
            if study_accession: 
                group_global_genomes = group_global_genomes - study_accession
        # If the count of global genomes and study genomes together is less than 4, skip
        # This normally will not happen because each study genome brings in 3 global genomes
        group_total_genomes_count = len(group_genomes) + len(group_global_genomes)

        if group_total_genomes_count < 4:
            continue
        else:
            # Find the reference genome with largest N50
            reference = max(group_n50_dict, key=group_n50_dict.get)
            reference_dict[group] = reference
            # Create the tab file
            with open(f'{tab_dir}/Group{group}.tab', 'w') as f:
                for genome in group_genomes:
                    if genome != reference:
                        f.write(f"{genome}\t{genome_path[genome]}\n")
                for global_genome in group_global_genomes:
                    f.write(f"{global_genome}\t{global_genome_path}/{global_genome}.fna\n")

# Create the reference.txt file summarizing the reference genome for each group in analysis_groups

with open(f'{tab_dir}/snippy_reference.txt', 'w') as f:
        f.write(f"Group\tReference\n")
        for group, reference in reference_dict.items():
            f.write(f"Group{group}\t{reference}\n")
# Create the snippy-multi command
with open(f'{script_dir}/snippy_multi.sh', 'w') as f:
    for group in reference_dict.keys():
        f.write(f"snippy-multi {tab_dir}/Group{group}.tab --ref {genome_path[reference_dict[group]]} --cpus {threads} --ram {memory} > Group{group}.sh\n")