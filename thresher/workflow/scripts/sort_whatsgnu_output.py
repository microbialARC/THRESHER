import glob
import pandas as pd
import os

#Import from snakemake
whatsgnu_output_dir = snakemake.params.whatsgnu_output_dir
script_path = snakemake.params.script_path
output_path = snakemake.params.output_path
study_genome_dict = snakemake.params.study_genome_dict
global_genome_dir = snakemake.params.global_genome_dir
# The Genbank accession ID of the study genomes
study_accession = snakemake.params.study_accession
# The actual downloaded top genomes list
actual_download_topgenomes_path = snakemake.input.actual_download_topgenomes
with open(actual_download_topgenomes_path, 'r') as f:
    actual_download_topgenomes = set(line.strip() for line in f if line.strip())

# get all files in the whatsgnu output directory
whatsgnu_output_file = glob.glob(os.path.join(whatsgnu_output_dir,"*","*_WhatsGNU_topgenomes.txt"))
os.makedirs(script_path, exist_ok=True)

for file in whatsgnu_output_file:
    study_genome = os.path.basename(file).replace("_WhatsGNU_topgenomes.txt", "")
    global_genomes_df = pd.read_csv(file, sep='\t', header=None, skiprows=1).iloc[:, 0]
    global_genomes = [genome for genome in global_genomes_df if 'GCA_' in genome]

    # Because for some whatsgnu database, the genome name may contain different naming styles
    # First split by "_", then find the position of "GCA".
    # Only take the "GCA" and the next part (9-digit accession)
    # Normalize genome names: extract the GCA accession and its following part if present
    extracted = set()
    for genome_entry in global_genomes:
        genome_name_parts = genome_entry.split("_")
        # find first part starting with GCA
        gca_index = [i for i, part in enumerate(genome_name_parts) if part.startswith("GCA")]
        gca_index_value = gca_index[0]
        extracted.add("_".join(genome_name_parts[gca_index_value:gca_index_value+2]))
    
    global_genomes = list(extracted)
    
    # Only keep the global genomes that are not the same as the study genome
    # Also ensure they are in the actual downloaded top genomes
    global_genomes = [genome for genome in global_genomes if genome in actual_download_topgenomes]
    global_genomes = [genome for genome in global_genomes if genome not in study_accession]
    
    with open(os.path.join(script_path, f'dnadiff_{study_genome}.sh'), 'w') as f:
        if len(global_genomes) > 0:
            f.write(f"mkdir -p {output_path}/{study_genome}\n")
            f.write(f"cd {output_path}/{study_genome}\n")
            for global_genome_entry in global_genomes:
                f.write(f"dnadiff {study_genome_dict[study_genome]} {global_genome_dir}/{global_genome_entry}.fna -p s_{study_genome}_q_{global_genome_entry} \n")
                f.write(f"sed -n '1p; /AlignedBases/p; /TotalSNPs/p; /TotalGSNPs/p' s_{study_genome}_q_{global_genome_entry}.report > sorted_s_{study_genome}_q_{global_genome_entry}.report\n")
                f.write(f"dnadiff {global_genome_dir}/{global_genome_entry}.fna {study_genome_dict[study_genome]} -p s_{global_genome_entry}_q_{study_genome} \n")
                f.write(f"sed -n '1p; /AlignedBases/p; /TotalSNPs/p; /TotalGSNPs/p' s_{global_genome_entry}_q_{study_genome}.report > sorted_s_{global_genome_entry}_q_{study_genome}.report\n")
            f.write(f"cat {output_path}/{study_genome}/sorted*.report > {output_path}/{study_genome}_concatenated.report\n")
            f.write(f"rm -rf {output_path}/{study_genome}/\n")
        else:
            f.write(f"echo 'No global genomes found for {study_genome} other than itself.'\n")
            # Generate an empty concatenated report
            f.write(f"echo '' > {output_path}/{study_genome}_concatenated.report\n")
