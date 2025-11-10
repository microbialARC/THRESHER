import glob
import pandas as pd
import os

input_dir = snakemake.params.input_dir
output_dir = snakemake.params.output_dir
study_accession = snakemake.params.study_accession
whatsgnu_output = glob.glob(os.path.join(input_dir, "**/*_WhatsGNU_topgenomes.txt"), recursive=True)

unique_topgenomes = set()

# Read and process WhatsGNU output files
for file in whatsgnu_output:
    topgenomes = pd.read_csv(file, sep='\t', header=None, skiprows=1, names=['genome','score'])
    topgenomes = topgenomes.genome.unique()
    unique_topgenomes.update(topgenomes)
# Exclude study genomes and non-Genbank genomes
unique_topgenomes = [genome for genome in unique_topgenomes if genome.startswith("GCA")]

study_accession.discard("new")
if study_accession:
    unique_topgenomes = unique_topgenomes - study_accession


# Write the filtered results
with open(os.path.join(output_dir, 'topgenomes.txt'), 'w') as f:
    f.write('\n'.join(unique_topgenomes))

os.makedirs(os.path.join(output_dir, "scripts"), exist_ok=True)
for genome in unique_topgenomes:
    with open(os.path.join(output_dir, 'scripts', f'datasets_{genome}.sh'), 'w') as f:
        f.write(f"datasets download genome accession {genome} --filename {os.path.join(output_dir, genome + '.zip')}\n")
        f.write(f"unzip {os.path.join(output_dir, genome + '.zip')} -d {os.path.join(output_dir, genome)}\n")
        f.write(f"mv {os.path.join(output_dir, genome, 'ncbi_dataset', 'data', genome, genome + '*.fna')} {os.path.join(output_dir, genome + '.fna')}\n")
        f.write(f"rm {os.path.join(output_dir, genome + '.zip')}\n")
        f.write(f"rm -rf {os.path.join(output_dir, genome)}\n")