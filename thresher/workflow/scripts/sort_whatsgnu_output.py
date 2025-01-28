import glob
import pandas as pd
import os

#Import from snakemake
whatsgnu_output_dir = snakemake.params.whatsgnu_output_dir
script_path = snakemake.params.script_path
output_path = snakemake.params.output_path
study_genome_dict = snakemake.params.stduy_genome_dict
global_genome_dir = snakemake.params.global_genome_dir

# get all files in the whatsgnu output directory
whatsgnu_output_file = glob.glob(os.path.join(whatsgnu_output_dir,"*","*_WhatsGNU_topgenomes.txt"))
os.makedirs(script_path, exist_ok=True)

for file in whatsgnu_output_file:
    study_genome = os.path.basename(file).replace("_WhatsGNU_topgenomes.txt", "")
    global_genomes = pd.read_csv(file, sep='\t', header=None, skiprows=1).iloc[:, 0]
    # Only keep the top 10 global genomes with GCA
    global_genomes = [genome for genome in global_genomes if 'GCA_' in genome][:10]
    # Only keep the genome 
    with open(os.path.join(script_path, f'dnadiff_{study_genome}.sh'), 'w') as f:
        f.write(f"mkdir -p {output_path}/{study_genome}\n")
        f.write(f"cd {output_path}/{study_genome}\n")
        for global_genome_entry in global_genomes:
            f.write(f"dnadiff {study_genome_dict[study_genome]} {global_genome_dir}/{global_genome_entry}.fna -p s_{study_genome}_q_{global_genome_entry} \n")
            f.write(f"sed -n '1p; /AlignedBases/p; /TotalSNPs/p; /TotalGSNPs/p' s_{study_genome}_q_{global_genome_entry}.report > sorted_s_{study_genome}_q_{global_genome_entry}.report\n")
            f.write(f"dnadiff {global_genome_dir}/{global_genome_entry}.fna {study_genome_dict[study_genome]} -p s_{global_genome_entry}_q_{study_genome} \n")
            f.write(f"sed -n '1p; /AlignedBases/p; /TotalSNPs/p; /TotalGSNPs/p' s_{global_genome_entry}_q_{study_genome}.report > sorted_s_{global_genome_entry}_q_{study_genome}.report\n")
        f.write(f"cat {output_path}/{study_genome}/sorted*.report > {output_path}/{study_genome}_concatenated.report\n")
        f.write(f"rm -rf {output_path}/{study_genome}/\n")