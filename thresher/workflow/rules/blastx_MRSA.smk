rule blastx_mrsa_cmd:
    conda:
        os.path.join(BASE_PATH,"envs/blast.yaml")
    input:
        genome_paths = list(genome_path_dict.values())
    params:
        genome_names = list(genome_path_dict.keys()),
        output_dir = os.path.join(config["output"],"blastx","mrsa","output","raw"),
        script_dir = os.path.join(config["output"],"blastx","mrsa","scripts"),
        db = os.path.join(BASE_PATH,"db/MRSA/sequences.fasta")
    output:
        blastx_mrsa_cmd = expand(os.path.join(config["output"],"blastx","mrsa","scripts","{genome_name}_blastx_mrsa.sh"), genome_name=genome_path_dict.keys())
    shell:
        """
        # Create output directory
        mkdir -p {params.output_dir}
        mkdir -p {params.script_dir}
        # Create the database
        makeblastdb -in {params.db} -parse_seqids -blastdb_version 5 -title "MRSA_db" -dbtype prot
        
        # Iterate over the genomes and make blastx commands
        genome_paths=({input.genome_paths})
        genome_names=({params.genome_names})

        for ((i=0; i<${{#genome_names[@]}}; i++)); do
            echo "blastx \
            -outfmt 6 \
            -query ${{genome_paths[$i]}} \
            -db {params.db} \
            -out {params.output_dir}/${{genome_names[$i]}}_blastx_mrsa.tsv" > {params.script_dir}/${{genome_names[$i]}}_blastx_mrsa.sh
            chmod +x {params.script_dir}/${{genome_names[$i]}}_blastx_mrsa.sh
        done
        """

rule blastx_mrsa_single:
    conda:
        os.path.join(BASE_PATH,"envs/blast.yaml")
    input:
        blastx_mrsa_single_cmd = os.path.join(config["output"],"blastx","mrsa","scripts","{genome_name}_blastx_mrsa.sh")
    output:
        blastx_mrsa_single_result = os.path.join(config["output"],"blastx","mrsa","output","raw","{genome_name}_blastx_mrsa.tsv")
    shell:
        """
        bash {input.blastx_mrsa_single_cmd}
        """

rule blastx_mrsa_results:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        blastx_mrsa = expand(os.path.join(config["output"],"blastx","mrsa","output","raw","{genome_name}_blastx_mrsa.tsv"), genome_name=genome_path_dict.keys()),
        plateau_strains_rds = os.path.join(config["output"],"thresher", "output", "plateau_strains.RDS"),
        global_strains_rds = os.path.join(config["output"], "thresher", "output", "global_strains.RDS"),
        peak_strains_rds = os.path.join(config["output"],"thresher", "output",  "peak_strains.RDS"),
        discrepancy_strains_rds = os.path.join(config["output"], "thresher", "output", "discrepancy_strains.RDS")
    output:
        strain_results = os.path.join(config["output"],"blastx","mrsa","output","summary","blastx_MRSA_strains.csv"),
        genome_results = os.path.join(config["output"],"blastx","mrsa","output","summary","blastx_MRSA_genomes.csv")
    params:
        raw_dir = os.path.join(config["output"],"blastx","mrsa","output","raw"),
        summary_dir = os.path.join(config["output"],"blastx","mrsa","output","summary"),
        endpoint = config["endpoint"]
    script:
        os.path.join(BASE_PATH,"scripts","blastx_MRSA.R")