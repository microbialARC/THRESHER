rule blastx_mrsa:
    conda:
        os.path.join(BASE_PATH,"envs/blast.yaml")
    input:
        genome_paths = list(genome_path_dict.values())
    threads:
        config["threads"]
    params:
        output_dir = os.path.join(config["output"],"blastx","mrsa","output"),
        script_dir = os.path.join(config["output"],"blastx","mrsa","scripts"),
        db = os.path.join(BASE_PATH,"db/MRSA/sequences.fasta")
    output:
        blastx_mrsa = [os.path.join(config["output"],"blastx","mrsa","output",f"{genome}_blastx_mrsa.tsv") for genome in list(genome_path_dict.keys())]
    shell:
        """
        # Create output directory
        mkdir -p {params.output_dir}
        mkdir -p {params.script_dir}
        # Create the database
        makeblastdb -in {params.db} -parse_seqids -blastdb_version 5 -title "MRSA_db" -dbtype prot
        # Iterate over the genomes and make blastx commands
        genome_paths=({input.genome_paths})

        for genome_path in ${{genome_paths[@]}}; do
            genome_name=$(basename "$genome_path" | sed 's#\\.[^.]*$##')
            echo "blastx -outfmt 6 -query ${{genome_path}} -db {params.db} -out {params.output_dir}/${{genome_name}}_blastx_mrsa.tsv" > {params.script_dir}/${{genome_name}}_blastx_mrsa.sh
        done

        # Parallel
        ls {params.script_dir}/*_blastx_mrsa.sh > {params.script_dir}/script_list.txt
        sed -i 's#//#/#g' {params.script_dir}/script_list.txt
        module load parallel
        parallel --jobs {threads} bash :::: {params.script_dir}/script_list.txt
        """

rule blastx_mrsa_results:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        blastx_mrsa = expand(os.path.join(config["output"],"blastx","mrsa","output","{genome_name}_blastx_mrsa.tsv"), genome_name=genome_path_dict.keys()),
        plateau_strains_rds = os.path.join(config["output"],"thresher", "output", "plateau_strains.RDS"),
        global_strains_rds = os.path.join(config["output"], "thresher", "output", "global_strains.RDS"),
        peak_strains_rds = os.path.join(config["output"],"thresher", "output",  "peak_strains.RDS"),
        discrepancy_strains_rds = os.path.join(config["output"], "thresher", "output", "discrepancy_strains.RDS")
    output:
        strain_results = os.path.join(config["output"],"blastx","mrsa","output","blastx_MRSA_strains.csv"),
        genome_results = os.path.join(config["output"],"blastx","mrsa","output","blastx_MRSA_genomes.csv")
    params:
        output_dir = os.path.join(config["output"],"blastx","mrsa","output"),
        endpoint = config["endpoint"]
    script:
        os.path.join(BASE_PATH,"scripts","blastx_MRSA.R")