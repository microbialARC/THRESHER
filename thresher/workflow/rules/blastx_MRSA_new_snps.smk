rule blastx_mrsa_new:
    conda:
        os.path.join(BASE_PATH,"envs/blast.yaml")
    input:
        genome_paths = list(new_genome_path_dict.values())
    threads:
        config["threads"]
    params:
        genome_names = list(new_genome_path_dict.keys()),
        output_dir = os.path.join(config["output"],"blastx","mrsa","output","raw"),
        script_dir = os.path.join(config["output"],"blastx","mrsa","scripts"),
        db = os.path.join(BASE_PATH,"db/MRSA/sequences.fasta")
    output:
        blastx_mrsa = [os.path.join(config["output"],"blastx","mrsa","output","raw",f"{genome}_blastx_mrsa.tsv") for genome in list(new_genome_path_dict.keys())]
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
        done

        # Parallel
        ls {params.script_dir}/*_blastx_mrsa.sh > {params.script_dir}/script_list.txt
        sed -i 's#//#/#g' {params.script_dir}/script_list.txt
        module load parallel
        parallel --silent --jobs {threads} bash :::: {params.script_dir}/script_list.txt
        """

rule blastx_mrsa_results_new_snps:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        blastx_mrsa = expand(os.path.join(config["output"],"blastx","mrsa","output","raw","{genome_name}_blastx_mrsa.tsv"), genome_name=new_genome_path_dict.keys()),
        plateau_strains_rds = os.path.join(config["output"],"thresher", "output", "new_plateau.RDS"),
        global_strains_rds = os.path.join(config["output"], "thresher", "output", "new_global.RDS"),
        peak_strains_rds = os.path.join(config["output"],"thresher", "output",  "new_peak.RDS"),
        discrepancy_strains_rds = os.path.join(config["output"], "thresher", "output", "new_discrepancy.RDS")
    output:
        # The updated genome results are universal input for all endpoints
        updated_genome_results = os.path.join(config["output"],"blastx","mrsa","output","summary","blastx_MRSA_genomes.csv"),
        plateau_strain_results = os.path.join(config["output"],"blastx","mrsa","output","summary","blastx_MRSA_plateau_strains.csv"),
        global_strain_results = os.path.join(config["output"],"blastx","mrsa","output","summary","blastx_MRSA_global_strains.csv"),
        peak_strain_results = os.path.join(config["output"],"blastx","mrsa","output","summary","blastx_MRSA_peak_strains.csv"),
        discrepancy_strain_results = os.path.join(config["output"],"blastx","mrsa","output","summary","blastx_MRSA_discrepancy_strains.csv")
    params:
        raw_dir = os.path.join(config["output"],"blastx","mrsa","output","raw"),
        summary_dir = os.path.join(config["output"],"blastx","mrsa","output","summary"),
        thresher_output = config["thresher_output"]
    script:
        os.path.join(BASE_PATH,"scripts","blastx_MRSA_new_snps.R")