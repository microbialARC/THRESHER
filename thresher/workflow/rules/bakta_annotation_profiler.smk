rule bakta_annotation_profiler:
    conda:
        os.path.join(BASE_PATH,"envs/bakta.yaml")
    input:
        genome_path = input_genome_path,
        bakta_db = [os.path.join(config["output"],"bakta_db","bakta.db")] if config["bakta_db_path"] == "None" else []
    threads:
        config["threads"]
    output:
        gff3 = os.path.join(config["output"], "bakta_annotation",genome_name,f"{genome_name}.gff3"),
        fna = os.path.join(config["output"], "bakta_annotation",genome_name,f"{genome_name}.fna"),
        faa = os.path.join(config["output"], "bakta_annotation",genome_name,f"{genome_name}.faa")
    params:
        genome_name = genome_name,
        species = config["species"],
        output_dir = config["output"],
        db_path = [config["bakta_db_path"] if config["bakta_db_path"] != "None" else os.path.join(config["output"],"bakta_db")]
    shell:
        """
        mkdir -p {params.output_dir}/bakta_annotation/
        mkdir -p {params.output_dir}/bakta_annotation/tmpdir

        genome_name=({params.genome_name})
        genome_path=({input.genome_path})

        # Define species
        if [ {params.species} == "sau" ]; then
            annotation_genus="Staphylococcus"
            annotation_species="aureus"
        elif [ {params.species} == "cdiff" ]; then
            annotation_genus="Clostridioides"
            annotation_species="difficile"
        elif [ {params.species} == "kp" ]; then
            annotation_genus="Klebsiella"
            annotation_species="pneumoniae"
        elif [ {params.species} == "sepi" ]; then
            annotation_genus="Staphylococcus"
            annotation_species="epidermidis"
        fi

        # Run Bakta annotation
        
        bakta --force --db {params.db_path} \
        --output {params.output_dir}/bakta_annotation/{params.genome_name}/ \
        --threads {threads} \
        --prefix {params.genome_name} \
            --tmp-dir {params.output_dir}/bakta_annotation/tmpdir \
            --genus ${{annotation_genus}} \
            --species ${{annotation_species}} \
            --locus-tag {params.genome_name} \
            {input.genome_path}
        
        rm -rf {params.output_dir}/bakta_annotation/tmpdir
        """