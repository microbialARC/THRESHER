rule bakta_annotation_new_full:
    conda:
        os.path.join(BASE_PATH,"envs/bakta.yaml")
    input:
        genome_paths = list(new_genome_path_dict.values()),
        bakta_db = [os.path.join(config["output"],"bakta_db","bakta.db")] if config["bakta_db_path"] == "None" else []
    threads:
        config["threads"]
    output:
        gff3 = [os.path.join(config["output"],"bakta_annotation",f"{genome}",f"{genome}.gff3") for genome in list(new_genome_path_dict.keys())],
        faa = [os.path.join(config["output"],"bakta_annotation",f"{genome}",f"{genome}.faa") for genome in list(new_genome_path_dict.keys())]
    params:
        genome_names = list(new_genome_path_dict.keys()),
        species = config["species"],
        output_dir = config["output"],
        script_dir = os.path.join(config["output"],"bakta_annotation","scripts"),
        db_path = [config["bakta_db_path"] if config["bakta_db_path"] != "None" else os.path.join(config["output"],"bakta_db")]
    shell:
        """
        mkdir -p {params.output_dir}/bakta_annotation/
        mkdir -p {params.output_dir}/bakta_annotation/tmpdir
        mkdir -p {params.script_dir}

        genome_names=({params.genome_names})
        genome_paths=({input.genome_paths})
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

        for ((i=0; i<${{#genome_names[@]}}; i++)); do
            echo "bakta --force --db {params.db_path} \
            --output {params.output_dir}/bakta_annotation/${{genome_names[$i]}}/ \
            --threads {threads} \
            --prefix ${{genome_names[$i]}} \
            --tmp-dir {params.output_dir}/bakta_annotation/tmpdir \
            --genus ${{annotation_genus}} \
            --species ${{annotation_species}} \
            --locus-tag ${{genome_names[$i]}} \
            ${{genome_paths[$i]}}" > {params.script_dir}/${{genome_names[$i]}}_bakta.sh
        done

        ls {params.script_dir}/*_bakta.sh > {params.script_dir}/script_list.txt
        sed -i 's#//#/#g' {params.script_dir}/script_list.txt
        module load parallel
        # Run the scripts in parallel using GNU parallel
        # The number of jobs is set to 5
        parallel --jobs 5 bash :::: {params.script_dir}/script_list.txt
        rm -rf {params.output_dir}/bakta_annotation/tmpdir
        """