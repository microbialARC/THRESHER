rule bakta_annotation_new_full_single:
    conda:
        os.path.join(BASE_PATH,"envs/bakta.yaml")
    input:
        new_genome_path = lambda wc: new_genome_path_dict[wc.genome_name],
        bakta_db = [os.path.join(config["output"],"bakta_db","bakta.db")] if config["bakta_db_path"] == "None" else []
    output:
        gff3 = os.path.join(config["output"], "bakta_annotation", "{genome_name}", "{genome_name}.gff3"),
        faa = os.path.join(config["output"], "bakta_annotation", "{genome_name}", "{genome_name}.faa")
    params:
        species = config["species"],
        output_dir = config["output"],
        db_path = [config["bakta_db_path"] if config["bakta_db_path"] != "None" else os.path.join(config["output"],"bakta_db")]
    shell:
        """
        mkdir -p {params.output_dir}/bakta_annotation/
        mkdir -p {params.output_dir}/bakta_annotation/tmpdir
        
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

        bakta --force --db {params.db_path} \
        --output {params.output_dir}/bakta_annotation/{wildcards.genome_name}/ \
        --threads 1 \
        --prefix {wildcards.genome_name} \
        --tmp-dir {params.output_dir}/bakta_annotation/tmpdir \
        --genus ${{annotation_genus}} \
        --species ${{annotation_species}} \
        --locus-tag {wildcards.genome_name} \
        {input.new_genome_path}
        """

rule bakta_annotation_new_full_all:
    input:
        gff3 = [os.path.join(config["output"],"bakta_annotation",f"{genome}",f"{genome}.gff3") for genome in list(new_genome_path_dict.keys())],
        faa = [os.path.join(config["output"],"bakta_annotation",f"{genome}",f"{genome}.faa") for genome in list(new_genome_path_dict.keys())]
    output:
        done = os.path.join(config["output"], "bakta_annotation", ".bakta_complete")
    shell:
        """
        touch {output.done}
        """