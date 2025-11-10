rule panaroo_new_full:
    conda:
        os.path.join(BASE_PATH,"envs/panaroo.yaml")
    input:
        new_gff3 = [os.path.join(config["output"],"bakta_annotation",f"{genome}",f"{genome}.gff3") for genome in list(new_genome_path_dict.keys())]
    params:
        thresher_output = config["thresher_output"],
        output_dir = config["output"]
    output:
        panaroo_new_graph = os.path.join(config["output"],"panaroo_new","final_graph.gml"),
        core_aln = os.path.join(config["output"],"panaroo","core_gene_alignment_filtered.aln")
    threads:
        config["threads"]
    shell:
        """
        # For this rule, we run the pangenome analysis only on the new genomes 
        # And then we use the 'Merge Panaroo graphs' function to combine the new analysis with the previous one
        # Details of merge function can be found at https://gthlab.au/panaroo/#/merge/merge_graphs
        # The MSA will be generated using the merged graph and panaroo-msa function

        # First, run Panaroo on new genomes only
        # For this step we don't need to get the MSA yet
        # Because we will use the 'Merge Panaroo graphs' function and get the combined MSA later

        mkdir -p {params.output_dir}/panaroo_new/
        mkdir -p {params.output_dir}/panaroo_new/input/
        cp {input.new_gff3} {params.output_dir}/panaroo_new/input/
        panaroo -i {params.output_dir}/panaroo_new/input/*.gff3 \
        -o {params.output_dir}/panaroo_new/ \
        --remove-invalid-genes \
        --clean-mode strict \
        --core_threshold 0.95 \
        --aligner none \
        --family_threshold 0.7 \
        --refind_prop_match 0.5 \
        --search_radius 5000 \
        --threads {threads}
        rm -rf {params.output_dir}/panaroo_new/input/

        # Merge Panaroo graphs

        panaroo-merge -d {params.thresher_output}/panaroo/ {params.output_dir}/panaroo_new/ \
        -o {params.output_dir}/panaroo/ \
        --threads {threads}

        # Generate MSA using the merged graph

        panaroo-msa -o {params.output_dir}/panaroo/ \
        --verbose \
        --alignment core \
        --core_threshold 0.95 \
        --aligner mafft \
        --threads {threads} \
        
        # Keep the panaroo_new directory for record
        """