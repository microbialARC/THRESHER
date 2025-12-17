rule panaroo:
    conda:
        os.path.join(BASE_PATH,"envs/panaroo.yaml")
    input:
        bakta_complete = os.path.join(config["output"], "bakta_annotation", ".bakta_complete")
    params:
        output_dir = config["output"],
        gff3 = [os.path.join(config["output"],"bakta_annotation",f"{genome}",f"{genome}.gff3") for genome in list(genome_path_dict.keys())],
        core_threshold = config["core_threshold"]
    output:
        core_aln = os.path.join(config["output"],"panaroo","core_gene_alignment_filtered.aln")
    threads:
        config["threads"]
    shell:
        """
        mkdir -p {params.output_dir}/panaroo/
        mkdir -p {params.output_dir}/panaroo/input
        cp {params.gff3} {params.output_dir}/panaroo/input/
        panaroo -i {params.output_dir}/panaroo/input/*.gff3 \
        -o {params.output_dir}/panaroo/ \
        --remove-invalid-genes \
        --clean-mode strict \
        --alignment core \
        --core_threshold {params.core_threshold} \
        --aligner mafft \
        --family_threshold 0.7 \
        --refind_prop_match 0.5 \
        --search_radius 5000 \
        --threads {threads}
        rm -rf {params.output_dir}/panaroo/input
        """