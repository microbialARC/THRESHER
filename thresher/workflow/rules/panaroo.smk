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
        --threads {threads} \
        > /dev/null 2>&1

        # If core_threshold is set too high
        # "No gene clusters were present above the core frequency threshold! Try adjusting the '--core_threshold' parameter" will be printed
        # In this case, we will fall back to core_threshold = 0.9 and rerun the MSA generation step
        if [ ! -f {params.output_dir}/panaroo/core_gene_alignment_filtered.aln ]; then
            echo "No gene clusters were present above the core frequency threshold set at {params.core_threshold}."
            echo "Falling back to core_threshold = 0.9 and rerunning the MSA generation step"
            panaroo-msa -o {params.output_dir}/panaroo/ \
            --verbose \
            --alignment core \
            --core_threshold 0.9 \
            --aligner mafft \
            --threads {threads} \
            > /dev/null 2>&1
        fi

        # If still no core gene alignment is generated, report and exit
        if [ ! -f {params.output_dir}/panaroo/core_gene_alignment_filtered.aln ]; then
            echo "No core gene alignment was generated after re-running panaroo-msa with core_threshold = 0.9."
            echo "Please check the input genomes."
            exit 1
        fi

        rm -rf {params.output_dir}/panaroo/input
        """