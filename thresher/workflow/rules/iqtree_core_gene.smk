rule iqtree_core_gene:
    conda:
        os.path.join(BASE_PATH,"envs/iqtree.yaml")
    input:
        core_aln = os.path.join(config["output"],"panaroo","core_gene_alignment_filtered.aln")
    params:
        output_dir = config["output"],
        bootstrap_method = config["core_bootstrap_method"],
        bootstrap_number = config["core_bootstrap_number"]
    output:
        tree = os.path.join(config["output"],"iqtree","core_gene_tree","core_gene_tree.treefile"),
        contree = os.path.join(config["output"],"iqtree","core_gene_tree","core_gene_tree.contree")
    threads:
        config["threads"]
    shell:
        """
        mkdir -p {params.output_dir}/iqtree/core_gene_tree
        cd {params.output_dir}/iqtree/core_gene_tree
        # if bootstrap method is "nonparametric", use "--boot " option; if "ultrafast", use "--ufboot" option
        # put -redo to overwrite previous runs
        
        if [ "{params.bootstrap_method}" == "nonparametric" ]; then
            iqtree -redo -T AUTO --threads-max {threads} -s {input.core_aln} -m GTR+R --boot {params.bootstrap_number} --bnni -alrt {params.bootstrap_number} --prefix core_gene_tree
        elif [ "{params.bootstrap_method}" == "ultrafast" ]; then
            iqtree -redo -T AUTO --threads-max {threads} -s {input.core_aln} -m GTR+R --ufboot {params.bootstrap_number} --bnni -alrt {params.bootstrap_number} --prefix core_gene_tree
        fi
        """
 