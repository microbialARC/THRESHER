rule iqtree_everything:
    conda:
        os.path.join(BASE_PATH,"envs/iqtree.yaml")
    input:
        core_aln = os.path.join(config["output"],"panaroo","core_gene_alignment_filtered.aln")
    params:
        output_dir = config["output"]
    output:
        tree = os.path.join(config["output"],"iqtree","everything_tree","everything_tree.treefile"),
        contree = os.path.join(config["output"],"iqtree","everything_tree","everything_tree.contree")
    threads:
        config["threads"]
    shell:
        """
        mkdir -p {params.output_dir}/iqtree/everything_tree
        cd {params.output_dir}/iqtree/everything_tree
        iqtree -T {threads} -s {input.core_aln} -m GTR+R --nmax 1000 -ninit 100 --ufboot 1000 --bnni -alrt 1000 --prefix everything_tree
        """
        
