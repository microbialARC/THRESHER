rule iqtree_comprehensive:
    conda:
        os.path.join(BASE_PATH,"envs/iqtree.yaml")
    input:
        core_aln = os.path.join(config["output"],"panaroo","core_gene_alignment_filtered.aln")
    params:
        output_dir = config["output"],
        bootstrap_method = config["core_bootstrap_method"],
        bootstrap_number = config["core_bootstrap_number"]
    output:
        tree = os.path.join(config["output"],"iqtree","comprehensive_tree","comprehensive_tree.treefile"),
        contree = os.path.join(config["output"],"iqtree","comprehensive_tree","comprehensive_tree.contree")
    threads:
        config["threads"]
    shell:
        """
        mkdir -p {params.output_dir}/iqtree/comprehensive_tree
        cd {params.output_dir}/iqtree/comprehensive_tree
        # if bootstrap method is "nonparametric", use "--boot " option; if "ultrafast", use "--ufboot" option
        # put -redo to overwrite previous runs
        
        if [ "{params.bootstrap_method}" == "nonparametric" ]; then
            iqtree -redo -T {threads} -s {input.core_aln} -m GTR+R --boot {params.bootstrap_number} --bnni -alrt {params.bootstrap_number} --prefix comprehensive_tree
        elif [ "{params.bootstrap_method}" == "ultrafast" ]; then
            iqtree -redo -T {threads} -s {input.core_aln} -m GTR+R --ufboot {params.bootstrap_number} --bnni -alrt {params.bootstrap_number} --prefix comprehensive_tree
        fi
        """
 