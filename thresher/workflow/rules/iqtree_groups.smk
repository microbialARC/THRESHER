rule iqtree_groups:
    conda:
        os.path.join(BASE_PATH,"envs/iqtree.yaml")
    input:
        snippy_aln = os.path.join(config["output"], "snippy","output","snippy_aln.txt")
    output:
        iqtree_group_path = os.path.join(config["output"], "iqtree","group_tree","iqtree_group.txt")
    params:
        iqtree_dir = os.path.join(config["output"], "iqtree","group_tree"),
        cleaned_aln_dir = os.path.join(config["output"], "snippy","output","cleaned_aln")
    threads:
        config["threads"]
    shell:
        """
        # Create and go to the output directory
        mkdir -p {params.iqtree_dir}
        cd {params.iqtree_dir}

        # Run iqtree for each group
        for aln in $(cat {input.snippy_aln}); do
            group_name=$(basename "$aln" | sed 's#_cleaned.aln##')
            iqtree -redo -s "${{aln}}" -m GTR+R -ninit 100 --boot 100 -T {threads} --prefix "${{group_name}}"
        done 

        # Check if every group has an iqtree output

        iqtree_path=($(ls {params.iqtree_dir}/*.contree))
        iqtree_group=($(basename -s ".contree" "${{iqtree_path[@]}}"))

        group_aln_path=($(ls {params.cleaned_aln_dir}/*_cleaned.aln))
        group_aln_list=($(basename -s "_cleaned.aln" "${{group_aln_path[@]}}"))
        sorted_group_aln_list=($(sort <<< "${{group_aln_list[@]}}"))

        if [ "${{iqtree_group[*]}}" == "${{sorted_group_aln_list[*]}}" ]; then
            echo "${{iqtree_path[@]}}" | tr ' ' '\n' >> {output.iqtree_group_path}
        fi
        """