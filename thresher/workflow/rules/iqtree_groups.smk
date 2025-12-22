rule iqtree_groups:
    conda:
        os.path.join(BASE_PATH,"envs/iqtree.yaml")
    input:
        snippy_aln = os.path.join(config["output"], "snippy","output","snippy_aln.txt")
    output:
        iqtree_group_path = os.path.join(config["output"], "iqtree","group_tree","iqtree_group.txt")
    params:
        iqtree_dir = os.path.join(config["output"], "iqtree","group_tree"),
        cleaned_aln_dir = os.path.join(config["output"], "snippy","output","cleaned_aln"),
        bootstrap_method = config["group_bootstrap_method"],
        bootstrap_number = config["group_bootstrap_number"]
    threads:
        config["threads"]
    shell:
        """
        # Create and go to the output directory
        mkdir -p {params.iqtree_dir}
        cd {params.iqtree_dir}

        # If bootstrap method is "nonparametric", use "--boot " option; if "ultrafast", use "--ufboot" option
        # Redo option is used to overwrite previous runs (can be removed if not desired)
        # Here instead of using parallel in shell, I use a for loop to go through each alignment file and use iqtree's "-T" option to use multiple threads

        for aln in $(cat {input.snippy_aln}); do
            group_name=$(basename "$aln" | sed 's#_cleaned.aln##')
            if [ "{params.bootstrap_method}" == "nonparametric" ]; then
                iqtree -redo -s "${{aln}}" -m GTR+R --boot {params.bootstrap_number} -T AUTO --threads-max {threads} --prefix "${{group_name}}"
            elif [ "{params.bootstrap_method}" == "ultrafast" ]; then
                iqtree -redo -s "${{aln}}" -m GTR+R --ufboot {params.bootstrap_number} -T AUTO --threads-max {threads} --prefix "${{group_name}}"
            fi
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