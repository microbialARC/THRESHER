rule snippy_groups:
    conda:
       os.path.join(BASE_PATH,"envs/snippy.yaml")
    input:
        snippy_multi = os.path.join(config["output"], "snippy","scripts","snippy_multi.sh"),
        snippy_reference = os.path.join(config["output"], "snippy","tabs","snippy_reference.txt")
    output:
        snippy_aln = os.path.join(config["output"], "snippy","output","snippy_aln.txt")
    params:
        script_dir = os.path.join(config["output"], "snippy","scripts"),
        output_dir = os.path.join(config["output"], "snippy","output"),
        cleaned_aln_dir = os.path.join(config["output"], "snippy","output","cleaned_aln")
    shell:
        """
        cd {params.script_dir}
        bash {input.snippy_multi}

        # Get list of scripts
        group_scripts=($(ls {params.script_dir}/Group*.sh))
        # Get list of group names
        group_list=($(basename -s ".sh" "${{group_scripts[@]}}"))

        # Read the snippy reference file
        snippy_reference=$(cat {input.snippy_reference})
        # Loop through each script
        mkdir -p {params.cleaned_aln_dir}
        for script in "${{group_scripts[@]}}"; do

            group_name=$(basename "$script" | sed 's#\\.[^.]*$##')
            reference=$(awk -v group="$group_name" 'NR>1 && $1 == group {{print $2}}' <<< "$snippy_reference")
            mkdir -p {params.output_dir}/"${{group_name}}"
            cd {params.output_dir}/"${{group_name}}"
            bash ${{script}}
            
            # Clean the alignment
            snippy-clean_full_aln {params.output_dir}/"${{group_name}}"/core.full.aln > {params.cleaned_aln_dir}/"${{group_name}}"_cleaned.aln
            # Change "Reference" in the alignment to the actual reference genome name
            sed -i "s#Reference#${{reference}}#g" {params.cleaned_aln_dir}/"${{group_name}}"_cleaned.aln

        done

        # Check the snippy output
        # The snippy scripts must match the cleaned alignment to generate the snippy_aln.txt file
        group_aln_path=($(ls {params.cleaned_aln_dir}/*_cleaned.aln))
        group_aln_list=($(basename -s "_cleaned.aln" "${{group_aln_path[@]}}"))
        
        sorted_group_list=($(sort <<< "${{group_list[@]}}"))
        sorted_group_aln_list=($(sort <<< "${{group_aln_list[@]}}"))

        if [ "${{sorted_group_list[*]}}" == "${{sorted_group_aln_list[*]}}" ]; then
            echo "${{group_aln_path[@]}}" | tr ' ' '\n' >> {output.snippy_aln}
        fi

        """