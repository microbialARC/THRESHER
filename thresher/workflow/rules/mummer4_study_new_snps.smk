rule mummer4_new_snps:
    conda:
        os.path.join(BASE_PATH,"envs/mummer4.yaml")
    input:
        original_genome_paths = list(original_genome_path_dict.values()),
        new_genome_paths = list(new_genome_path_dict.values())
    params:
        original_genome_names = list(original_genome_path_dict.keys()),
        new_genome_names = list(new_genome_path_dict.keys()),
        output_dir = config["output"]
    threads:
        config["threads"]
    output:
        expand(os.path.join(config["output"], "mummer4", "{genome_name}_concatenated.report"), genome_name=new_genome_path_dict.keys())
    shell:
        """
        mkdir -p {params.output_dir}/mummer4/
        mkdir -p {params.output_dir}/mummer4/scripts

        # New genomes
        new_genome_names=({params.new_genome_names})
        new_genome_paths=({input.new_genome_paths})

        # Original genomes
        original_genome_names=({params.original_genome_names})
        original_genome_paths=({input.original_genome_paths})

        # The new SNP matrix will include SNP comparisons of:
        # 1. New genomes vs Original genomes
        # 2. New genomes vs New genomes
        # For each comparison, the two genomes will be compared in both directions (A vs B and B vs A)
        # Meaning that each genome will be used as reference once and query once 

        # Create the shell scripts with mummer4(dnadiff) commands for each of the new genomes
        for ((i=0; i<${{#new_genome_names[@]}}; i++)); do
        ref_genome=${{new_genome_names[$i]}}
        ref_path=${{new_genome_paths[$i]}}

        echo "mkdir -p {params.output_dir}/mummer4/${{ref_genome}}" > {params.output_dir}/mummer4/scripts/dnadiff_${{ref_genome}}.sh
        echo "cd {params.output_dir}/mummer4/${{ref_genome}}" >> {params.output_dir}/mummer4/scripts/dnadiff_${{ref_genome}}.sh

        # Note: 
        # If the query is from original genomes, then the query is used as reference once and query once
        # If the query is from the new genomes, no need to repeat the comparison again since it will be done when that genome is used as reference

        # First compare with original genomes

            for ((j=0; j<${{#original_genome_names[@]}}; j++)); do
                # This new genome as reference, original genome as query

                qry_genome=${{original_genome_names[$j]}}
                qry_path=${{original_genome_paths[$j]}}
                echo "dnadiff $ref_path $qry_path -p s_${{ref_genome}}_q_${{qry_genome}}" >> {params.output_dir}/mummer4/scripts/dnadiff_${{ref_genome}}.sh
                echo "sed -n '1p; /AlignedBases/p; /TotalSNPs/p; /TotalGSNPs/p' s_${{ref_genome}}_q_${{qry_genome}}.report > sorted_s_${{ref_genome}}_q_${{qry_genome}}.report" >> {params.output_dir}/mummer4/scripts/dnadiff_${{ref_genome}}.sh

                # Now the original genome as reference, this new genome as query
                echo "dnadiff $qry_path $ref_path -p s_${{qry_genome}}_q_${{ref_genome}}" >> {params.output_dir}/mummer4/scripts/dnadiff_${{ref_genome}}.sh
                echo "sed -n '1p; /AlignedBases/p; /TotalSNPs/p; /TotalGSNPs/p' s_${{qry_genome}}_q_${{ref_genome}}.report > sorted_s_${{qry_genome}}_q_${{ref_genome}}.report" >> {params.output_dir}/mummer4/scripts/dnadiff_${{ref_genome}}.sh
            done
            
        # Then compare with new genomes

            for ((j=0; j<${{#new_genome_names[@]}}; j++)); do
                qry_genome=${{new_genome_names[$j]}}
                qry_path=${{new_genome_paths[$j]}}
                if [ "${{ref_genome}}" != "${{qry_genome}}" ]; then
                    echo "dnadiff $ref_path $qry_path -p s_${{ref_genome}}_q_${{qry_genome}}" >> {params.output_dir}/mummer4/scripts/dnadiff_${{ref_genome}}.sh
                    echo "sed -n '1p; /AlignedBases/p; /TotalSNPs/p; /TotalGSNPs/p' s_${{ref_genome}}_q_${{qry_genome}}.report > sorted_s_${{ref_genome}}_q_${{qry_genome}}.report" >> {params.output_dir}/mummer4/scripts/dnadiff_${{ref_genome}}.sh
                fi
            done

            echo "cat {params.output_dir}/mummer4/${{ref_genome}}/sorted_* > {params.output_dir}/mummer4/${{ref_genome}}_concatenated.report" >> {params.output_dir}/mummer4/scripts/dnadiff_${{ref_genome}}.sh
            echo "rm -rf {params.output_dir}/mummer4/${{ref_genome}}" >> {params.output_dir}/mummer4/scripts/dnadiff_${{ref_genome}}.sh
            # Make the script executable
            chmod +x {params.output_dir}/mummer4/scripts/dnadiff_${{ref_genome}}.sh
        done

        # Create a list of all script files
        ls {params.output_dir}/mummer4/scripts/dnadiff_*.sh > {params.output_dir}/mummer4/scripts/script_list.txt
        # This is critical 
        # No // in the file containing paths to the scripts otherwise there would be error!!!
        sed -i 's#//#/#g' {params.output_dir}/mummer4/scripts/script_list.txt
        # Run the scripts in parallel using GNU parallel
        module load parallel
        parallel --silent --jobs {threads} bash :::: {params.output_dir}/mummer4/scripts/script_list.txt
        rm {params.output_dir}/mummer4/scripts/script_list.txt
        """
    