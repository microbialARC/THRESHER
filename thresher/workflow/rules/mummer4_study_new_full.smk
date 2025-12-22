rule mummer4_new_full_cmd:
    input:
        original_genome_paths = list(original_genome_path_dict.values()),
        new_genome_paths = list(new_genome_path_dict.values())
    params:
        original_genome_names = list(original_genome_path_dict.keys()),
        new_genome_names = list(new_genome_path_dict.keys()),
        output_dir = config["output"]
    output:
        expand(os.path.join(config["output"], "mummer4_study", "{genome_name}_concatenated.report"), genome_name=new_genome_path_dict.keys())
    shell:
        """
        mkdir -p {params.output_dir}/mummer4_study/
        mkdir -p {params.output_dir}/mummer4_study/scripts

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

        echo "mkdir -p {params.output_dir}/mummer4_study/${{ref_genome}}" > {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh
        echo "cd {params.output_dir}/mummer4_study/${{ref_genome}}" >> {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh

        # Note: 
        # If the query is from original genomes, then the query is used as reference once and query once
        # If the query is from the new genomes, no need to repeat the comparison again since it will be done when that genome is used as reference

        # First compare with original genomes

            for ((j=0; j<${{#original_genome_names[@]}}; j++)); do
                # This new genome as reference, original genome as query

                qry_genome=${{original_genome_names[$j]}}
                qry_path=${{original_genome_paths[$j]}}
                echo "dnadiff $ref_path $qry_path -p s_${{ref_genome}}_q_${{qry_genome}}" >> {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh
                echo "sed -n '1p; /AlignedBases/p; /TotalSNPs/p; /TotalGSNPs/p' s_${{ref_genome}}_q_${{qry_genome}}.report > sorted_s_${{ref_genome}}_q_${{qry_genome}}.report" >> {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh

                # Now the original genome as reference, this new genome as query
                echo "dnadiff $qry_path $ref_path -p s_${{qry_genome}}_q_${{ref_genome}}" >> {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh
                echo "sed -n '1p; /AlignedBases/p; /TotalSNPs/p; /TotalGSNPs/p' s_${{qry_genome}}_q_${{ref_genome}}.report > sorted_s_${{qry_genome}}_q_${{ref_genome}}.report" >> {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh
            done
            
        # Then compare with new genomes

            for ((j=0; j<${{#new_genome_names[@]}}; j++)); do
                qry_genome=${{new_genome_names[$j]}}
                qry_path=${{new_genome_paths[$j]}}
                if [ "${{ref_genome}}" != "${{qry_genome}}" ]; then
                    echo "dnadiff $ref_path $qry_path -p s_${{ref_genome}}_q_${{qry_genome}}" >> {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh
                    echo "sed -n '1p; /AlignedBases/p; /TotalSNPs/p; /TotalGSNPs/p' s_${{ref_genome}}_q_${{qry_genome}}.report > sorted_s_${{ref_genome}}_q_${{qry_genome}}.report" >> {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh
                fi
            done

            echo "cat {params.output_dir}/mummer4_study/${{ref_genome}}/sorted_* > {params.output_dir}/mummer4_study/${{ref_genome}}_concatenated.report" >> {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh
            echo "rm -rf {params.output_dir}/mummer4_study/${{ref_genome}}" >> {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh
            # Make the script executable
            chmod +x {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh
        done
        """

rule mummer4_new_full_single:
    conda:
        os.path.join(BASE_PATH,"envs/mummer4.yaml")
    input:
        mummer4_new_full_single_cmd = os.path.join(config["output"], "mummer4_study", "scripts", "dnadiff_{genome_name}.sh")
    output:
        mummer4_new_full_single_result = os.path.join(config["output"], "mummer4_study", "{genome_name}_concatenated.report")
    params:
        output_dir = config["output"]


rule mummer4_new_full_all:
    input:
        mummer4_results = expand(os.path.join(config["output"], "mummer4_study", "{genome_name}_concatenated.report"),genome_name=new_genome_path_dict.keys())
    output:
        mummer4_new_full_done = os.path.join(config["output"], "mummer4_study", ".mummer4_study_complete")
    shell:
        """
        touch {output.mummer4_new_full_done}
        """

    