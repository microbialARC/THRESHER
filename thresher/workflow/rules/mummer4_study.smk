rule mummer4_study:
    conda:
        os.path.join(BASE_PATH,"envs/mummer4.yaml")
    input:
        genome_paths = list(genome_path_dict.values())
    params:
        genome_names = list(genome_path_dict.keys()),
        output_dir = config["output"]
    threads:
        config["threads"]
    output:
        reports = expand(os.path.join(config["output"], "mummer4_study", "{genome_name}_concatenated.report"), genome_name=genome_path_dict.keys())
    shell:
        """
        mkdir -p {params.output_dir}/mummer4_study/
        mkdir -p {params.output_dir}/mummer4_study/scripts
        genome_names=({params.genome_names})
        genome_paths=({input.genome_paths})

        # Create the shell scripts with mummer4(dnadiff) commands for each of the genomes
        for ((i=0; i<${{#genome_names[@]}}; i++)); do
        ref_genome=${{genome_names[$i]}}
        ref_path=${{genome_paths[$i]}}
        echo "mkdir -p {params.output_dir}/mummer4_study/${{ref_genome}}" > {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh
        echo "cd {params.output_dir}/mummer4_study/${{ref_genome}}" >> {params.output_dir}/mummer4_study/scripts/dnadiff_${{ref_genome}}.sh
            for ((j=0; j<${{#genome_names[@]}}; j++)); do
                qry_genome=${{genome_names[$j]}}
                qry_path=${{genome_paths[$j]}}
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

        # Create a list of all script files
        ls {params.output_dir}/mummer4_study/scripts/dnadiff_*.sh > {params.output_dir}/mummer4_study/scripts/script_list.txt
        # This is critical 
        # No // in the file containing paths to the scripts otherwise there would be error!!!
        sed -i 's#//#/#g' {params.output_dir}/mummer4_study/scripts/script_list.txt
        # Run the scripts in parallel using GNU parallel
        module load parallel
        parallel --jobs {threads} bash :::: {params.output_dir}/mummer4_study/scripts/script_list.txt
        rm {params.output_dir}/mummer4_study/scripts/script_list.txt
        """
    