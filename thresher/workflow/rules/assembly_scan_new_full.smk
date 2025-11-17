rule assembly_scan_new_full:
    conda:
        os.path.join(BASE_PATH,"envs/assembly_scan.yaml")
    input:
        genome_paths = list(new_genome_path_dict.values())
    output:
        assembly_scan_output = expand(os.path.join(config["output"], "assembly_scan","{genome_name}_assembly_scan.txt"),genome_name=new_genome_path_dict.keys())
    params:
        output_dir=os.path.join(config["output"], "assembly_scan"),
        genome_names=list(new_genome_path_dict.keys())
    shell:
        """
        mkdir -p {params.output_dir}
        genome_paths=({input.genome_paths})
        genome_names=({params.genome_names})

        for i in "${{!genome_names[@]}}"; do
            genome_name_entry="${{genome_names[$i]}}"
            genome_path_entry="${{genome_paths[$i]}}"
            assembly-scan --transpose "$genome_path_entry" > {params.output_dir}/"${{genome_name_entry}}"_assembly_scan.txt
        done
        """