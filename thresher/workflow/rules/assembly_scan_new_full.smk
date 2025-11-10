rule assembly_scan_new_full:
    conda:
        os.path.join(BASE_PATH,"envs/assembly_scan.yaml")
    input:
        genome_paths = list(new_genome_path_dict.values())
    output:
        assembly_scan_output = expand(os.path.join(config["output"], "assembly_scan","{genome_name}_assembly_scan.txt"),genome_name=new_genome_path_dict.keys())
    params:
        output_dir=os.path.join(config["output"], "assembly_scan")
    shell:
        """
        mkdir -p {params.output_dir}
        for genome in {input.genome_paths}; do
            assembly-scan --transpose $genome > {params.output_dir}/$(basename $genome | sed 's#\\.[^.]*$##')_assembly_scan.txt
        done
        """