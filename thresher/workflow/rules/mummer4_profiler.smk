rule mummer4_profiler_cmd:
    input:
        topgenomes_filtered = os.path.join(config["output"], "fastani",f"{genome_name}_topgenomes_filtered.txt")
    output:
        script_list = os.path.join(config["output"], "mummer4", "scripts", f"{genome_name}_cmd_list.txt")
    params:
        reference_genome_path = input_genome_path,
        reference_genome_name = genome_name,
        output_dir = config["output"]
    script:
        os.path.join(BASE_PATH,"scripts","mummer4_profiler_cmd.py")

rule mummer4_profiler_exec:
    conda:
        os.path.join(BASE_PATH,"envs/mummer4.yaml")
    input:
        script_list = os.path.join(config["output"], "mummer4", "scripts", f"{genome_name}_cmd_list.txt")
    params:
        output_dir = os.path.join(config["output"], "mummer4","output")
    threads:
        config["threads"]
    output:
        coords_list = os.path.join(config["output"], "mummer4", "output", f"{genome_name}_coords_list.txt"),
        snps_list = os.path.join(config["output"], "mummer4", "output", f"{genome_name}_snps_list.txt")
    shell:
        """
        # Import from snakemake params
        output_dir={params.output_dir}
        coords_list_path={output.coords_list}
        snps_list_path={output.snps_list}
        
        # Parallel execution of MUMmer4 commands from the script list 
        module load parallel    
        parallel --jobs {threads} bash :::: {input.script_list}
        
        # Generate the coords_list and snps_list files as snakemake output checkpoints
        ls $output_dir/*.coords > $coords_list_path
        ls $output_dir/*.snps > $snps_list_path
        """