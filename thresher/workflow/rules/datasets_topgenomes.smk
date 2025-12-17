rule dataset_topgenomes:
    conda:
        os.path.join(BASE_PATH,"envs/datasets.yaml")
    input:
        expected_download_topgenomes = os.path.join(config["output"], "datasets_topgenomes","expected_download_topgenomes.txt")
    output:
        actual_download_topgenomes = os.path.join(config["output"],"datasets_topgenomes","actual_download_topgenomes.txt")
    params:
        output_dir = os.path.join(config["output"],"datasets_topgenomes")
    shell:
        """
        ls {params.output_dir}/scripts/datasets_*.sh > {params.output_dir}/scripts/script_list.txt
        # This is critical 
        # No // in the file containing paths to the scripts otherwise there would be error!!!
        sed -i 's#//#/#g' {params.output_dir}/scripts/script_list.txt
        # Run the scripts in parallel using GNU parallel
        module load parallel
        parallel --silent --jobs {threads} bash :::: {params.output_dir}/scripts/script_list.txt
        # Get the list of actually downloaded global genomes
        rm -rf {params.output_dir}/scripts
        find {params.output_dir} -name "*.fna" | awk -F"/" '{{print $NF}}' | sed 's/\\.fna$//g' > {output.actual_download_topgenomes}
        """