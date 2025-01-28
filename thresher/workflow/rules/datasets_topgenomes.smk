rule dataset_topgenomes:
    conda:
        os.path.join(BASE_PATH,"envs/datasets.yaml")
    input:
        topgenomes = os.path.join(config["output"], "datasets_topgenomes","topgenomes.txt")
    output:
        result_check = os.path.join(config["output"],"datasets_topgenomes","result_check.txt")
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
        parallel --jobs {threads} bash :::: {params.output_dir}/scripts/script_list.txt
        rm -rf {params.output_dir}/scripts
        echo "All global genomes downloaded" > {output.result_check}
        """
