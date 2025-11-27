rule dataset_topgenomes_profiler:
    conda:
        os.path.join(BASE_PATH,"envs/datasets.yaml")
    input:
        topgenomes_result = os.path.join(config["output"], "whatsgnu", genome_name, f"{genome_name}_WhatsGNU_topgenomes.txt")
    output:
        actual_download_topgenomes = os.path.join(config["output"],"datasets_topgenomes","actual_download_topgenomes.txt")
    params:
        output_dir = os.path.join(config["output"],"datasets_topgenomes")
    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p {params.output_dir}/scripts
        topgenomes_file="{input.topgenomes_result}"

        # Read the topgenomes from WhatsGNU output
        # Compared to the rules and scripts used in Strain Identifer
        # Instead of using a python script,
        # this rule directly uses shell commands to create the command for downloading top genomes determined by WhatsGNU
        # Skip the first row(tail -n +2) of the whatsgnu topgenomes file and extract the genome accessions at column 1, separated by tab(cut -f1)
        topgenomes=$(tail -n +2 $topgenomes_file | cut -f1)
        # For each genome accession, create a script to download the genome using datasets command
        while IFS= read -r genome_entry; do
            # Create a script to download each genome using datasets command
            echo "#!/bin/bash
            datasets download genome accession ${{genome_entry}} --filename {params.output_dir}/${{genome_entry}}.zip
            unzip -o {params.output_dir}/${{genome_entry}}.zip -d {params.output_dir}/${{genome_entry}}/
            mv {params.output_dir}/${{genome_entry}}/ncbi_dataset/data/*/*.fna {params.output_dir}/${{genome_entry}}.fna
            rm {params.output_dir}/${{genome_entry}}.zip
            rm -rf {params.output_dir}/${{genome_entry}}" > {params.output_dir}/scripts/datasets_${{genome_entry}}.sh
        done < <(echo "$topgenomes")

        ls {params.output_dir}/scripts/datasets_*.sh > {params.output_dir}/scripts/script_list.txt
        # This is critical 
        # No // in the file containing paths to the scripts otherwise there would be error!!!
        sed -i 's#//#/#g' {params.output_dir}/scripts/script_list.txt
        # Run the scripts in parallel using GNU parallel
        module load parallel
        parallel --jobs {threads} bash :::: {params.output_dir}/scripts/script_list.txt
        rm -rf {params.output_dir}/scripts
        find {params.output_dir} -name "*.fna" | awk -F"/" '{{print $NF}}' | sed 's/\\.fna$//g' > {output.actual_download_topgenomes}
        """
