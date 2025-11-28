rule mummer4_global_cmd:
    input:
        whatsgnu_output = [os.path.join(config["output"],"whatsgnu",f"{genome}",f"{genome}_WhatsGNU_topgenomes.txt") for genome in list(genome_path_dict.keys())],
        actual_download_topgenomes = os.path.join(config["output"],"datasets_topgenomes","actual_download_topgenomes.txt")
    output:
        [os.path.join(config["output"],"mummer4_global","scripts",f"dnadiff_{genome}.sh") for genome in list(genome_path_dict.keys())]
    params:
        study_accession = genome_accession,
        whatsgnu_output_dir = os.path.join(config["output"],"whatsgnu"),
        global_genome_dir = os.path.join(config["output"],"datasets_topgenomes"),
        study_genome_dict = genome_path_dict,
        script_path = os.path.join(config["output"],"mummer4_global","scripts"),
        output_path = os.path.join(config["output"],"mummer4_global")
    script:
        os.path.join(BASE_PATH,"scripts","sort_whatsgnu_output.py")

rule mummer4_global:
    conda:
        os.path.join(BASE_PATH,"envs/mummer4.yaml")
    input:
        script = [os.path.join(config["output"],"mummer4_global","scripts",f"dnadiff_{genome}.sh") for genome in list(genome_path_dict.keys())]
    output:
        reports = expand(os.path.join(config["output"], "mummer4_global","{genome_name}_concatenated.report"), genome_name=genome_path_dict.keys())
    params:
        output_dir = config["output"]
    threads:
        config["threads"]
    shell:
        """
        chmod +x {params.output_dir}/mummer4_global/scripts/dnadiff*.sh
        ls {params.output_dir}/mummer4_global/scripts/dnadiff_*.sh > {params.output_dir}/mummer4_global/scripts/script_list.txt
        sed -i 's#//#/#g' {params.output_dir}/mummer4_global/scripts/script_list.txt
        module load parallel
        parallel --jobs {threads} bash :::: {params.output_dir}/mummer4_global/scripts/script_list.txt
        rm {params.output_dir}/mummer4_global/scripts/script_list.txt
        """