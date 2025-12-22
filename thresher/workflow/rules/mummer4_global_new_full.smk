rule mummer4_global_new_full_cmd:
    input:
        whatsgnu_output = [os.path.join(config["output"],"whatsgnu","whatsgnu_results",f"{genome}_WhatsGNU_topgenomes.txt") for genome in list(new_genome_path_dict.keys())],
        actual_download_topgenomes = os.path.join(config["output"],"datasets_topgenomes","actual_download_topgenomes.txt")
    output:
        [os.path.join(config["output"],"mummer4_global","scripts",f"dnadiff_{genome}.sh") for genome in list(new_genome_path_dict.keys())]
    params:
        study_accession = new_genome_accession,
        whatsgnu_output_dir = os.path.join(config["output"],"whatsgnu"),
        global_genome_dir = os.path.join(config["output"],"datasets_topgenomes"),
        study_genome_dict = new_genome_path_dict,
        script_path = os.path.join(config["output"],"mummer4_global","scripts"),
        output_path = os.path.join(config["output"],"mummer4_global")
    script:
        os.path.join(BASE_PATH,"scripts","sort_whatsgnu_output.py")

rule mummer4_global_new_full_single:
    conda:
        os.path.join(BASE_PATH,"envs/mummer4.yaml")
    input:
        script = os.path.join(config["output"], "mummer4_global", "scripts", "dnadiff_{genome_name}.sh")
    output:
        reports = os.path.join(config["output"], "mummer4_global","{genome_name}_concatenated.report")
    params:
        output_dir = config["output"]
    shell:
        """
        chmod +x {input.script}
        bash {input.script}
        """

rule mummer4_global_new_full_all:
    input:
        reports = expand(os.path.join(config["output"], "mummer4_global","{genome_name}_concatenated.report"), genome_name=new_genome_path_dict.keys())
    output:
        done = os.path.join(config["output"], "mummer4_global", ".mummer4_global_complete")
    shell:
        """
        touch {output.done}
        """
