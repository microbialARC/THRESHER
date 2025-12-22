rule assembly_scan_new_full_single:
    conda:
        os.path.join(BASE_PATH,"envs/assembly_scan.yaml")
    input:
        genome_path = lambda wc: new_genome_path_dict[wc.genome_name]
    output:
        assembly_scan_output = os.path.join(config["output"], "assembly_scan","{genome_name}_assembly_scan.txt")
    shell:
        """
        assembly-scan --transpose "{input.genome_path}" > {output.assembly_scan_output}
        """

rule assembly_scan_new_full_all:
    input:
        assembly_scan_outputs = expand(os.path.join(config["output"], "assembly_scan","{genome_name}_assembly_scan.txt"), genome_name=list(new_genome_path_dict.keys()))
    output:
        assembly_scan_complete = os.path.join(config["output"], "assembly_scan",".assembly_scan_new_full_complete")
    shell:
        """
        touch {output.assembly_scan_complete}
        """
