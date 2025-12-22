rule assembly_scan_single:
    conda:
        os.path.join(BASE_PATH,"envs/assembly_scan.yaml")
    input:
        genome_path = lambda wc: genome_path_dict[wc.genome_name]
    output:
        assembly_scan_single = os.path.join(config["output"], "assembly_scan","{genome_name}_assembly_scan.txt")
    shell:
        """
        assembly-scan --transpose {input.genome_path} > {output.assembly_scan_single}
        """

rule assembly_scan_all:
    input:
        assembly_scan_outputs = expand(os.path.join(config["output"], "assembly_scan","{genome_name}_assembly_scan.txt"), genome_name=list(genome_path_dict.keys()))
    output:
        assembly_scan_complete = os.path.join(config["output"], "assembly_scan",".assembly_scan_complete")
    shell:
        """
        touch {output.assembly_scan_complete}
        """
