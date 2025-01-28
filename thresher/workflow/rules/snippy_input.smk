rule snippy_input:
    input:
        hc_groups_csv = os.path.join(config["output"], "thresher", "input", "hierarchical_clustering_groups_simplified.csv"),
        whatsgnu_results = expand(os.path.join(config["output"], "whatsgnu","{genome_name}","{genome_name}_WhatsGNU_topgenomes.txt"),genome_name=genome_path_dict.keys()),
        assembly_scan_results = expand(os.path.join(config["output"], "assembly_scan","{genome_name}_assembly_scan.txt"),genome_name=genome_path_dict.keys())
    output:
        snippy_multi = os.path.join(config["output"], "snippy","scripts","snippy_multi.sh"),
        snippy_reference = os.path.join(config["output"], "snippy","tabs","snippy_reference.txt")
    params:
        script_dir = os.path.join(config["output"], "snippy","scripts"),
        tab_dir = os.path.join(config["output"], "snippy","tabs"),
        genome_path = genome_path_dict,
        global_genome_path = os.path.join(config["output"],"datasets_topgenomes"),
        study_accession = genome_accession
    threads:
        config["threads"]
    resources:
        mem_mb = int(config["threads"])*1000
    script:
        os.path.join(BASE_PATH,"scripts","snippy_input.py")
