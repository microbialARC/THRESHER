rule snippy_input_new_full:
    input:
        hc_groups_csv = os.path.join(config["output"], "thresher", "input", "hierarchical_clustering_groups_simplified.csv"),
        new_whatsgnu_results = expand(os.path.join(config["output"], "whatsgnu","whatsgnu_results","{genome_name}_WhatsGNU_topgenomes.txt"),genome_name=new_genome_path_dict.keys()),
        new_assembly_scan_results = expand(os.path.join(config["output"], "assembly_scan","{genome_name}_assembly_scan.txt"),genome_name=new_genome_path_dict.keys()),
        actual_download_topgenomes = os.path.join(config["output"],"datasets_topgenomes","actual_download_topgenomes.txt")
    output:
        snippy_multi = os.path.join(config["output"], "snippy","scripts","snippy_multi.sh"),
        snippy_reference = os.path.join(config["output"], "snippy","tabs","snippy_reference.txt")
    params:
        original_whatsgnu_results = expand(os.path.join(config["thresher_output"], "whatsgnu","whatsgnu_results","{genome_name}_WhatsGNU_topgenomes.txt"),genome_name=original_genome_path_dict.keys()),
        original_assembly_scan_results = expand(os.path.join(config["thresher_output"], "assembly_scan","{genome_name}_assembly_scan.txt"),genome_name=original_genome_path_dict.keys()),
        script_dir = os.path.join(config["output"], "snippy","scripts"),
        tab_dir = os.path.join(config["output"], "snippy","tabs"),
        new_genome_path = new_genome_path_dict,
        original_genome_path = original_genome_path_dict,
        new_global_genome_path = os.path.join(config["output"],"datasets_topgenomes"),
        original_global_genome_path = os.path.join(config["thresher_output"],"datasets_topgenomes"),
        new_study_accession = new_genome_accession,
        original_study_accession = original_genome_accession
    threads:
        config["threads"]
    script:
        os.path.join(BASE_PATH,"scripts","snippy_input_new_full.py")
