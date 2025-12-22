rule mlst_output_single:
    conda:
        os.path.join(BASE_PATH,"envs/mlst.yaml")
    input:
        genome_path = lambda wc: new_genome_path_dict[wc.genome_name]
    output:
        mlst_output_single = os.path.join(config["output"], "mlst","raw","{genome_name}_mlst.csv")
    params:
        output_dir=os.path.join(config["output"], "mlst","raw")
    shell:
        """
        mlst {input.genome_path} > {output.mlst_output_single}
        """

rule mlst_output_all:
    input:
        mlst_outputs = expand(os.path.join(config["output"], "mlst","raw","{genome_name}_mlst.csv"), genome_name=list(new_genome_path_dict.keys()))
    output:
        mlst_output = os.path.join(config["output"], "mlst","raw","mlst_raw.csv")
    shell:
        """
        cat {input.mlst_outputs} > {output.mlst_output}
        """

rule mlst_results: 
    input:
        mlst_output = os.path.join(config["output"], "mlst","raw","mlst_raw.csv")
    output:
        mlst_results = os.path.join(config["output"], "mlst","summary","mlst_results.csv")
    params:
        species = config["species"],
        analysis_mode = config["mode"],
        metadata = config["new_metadata"],
        mlst_sau_db = os.path.join(BASE_PATH,"db/mlst/Sau_mlst.txt"),
        mlst_cdiff_db = os.path.join(BASE_PATH,"db/mlst/Cdiff_mlst.txt"),
        mlst_kp_db = os.path.join(BASE_PATH,"db/mlst/Kp_mlst.txt")
    script:
        os.path.join(BASE_PATH,"scripts","mlst_results.py")