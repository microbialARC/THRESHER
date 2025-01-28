rule mlst_output:
    conda:
        os.path.join(BASE_PATH,"envs/mlst.yaml")
    input:
        genome_paths = list(genome_path_dict.values())
    output:
        mlst_output = os.path.join(config["output"], "mlst","mlst_raw.csv")
    params:
        output_dir=os.path.join(config["output"], "mlst")
    shell:
        """
        mkdir -p {params.output_dir}
        for genome in {input.genome_paths}; do
            mlst --csv $genome > {params.output_dir}/$(basename $genome)_mlst.csv
        done
        cat {params.output_dir}/*_mlst.csv > {output.mlst_output}
        """

rule mlst_results: 
    input:
        mlst_output = os.path.join(config["output"], "mlst","mlst_raw.csv")
    output:
        mlst_results = os.path.join(config["output"], "mlst","mlst_results.csv")
    params:
        species = config["species"],
        mlst_sau_db = os.path.join(BASE_PATH,"db/mlst/Sau_mlst.txt"),
        mlst_cdiff_db = os.path.join(BASE_PATH,"db/mlst/Cdiff_mlst.txt"),
        mlst_kp_db = os.path.join(BASE_PATH,"db/mlst/Kp_mlst.txt")
    script:
        os.path.join(BASE_PATH,"scripts","mlst_results.py")