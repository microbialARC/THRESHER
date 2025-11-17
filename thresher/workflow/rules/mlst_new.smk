rule mlst_output_new:
    conda:
        os.path.join(BASE_PATH,"envs/mlst.yaml")
    input:
        new_genome_paths = list(new_genome_path_dict.values())
    output:
        mlst_output = os.path.join(config["output"], "mlst","raw","mlst_raw.csv")
    params:
        output_dir=os.path.join(config["output"], "mlst","raw")
    shell:
        """
        mkdir -p {params.output_dir}
        for genome in {input.new_genome_paths}; do
            mlst --csv $genome > {params.output_dir}/$(basename $genome)_mlst.csv
        done
        cat {params.output_dir}/*_mlst.csv > {output.mlst_output}
        """

rule mlst_results_new: 
    input:
        mlst_output = os.path.join(config["output"], "mlst","raw","mlst_raw.csv")
    output:
        mlst_results = os.path.join(config["output"], "mlst","summary","mlst_results.csv")
    params:
        species = config["species"],
        metadata = config["new_metadata"],
        mlst_sau_db = os.path.join(BASE_PATH,"db/mlst/Sau_mlst.txt"),
        mlst_cdiff_db = os.path.join(BASE_PATH,"db/mlst/Cdiff_mlst.txt"),
        mlst_kp_db = os.path.join(BASE_PATH,"db/mlst/Kp_mlst.txt")
    script:
        os.path.join(BASE_PATH,"scripts","mlst_results.py")