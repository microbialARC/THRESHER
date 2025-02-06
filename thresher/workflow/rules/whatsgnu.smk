rule whatsgnu:
    conda:
        os.path.join(BASE_PATH,"envs/whatsgnu.yaml")
    input: 
        faa = [os.path.join(config["output"],"bakta_annotation",f"{genome}",f"{genome}.faa") for genome in list(genome_path_dict.keys())]
    params:
        species = config["species"],
        genome_names = list(genome_path_dict.keys()),
        db_path = config["whatsgnu_db_path"],
        output_dir = config["output"]
    output:
        whatsgnu_output = [os.path.join(config["output"],"whatsgnu",f"{genome}",f"{genome}_WhatsGNU_topgenomes.txt") for genome in list(genome_path_dict.keys())]
    shell:
        """
        genome_names=({params.genome_names})
        input_faa=({input.faa}) 
        mkdir -p {params.output_dir}/whatsgnu/
        db_path=""
        if [ "{params.db_path}" == "None" ]; then
            mkdir -p {params.output_dir}/whatsgnu/db
            if [ "{params.species}" == "sau" ]; then
                wget -O {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz https://zenodo.org/records/13388052/files/Saureus.tar.gz?download=1 
                tar -xvzf {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz -C {params.output_dir}/whatsgnu/db
                db_path="{params.output_dir}/whatsgnu/db/Saureus/Sau_WhatsGNU_Ortholog_db.pickle"
                whatsgnu_mode="ortholog"
            elif [ "{params.species}" == "cdiff" ]; then
                wget -O {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz https://zenodo.org/records/13387715/files/Cdiff.tar.gz?download=1
                tar -xvzf {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz -C {params.output_dir}/whatsgnu/db
                db_path="{params.output_dir}/whatsgnu/db/Cdiff/Cdiff_WhatsGNU_Ortholog_db.pickle"
                whatsgnu_mode="ortholog"
            elif [ "{params.species}" == "kp" ]; then
                wget -O {params.output_dir}/whatsgnu/db/whatsgnu_db.zip https://zenodo.org/record/7812697/files/Kp.zip?download=1
                unzip {params.output_dir}/whatsgnu/db/whatsgnu_db.zip -d {params.output_dir}/whatsgnu/db
                db_path="{params.output_dir}/whatsgnu/db/Kp_Ortholog_8752.pickle"
                whatsgnu_mode="ortholog"
            elif [ "{params.species}" == "sepi" ]; then
                wget -O {params.output_dir}/whatsgnu/db/whatsgnu_db.gz https://zenodo.org/records/14751549/files/Sepi_WhatsGNU_basic.txt.gz?download=1
                mkdir -p {params.output_dir}/whatsgnu/db
                gunzip -c {params.output_dir}/whatsgnu/db/whatsgnu_db.gz > {params.output_dir}/whatsgnu/db/Sepi_WhatsGNU_basic.txt
                db_path="{params.output_dir}/whatsgnu/db/Sepi_WhatsGNU_basic.txt"
                whatsgnu_mode="basic"
            fi
        elif [ "{params.db_path}" != "None" ]; then
            db_path="{params.db_path}"
            if [ "{params.species}" == "sepi" ]; then
                whatsgnu_mode="basic"
            else
                whatsgnu_mode="ortholog"
            fi
        fi
        
        cd {params.output_dir}/whatsgnu/

        for ((i=0; i<${{#genome_names[@]}}; i++)); do
            WhatsGNU_main.py \
            -d "${{db_path}}" \
            -dm "${{whatsgnu_mode}}" \
            --force \
            -t \
            -tn 20 \
            -o "${{genome_names[$i]}}" \
            "${{input_faa[$i]}}"
        done
        """