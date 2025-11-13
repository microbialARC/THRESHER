rule whatsgnu_profiler:
    conda:
        os.path.join(BASE_PATH,"envs/whatsgnu.yaml")
    input: 
        faa = os.path.join(config["output"], "bakta_annotation",genome_name,f"{genome_name}.faa")
    params:
        species = config["species"],
        genome_name = genome_name,
        db_path = config["whatsgnu_db_path"],
        output_dir = config["output"],
        top_genome = config["top_genomes"]
    output:
        whatsgnu_output = os.path.join(config["output"], "whatsgnu", genome_name, f"{genome_name}_WhatsGNU_topgenomes.txt")
    shell:
        """
        genome_name="{params.genome_name}"
        input_faa_path="{input.faa}"
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

        # Import the top genomes number 
        top_genomes="{params.top_genome}"

        WhatsGNU_main.py \
        -d "${{db_path}}" \
        -dm "${{whatsgnu_mode}}" \
        --force \
        -t \
        -tn "${{top_genomes}}" \
        -o "${{genome_name}}" \
        "${{input_faa_path}}"
        """