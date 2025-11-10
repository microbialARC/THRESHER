rule bakta_db:
    conda:
        os.path.join(BASE_PATH,"envs/bakta.yaml")
    params:
        bakta_db_type = config["bakta_db_type"],
        bakta_db_path = config["bakta_db_path"],
    output:
        os.path.join(config["output"],"bakta_db","bakta.db")
    shell:
        """
        # If the bakta database already exists, skip the download
        if [ -e "{params.bakta_db_path}/bakta.db" ]; then
            echo "Bakta database already exists, skipping download"
        else
            # Download the bakta database if not exists
            if [ "{params.bakta_db_type}" != "full" ] && [ "{params.bakta_db_type}" != "light" ]; then
                echo "Error: Unknown bakta_db type: {params.bakta_db_type}" >&2
                # If the bakta_db_type is not valid, continue with default "full" version of bakta_db
                params.bakta_db_type="full"
            fi

            echo "Creating the directory for {params.bakta_db_type} database"
            mkdir -p {params.bakta_db_path}

            echo "Downloading {params.bakta_db_type} database"
            bakta_db download --output {params.bakta_db_path} --type {params.bakta_db_type}
            mv {params.bakta_db_path}/*/* {params.bakta_db_path}
            rm -rf {params.bakta_db_path}/db/
            rm -rf {params.bakta_db_path}/light-db/
            echo "{params.bakta_db_type} database downloaded"
        """
