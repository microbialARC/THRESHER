if config.get("bakta_db_path") == "None":
    config["bakta_db_path"] = os.path.join(config["output"],"bakta_db/")
    rule bakta_db:
            conda:
                os.path.join(BASE_PATH,"envs/bakta.yaml")
            params:
                bakta_db_type = config["bakta_db"]
            output:
                os.path.join(config["bakta_db_path"],"bakta.db")
            shell:
                """
                if [ "{params.bakta_db_type}" != "full" ] && [ "{params.bakta_db_type}" != "light" ]; then
                    echo "Error: Unknown bakta_db type: {params.bakta_db_type}" >&2
                    exit 1
                fi

                echo "Creating the directory for database"
                mkdir -p output/bakta_db/

                echo "Downloading database"
                bakta_db download --output output/bakta_db/ --type {params.bakta_db_type}
                """
