rule bakta_db:
    conda:
        os.path.join(BASE_PATH,"envs/bakta.yaml")
    params:
        output_dir = config["output"],
        bakta_db_type = config["bakta_db_type"],
        bakta_db_path = config["bakta_db_path"],
    output:
        os.path.join(config["output"],"bakta_db","bakta.db")
    shell:
        """
        set -euo pipefail

        output_dir="{params.output_dir}"
        bakta_db_dir="{params.bakta_db_path}"
        bakta_db_type="{params.bakta_db_type}"

        # If bakta_db_path is set to "None", use output_dir/bakta_db
        if [ "${{bakta_db_dir}}" = "None" ]; then
            echo "No bakta_db_path provided, using default path: ${{output_dir}}/bakta_db"
            bakta_db_dir="${{output_dir}}/bakta_db"
            echo "Creating the directory for ${{bakta_db_type}} database: ${{bakta_db_dir}}"
            mkdir -p "${{bakta_db_dir}}"

            echo "Downloading ${{bakta_db_type}} database to ${{bakta_db_dir}}"
            bakta_db download --output "${{bakta_db_dir}}" --type "${{bakta_db_type}}"

            # Move files up and clean folders
            
            if [ -d "${{bakta_db_dir}}/db" ] || [ -d "${{bakta_db_dir}}/db-light" ]; then
                mv "${{bakta_db_dir}}"/*/* "${{bakta_db_dir}}" || true
                rm -rf "${{bakta_db_dir}}/db" "${{bakta_db_dir}}/db-light"
                echo "${{bakta_db_type}} database downloaded to ${{bakta_db_dir}}"
            fi
        else
            echo "Using provided bakta_db_path: ${{bakta_db_dir}}"
        fi

        """
