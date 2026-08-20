rule bakta_db:
    conda:
        os.path.join(BASE_PATH,"envs/bakta.yaml")
    params:
        output_dir = config["output"],
        bakta_db_type = config["bakta_db_type"],
        bakta_db_path = config["bakta_db_path"],
        # Hardcoded AMRFinderPlus database version to ensure compatibility
        amrfinder_db_version = "2025-07-16.1",
        amrfinder_db_format = "4.0"
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

        # AMRFinderPlus database version check

        amrfinder_db_root="${{bakta_db_dir}}/amrfinderplus-db"
        required_db="{params.amrfinder_db_version}"
        db_format="{params.amrfinder_db_format}"
        ftp_url="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/${{db_format}}/${{required_db}}/"

        if [ -f "${{amrfinder_db_root}}/latest/version.txt" ]; then
            current_db="$(tr -d '[:space:]' < "${{amrfinder_db_root}}/latest/version.txt")"
        else
            current_db="none"
        fi
        
        # Print the current and required database versions
        echo "Current AMRFinderPlus database version: ${{current_db}}"
        echo "Required AMRFinderPlus database version: ${{required_db}}"

        if [ "${{current_db}}" = "${{required_db}}" ]; then
            echo "AMRFinderPlus database OK."
        else
            echo "Version mismatch, fetching ${{required_db}} from ${{ftp_url}}"
            target_dir="${{amrfinder_db_root}}/${{required_db}}"
            rm -rf "${{target_dir}}"
            mkdir -p "${{target_dir}}"
            wget --recursive --no-parent --no-directories \
                 --no-verbose --tries=3 --timeout=60 \
                 --execute robots=off \
                 --reject "index.html*" \
                 --directory-prefix "${{target_dir}}" \
                 "${{ftp_url}}"
            
            # NCBI ships flat data files only. So we need to locally build BLAST/HMMER indexes
            echo "Building BLAST/HMMER indexes for ${{target_dir}}"
            amrfinder_index "${{target_dir}}"

            # Repoint 'latest' at the pinned release
            if [ -d "${{amrfinder_db_root}}/latest" ] && [ ! -L "${{amrfinder_db_root}}/latest" ]; then
                rm -rf "${{amrfinder_db_root}}/latest"
            fi
            ln -sfn "${{required_db}}" "${{amrfinder_db_root}}/latest"

            if [ -f "${{amrfinder_db_root}}/latest/version.txt" ]; then
                current_db="$(tr -d '[:space:]' < "${{amrfinder_db_root}}/latest/version.txt")"
            else
                current_db="none"
            fi

            if [ "${{current_db}}" != "${{required_db}}" ]; then
                echo "ERROR: AMRFinderPlus database is '${{current_db}}' after download, expected '${{required_db}}'" >&2
                exit 1
            fi
            echo "AMRFinderPlus database pinned at ${{current_db}}"
        fi

        """
