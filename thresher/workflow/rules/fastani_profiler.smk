rule fastani_profiler_raw:
    conda:
        os.path.join(BASE_PATH,"envs/fastani.yaml")
    input:
        genome_path = input_genome_path,
        topgenomes_result = os.path.join(config["output"], "whatsgnu", genome_name, f"{genome_name}_WhatsGNU_topgenomes.txt"),
        result_check = os.path.join(config["output"],"datasets_topgenomes","result_check.txt")
    threads:
        config["threads"]
    output:
        topgenomes_all = os.path.join(config["output"], "fastani",f"{genome_name}_topgenomes_all.txt"),
        fastani_result = os.path.join(config["output"], "fastani",f"{genome_name}_fastani.csv"),
        topgenomes_filtered = os.path.join(config["output"], "fastani",f"{genome_name}_topgenomes_filtered.txt")
    params:
        ani_threshold = config["ani_threshold"],
        topgenomes_dir = os.path.join(config["output"], "datasets_topgenomes"),
        output_dir = os.path.join(config["output"], "fastani")
    shell:
        """
        # Import from Snakemake config
        # Input
        output_dir="{params.output_dir}"
        reference_genome_path="{input.genome_path}"
        topgenomes_dir="{params.topgenomes_dir}"
        topgenomes_result_path="{input.topgenomes_result}"
        ani_threshold={params.ani_threshold}

        # Output
        topgenomes_all_path="{output.topgenomes_all}"
        topgenomes_filtered_path="{output.topgenomes_filtered}"
        fastani_result="{output.fastani_result}"


        # Create output directory if it doesn't exist
        mkdir -p "$output_dir"

        # Read the input.result_check to confirm all global genomes have been downloaded
        # If the file content is "All global genomes downloaded", proceed
        check_message=$(cat "{input.result_check}")
        if [ "$check_message" != "All global genomes downloaded" ]; then
            echo "Error: Global genomes not downloaded properly."
            exit 1
        fi
        
        # Parse WhatsGNU output to obtain top genome accession and prepend the genomes directory,
        # yielding absolute file paths for use in fastANI

        tail -n +2 "$topgenomes_result_path" | while read -r genome_entry _; do
            echo "$topgenomes_dir/${{genome_entry}}.fna"
        done > $topgenomes_all_path

        # Run FastANI
        fastANI -q "$reference_genome_path" \
                --rl "$topgenomes_all_path" \
                --threads {threads} \
                --output "$fastani_result"

        # Filter FastANI result based on ANI threshold
        # Column 3 in FastANI output contains the ANI values
        # Only return the Column 2 (top genome path) to the filtered output file
        awk -v threshold="$ani_threshold" '$3 >= threshold {{print $2}}' "$fastani_result" > "$topgenomes_filtered_path"
        """
