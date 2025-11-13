# Resuming Interrupted Runs
Given the pipeline's complexity and potentially long runtime, THRESHER utilizes Snakemake's `--rerun-incomplete` functionality to handle interruptions. If the Strain Identifier pipeline is interrupted, Snakemake can resume from the last successful step without restarting the entire analysis.

To resume:
1. First run the unlock script: `bash <output_directory>/resume/unlock_snakemake_<output_prefix>.sh`
2. Then run the resume script: `bash <output_directory>/resume/resume_snakemake_<output_prefix>.sh`

This approach preserves completed work and saves computational resources.
