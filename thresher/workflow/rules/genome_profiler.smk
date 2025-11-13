rule genome_profiler:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        coords_list = os.path.join(config["output"], "mummer4", "output",f"{genome_name}_coords_list.txt"),
        snps_list = os.path.join(config["output"], "mummer4", "output",f"{genome_name}_snps_list.txt"),
        fna_path = os.path.join(config["output"], "bakta_annotation",genome_name,f"{genome_name}.fna"),
        gff3_path = os.path.join(config["output"], "bakta_annotation",genome_name,f"{genome_name}.gff3"),
    params:
        snp_dir = os.path.join(config["output"], "mummer4","output"),
        output_dir = os.path.join(config["output"], "genome_profiler"),
        genome_name = genome_name,
    threads:
        config["threads"]
    output:
        bin_summary = os.path.join(config["output"], "genome_profiler",f"{genome_name}_bin_summary.csv"),
        chr_bins_rds = os.path.join(config["output"], "genome_profiler",f"{genome_name}_chr_bins.RDS"),
        concat_fasta = os.path.join(config["output"], "genome_profiler",f"{genome_name}_concat.fasta"),
        entropy_csv = os.path.join(config["output"], "genome_profiler",f"{genome_name}_entropy.csv"),
        entropy_rds = os.path.join(config["output"], "genome_profiler",f"{genome_name}_entropy.RDS"),
        mges_seq_fasta = os.path.join(config["output"], "genome_profiler",f"{genome_name}_mges_seq.fasta"),
        mges_csv = os.path.join(config["output"], "genome_profiler",f"{genome_name}_mges.csv"),
        mges_rds = os.path.join(config["output"], "genome_profiler",f"{genome_name}_mges.RDS"),
        new_pos_rds = os.path.join(config["output"], "genome_profiler",f"{genome_name}_new_pos.RDS"),
        plot_df_rds = os.path.join(config["output"], "genome_profiler",f"{genome_name}_plot_df.RDS"),
        position_coverage_rds = os.path.join(config["output"], "genome_profiler",f"{genome_name}_position_coverage.RDS"),
        profiler_plot_pdf = os.path.join(config["output"], "genome_profiler",f"{genome_name}_profiler_plot.pdf"),
        snps_sum_rds = os.path.join(config["output"], "genome_profiler",f"{genome_name}_snps_sum.RDS")
    script:
        # The genome_profiler.R script is a command-line R script that accepts arguments.
        # It can be run independently of Snakemake as a standalone tool.
        # When used as standalone tool, it performs checks and will skip steps if the expected output files already exist.
        # Command-line Rscript tool. Accepts these flags when run standalone:
        #   --input_genome <path>   : reference assembly (FASTA; .fasta/.fna). Required.
        #   --gff <path>            : annotation file (GFF3). Required.
        #   --snp <dir>             : directory with SNP/coord outputs (expects *.snps and *.coords from MUMmer). Required.
        #   --cpus <int>            : number of CPU cores to use (optional; defaults to detectCores()).
        #   --output <dir>          : directory where all outputs (RDS/CSV/PDF) will be written. Required.
        #
        # Snakemake integration:
        #   When run inside a Snakemake rule, Snakemake should invoke Rscript with the same flags
        #   (see the pipeline's genome_profiler.smk rule). The script reads commandArgs(trailingOnly=TRUE).
        #
        # Outputs (written to --output):
        #   *_concat.fasta, *_new_pos.RDS, *_position_coverage.RDS, *_snps_sum.RDS, *_chr_bins.RDS,
        #   *_entropy.RDS / .csv, *_mges.RDS / .csv / _mges_seq.fasta, *_bin_summary.csv, *_profiler_plot.pdf
        os.path.join(BASE_PATH,"scripts","genome_profiler.R")