# Genome Profiler
Infer substitution probabilities and mobile genetic element (MGE) dynamics leveraging publicly available genomes.
```
thresher genome_profiler -h
options:
  -h, --help            show this help message and exit
  --input_genome INPUT_GENOME
                        Path to genome assembly in FASTA format
  --output OUTPUT       Output directory of the profiling results. If not provided, defaults to thresher_genome_profiler_output_<YYYY_MM_DD_HHMMSS> under the current working directory.
  --species {sau,sepi,cdiff,kp}
                        Bacteria species.
                        Available options: [sau, sepi, cdiff, kp]
                        sau: Staphylococcus aureus
                        sepi: Staphylococcus epidermidis
                        cdiff: Clostridium difficile
                        kp: Klebsiella pneumoniae
  --top_genomes TOP_GENOMES
                        Number of initial top genomes for profiling before applying the ANI exclusion threshold (default: 1000).
  --ani_threshold ANI_THRESHOLD
                        Average Nucleotide Identity (ANI) exclusion threshold (default: 99.5).
                        Genomes with ANI below this value will be removed from the profiling analysis.
  --bakta_db_type BAKTA_DB_TYPE
                        Bakta database.
                        Available options: [full, light]
                        Default is full
  --bakta_db_path BAKTA_DB_PATH
                        The path of the directory where the existing Bakta database locates.
                        If provided, the Bakta database will not be downloaded.
                        If not provided, defaults to <OUTPUT>/bakta/db
  --whatsgnu_db_path WHATSGNU_DB_PATH
                        The path to the existing WhatsGNU database.
                        If provided, the WhatsGNU database will not be downloaded.
                        If not provided, defaults to <OUTPUT>/whatsgnu/db.
  -t THREADS, --threads THREADS
                        Thread number. Default is the maximum available.
  --prefix PREFIX       Prefix for config file, output files, and analysis naming. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS
  --conda_prefix CONDA_PREFIX
                        Directory for conda environments needed for this analysis. If not provided, defaults to OUTPUT/conda_envs_<YYYY_MM_DD_HHMMSS>
  --force               Bypass system compatibility checks (operating system and available RAM) and force execution of the pipeline.
                        This may cause instability or failures.
```

## Required Input
1. **Input Genome:**
    
    Path to the genome assembly in FASTA format for profiling.

2. **Species:**

    Indicate the species being analyzed. Currently supported species include:
      - `sau` (*Staphylococcus aureus*)
      - `sepi` (*Staphylococcus epidermidis*)
      - `cdiff` (*Clostridium difficile*)
      - `kp` (*Klebsiella pneumoniae*)

## Optional Input
1. **Output Directory:**
    - Path to the output directory. If not provided, defaults to `thresher_genome_profiler_output_<YYYY_MM_DD_HHMMSS>` under the current working directory.
2. **Top Genomes:**
    - Number of initial top genomes for profiling before applying the ANI exclusion threshold (default: 1000).
3. **ANI Exclusion Threshold:**
    - Average Nucleotide Identity (ANI) exclusion threshold (default: 99.5).
    - Genomes with ANI below this value will be removed from the profiling analysis.
4. **Bakta Database Type and Path:**
    - Specify the type of Bakta database (`full` or `light`, default is `full`).
    - Path to an existing Bakta database. If not provided, the database will be downloaded to `<OUTPUT>/bakta/db`.
5. **WhatsGNU Database Path:**

    Path to an existing WhatsGNU database. If not provided, the database will be downloaded to `<OUTPUT>/whatsgnu/db`.
6. **Threads:**

    Number of threads to use (default is the maximum available).
7. **Prefix:**

    Prefix for config file, output files, and analysis naming. If not provided, defaults to a timestamp in the format `YYYY_MM_DD_HHMMSS`.
8. **Conda Prefix:**

    Directory for conda environments needed for this analysis. If not provided, defaults to `<OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>`.
9. **Force Execution:**

    Bypass system compatibility checks (operating system and available RAM) and force execution of the pipeline. This may cause instability or failures.  

## Output
1. **Bakta Annotation:**
  - `bakta_annotation/{genome_name}/`
2. **WhatsGNU:**
  - `whatsgnu/{genome_name}/`
3. **FastANI:**
 - Path to all query genomes used in FastANI:
 
    `fastani/{genome_name}_topgenomes_all.txt`
 - FastANI results between the input genome and all query genomes by top_genomes option:
 
    `fastani/{genome_name}_fastani.csv `
 - Path to filtered query genomes after applying ANI exclusion threshold:

    `fastani/{genome_name}_topgenomes_filtered.txt`

4. **Pair-wise comparison using Mummer4:**
 - Scripts for running Mummer4 pair-wise comparison: `mummer4/scripts/`
 - Raw Mummer4 Output: `mummer4/output/`

5. **Global genomes downloaded using Datasets:**
 - `datasets_topgenomes/*.fna`

6. **Genome Profiler Results:**
 - Bin summary of concatenated input genome:

    `genome_profiler/{genome_name}_bin_summary.csv`
 - Summary of MGEs identified in the input genome:

    `genome_profiler/{genome_name}_mges.csv`
 - Plots indicating the MGEs and substitution probabilities across the genome:

    `genome_profiler/{genome_name}_profiler_plot.pdf`

    `genome_profiler/{genome_name}_profiler_plot.RDS` (R object)

- Substitution probability table across each site within the genome:

    `genome_profiler/{genome_name}_entropy.csv`

    `genome_profiler/{genome_name}_entropy.RDS` (R object)