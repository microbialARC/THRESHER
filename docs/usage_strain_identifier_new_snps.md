# Strain Identifier - New SNPs: Update existing strain/transmission compositions with new genomes using predefined phylothresholds with previous endpoint method.
```
thresher strain_identifier new-snps -h
options:
  -h, --help            show this help message and exit
  --original_metadata ORIGINAL_METADATA
                        Path to the input metadata file containing original genomes.
                        The input file should be a tab-delimited file with 3 or 5 columns.
                        3 columns: genome_name, genome_accession, genome_path
                        5 columns: genome_name, genome_accession, genome_path, patient_id, collection_date
  --new_metadata NEW_METADATA
                        Path to the input metadata file containing new genomes.
                        At least 1 new genome should be provided to perform new-snps analysis.
                        The description of the file format is the same as that for --original_metadata.
  --thresher_output THRESHER_OUTPUT
                        Path to the existing THRESHER strain_identifier directory.
                        The existing analysis directory should contain the previous transmission cluster results.
                        (analysis_mode should be 'full' in the prior run)
  --output OUTPUT       Path to output directory of new-snps analysis.
                        If not provided, defaults to thresher_strain_identifier_new_snps_<YYYY_MM_DD_HHMMSS> under the current working directory.
  --species {sau,sepi,cdiff,kp}
                        Bacteria species.
                        Available options: [sau, sepi, cdiff, kp]
                        sau: Staphylococcus aureus
                        sepi: Staphylococcus epidermidis
                        cdiff: Clostridium difficile
                        kp: Klebsiella pneumoniae
  --snp_coverage_threshold SNP_COVERAGE_THRESHOLD
                        Minimum alignment coverage (0-100) required for pairwise SNP distances to be included in analysis.
                        Low-coverage alignments can yield unreliable SNP count.
                        SNP distance below this threshold are excluded. Default: 80.
  -t THREADS, --threads THREADS
                        Thread number. Default is 1.
  --prefix PREFIX       Prefix for config file. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS
  --conda_prefix CONDA_PREFIX
                        Directory for conda environments needed for this analysis. If not provided, defaults to <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>
```
## Required Input
1. **Original Metadata File(--original_metadata):** 

    A tab-delimited file containing information about the original genomes used in the previous THRESHER analysis.
    - Example metadata file for full mode: 
      
      [Example Input Metadata File](example/example_metadata.txt)
2. **New Metadata File(--new_metadata):** 

    A tab-delimited file containing information about the new genomes to be added to the analysis.
    - Example metadata file for full mode: 
      
      [Example Input Metadata File](example/example_metadata.txt)
    - How do I know the Genbank accession for my genomes? See the [Genbank Accession](genbank_accession.md) page.
3. **Existing THRESHER Output Directory(--thresher_output):** 

    The path to the existing THRESHER Strain Identifier full-pipeline or new-full output directory containing the previous analysis results.
4. **Species(--species):** 

    The bacterial species being analyzed (e.g., sau, sepi, cdiff, kp).

## Optional Input
1. **SNP Coverage Threshold(--snp_coverage_threshold):**
    - Minimum alignment coverage percentage (0-100) required for pairwise SNP distances to be included in analysis (default: 80). Genome pairs with alignment coverage below this threshold are excluded from downstream cladebreaker analysis.
    - Low-coverage alignments can produce artificially low SNP counts, as unaligned regions are not compared and potential variants in those regions go undetected. This can lead to falsely inflated genomic similarity between genomically distantly related genomes. The default threshold of 80% balances sensitivity with reliability.
2. **Output Directory(--output):** 

    Path to the output directory for the new-snps analysis. If not provided, defaults to `thresher_strain_identifier_new_snps_<YYYY_MM_DD_HHMMSS>` under the current working directory.
3. **Thread Number(--threads / -t):**

    Number of threads to use for the analysis. Default is 1. It is highly recommended to increase the thread count using this option to improve performance, regardless of dataset size.
4. **Prefix(--prefix):**

    Prefix for config file. If not provided, defaults to timestamp: `YYYY_MM_DD_HHMMSS`.
5. **Conda Environment Directory(--conda_prefix):**
    Directory for conda environments needed for this analysis. If not provided, defaults to `<OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>`.  You can reuse the conda environments from previous THRESHER runs to save time and disk space.

## Output
1. **Updated Genomes:**
- Plateau: `thresher/output/new_plateau_genomes.csv`
- Peak: `thresher/output/new_peak_genomes.csv`
- Discrepancy: `thresher/output/new_discrepancy_genomes.csv`
- Global: `thresher/output/new_global_genomes.csv`
2. **Updated Strains:**
- Plateau:
  - `thresher/output/new_plateau.RDS` (R object)
  - `thresher/output/new_plateau_strains.csv`
- Peak:
  - `thresher/output/new_peak.RDS` (R object)
  - `thresher/output/new_peak_strains.csv`
- Discrepancy:
  - `thresher/output/new_discrepancy.RDS` (R object)
  - `thresher/output/new_discrepancy_strains.csv`
- Global:
  - `thresher/output/new_global.RDS` (R object)
  - `thresher/output/new_global_strains.csv`
3. **Updated Clusters:**
- RDS: `thresher/output/clusters_summary_new_snps.RDS`
- CSV: `thresher/output/clusters_summary_new_snps.csv`
- Summary table listing each strain/cluster and its composition, including original genomes and newly added genomes: `thresher/output/genomes_summary_new_snps.csv`
  - Column 1 `genome_name`: Name of the genome
  - Column 2 `genome_category`: "original" or "new" indicating whether the genome was part of the original analysis or newly added
  - Column 3 `original_strain`: ID of the strain the genome belongs to in the original analysis ("New-Genome" for new genomes)
  - Column 4 `new_strain`: ID of the strain the genome belongs to in the updated analysis
  - Column 5 `original_cluster`: ID of the cluster the genome belongs to in the original analysis ("New-Cluster" for new genomes, "Non-Cluster" if not belonging to any cluster)
  - Column 6 `new_cluster`: ID of the cluster the genome belongs to in the updated analysis ("Non-Cluster" if not belonging to any cluster)
