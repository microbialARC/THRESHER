# Strain Identifier - New SNPs: Rerun the full pipeline to update the phylothresholds, and update strain/transmission compositions with new genomes.
```
thresher strain_identifier new-full -h
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
  --output OUTPUT       Path to output directory of new-full analysis.
                        If not provided, defaults to thresher_strain_identifier_new_full_<YYYY_MM_DD_HHMMSS> under the current working directory.
  --species {sau,sepi,cdiff,kp}
                        Bacteria species.
                        Available options: [sau, sepi, cdiff, kp]
                        sau: Staphylococcus aureus
                        sepi: Staphylococcus epidermidis
                        cdiff: Clostridium difficile
                        kp: Klebsiella pneumoniae
  --whatsgnu_db_path WHATSGNU_DB_PATH
                        The path to the existing WhatsGNU database.
                        If provided, the WhatsGNU database will not be downloaded.
                        If not provided, defaults to <OUTPUT>/whatsgnu/db
  --bakta_db_type BAKTA_DB_TYPE
                        The type of Bakta database.
                        Available options: [full, light]
                        Default is full
  --bakta_db_path BAKTA_DB_PATH
                        The path of the directory where the existing Bakta database locates.
                        If provided, the Bakta database will not be downloaded.
                        If not provided, defaults to <OUTPUT>/bakta/db
  --snp_coverage_threshold SNP_COVERAGE_THRESHOLD
                        Minimum alignment coverage (0-100) required for pairwise SNP distances to be included in analysis.
                        Low-coverage alignments can yield unreliable SNP count.
                        SNP distance below this threshold are excluded. Default: 80.
  --core_threshold CORE_THRESHOLD
                        Panaroo Core-genome sample threshold. The frequency of a gene in your sample required to classify it as 'core'.
                        Range is 0.0 to 1.0. Default is 0.95.
  --core_bootstrap_method {ultrafast,nonparametric}
                        The bootstrap method for core genome phylogeny used for hierarchical clustering.
                        Available options: [ultrafast, nonparametric]
                        ultrafast: Use ultrafast bootstrap method.
                        nonparametric: Use standard nonparametric bootstrap method.
                        Default is ultrafast
  --core_bootstrap_number CORE_BOOTSTRAP_NUMBER
                        The number of bootstrap number for core genome phylogeny.
                        If method is ultrafast, default is 1000.
                        If method is nonparametric, default is 100.
  --group_bootstrap_method {ultrafast,nonparametric}
                        The bootstrap method for phylogeny of each hierarchical group.
                        Available options: [ultrafast, nonparametric]
                        ultrafast: Use ultrafast bootstrap method.
                        nonparametric: Use standard nonparametric bootstrap method.
                        Default is ultrafast
  --group_bootstrap_number GROUP_BOOTSTRAP_NUMBER
                        The number of bootstrap replicates for phylogeny of each hierarchical group.
                        If method is ultrafast, default is 1000.
                        If method is nonparametric, default is 100.
  --use_cladebreaker USE_CLADEBREAKER
                        Use CladeBreaker to restrain the strain composition.
                        Options are [True, False].
                        Default is True.
  --threshold_ceiling THRESHOLD_CEILING
                        The ceiling of the range of SNP distances to search for the optimal phylothreshold.
                        Default is 500.
  --singleton_threshold SINGLETON_THRESHOLD
                        The SNP distance threshold above which every genome in the group is considered a singleton.
                        Default is 100.
  --correction_bootstrap CORRECTION_BOOTSTRAP
                        Minimum bootstrap support threshold for applying phylogenetic corrections to strain composition (default: 0).
  --endpoint ENDPOINT   The endpoint method to use for determing clusters and making plots.
                        Available Options: [plateau, peak, discrepancy, global]
                        plateau : Phylothreshold set at a plateau where further increases no longer change the number or composition of strains within the group
                        peak: Phylothreshold set at the peak number of clones defined within the group.
                        discrepancy: Phylothreshold set at the point where the discrepancy is minimized within the group.
                        global: Phylothreshold set at the first time a global genome is included in any strain within the group.
                        Default is plateau.
  --plateau_length PLATEAU_LENGTH
                        The plateau length for the plateau endpoint method.
                        Only used when endpoint method is 'plateau'. Default is 15
  -t THREADS, --threads THREADS
                        Thread number. Default is 1.
  --prefix PREFIX       Prefix for config file. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS
  --conda_prefix CONDA_PREFIX
                        Directory for conda environments needed for this analysis. If not provided, defaults to <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>
  --force               Bypass system compatibility checks (operating system and available RAM) and force execution of the pipeline.
                        This may cause instability or failures.
```
## Required Input
1. **Original Metadata File(--original_metadata):**

    - Path to the input metadata file containing original genomes.
    - Example metadata file for full mode: 
      
      [Example Input Metadata File](example/example_metadata.txt)
2. **New Metadata File(--new_metadata):**

    - Path to the input metadata file containing new genomes.
    - Example metadata file for full mode: 
      
      [Example Input Metadata File](example/example_metadata.txt)

    - How do I know the Genbank accession for my genomes? See the [Genbank Accession](genbank_accession.md) page.
   - Ensure the new genome assemblies have undergone quality control prior to use as input.
   - All new genome assemblies must belong to the same species; otherwise, Panaroo will fail to generate a core genome alignment due to insufficient shared core genes, and the pipeline will terminate at this step.
3. **Existing THRESHER Strain Identifier Directory(--thresher_output):**

    - Path to the existing THRESHER strain_identifier directory.
    - The existing analysis directory should contain the previous transmission cluster results.
    - (analysis_mode should be 'full' in the prior run)

4. **Species(--species):**
    - Indicate the species being analyzed. Currently supported species include:
        - `sau` (*Staphylococcus aureus*)
        - `sepi` (*Staphylococcus epidermidis*)
        - `cdiff` (*Clostridium difficile*)
        - `kp` (*Klebsiella pneumoniae*)

## Optional Input
- **WhatsGNU Database Path(--whatsgnu_db_path):**
    - The path to the existing WhatsGNU database.
    - If provided, the WhatsGNU database will not be downloaded.
    - If not provided, defaults to `<OUTPUT>/whatsgnu/db`
- **Bakta Database Type(--bakta_db_type):**
    - The type of Bakta database.
    - Available options: `full`, `light`
    - Default is `full`
- **Bakta Database Path(--bakta_db_path):**
    - The path of the directory where the existing Bakta database locates.
    - If provided, the Bakta database will not be downloaded.
    - If not provided, defaults to `<OUTPUT>/bakta/db`
- **SNP Coverage Threshold(--snp_coverage_threshold):**
    - Minimum alignment coverage percentage (0-100) required for pairwise SNP distances to be included in analysis (default: 80). Genome pairs with alignment coverage below this threshold are excluded from downstream cladebreaker analysis.
    - Low-coverage alignments can produce artificially low SNP counts, as unaligned regions are not compared and potential variants in those regions go undetected. This can lead to falsely inflated genomic similarity between genomically distantly related genomes. The default threshold of 80% balances sensitivity with reliability.
- **Core Threshold(--core_threshold):**
    - Panaroo Core-genome sample threshold. The frequency of a gene in your sample required to classify it as 'core'.
    - Range is 0.0 to 1.0. Default is 0.95. 
- **Core Bootstrap Method(--core_bootstrap_method):**
    - The bootstrap method for core genome phylogeny used for hierarchical clustering.
    - Available options: `ultrafast`, `nonparametric`
    - `ultrafast`: Use ultrafast bootstrap method.
    - `nonparametric`: Use standard nonparametric bootstrap method.
    - Default is `ultrafast`
- **Core Bootstrap Number(--core_bootstrap_number):**
    - The number of bootstrap number for core genome phylogeny.
    - If method is `ultrafast`, default is 1000.
    - If method is `nonparametric`, default is 100.
- **Group Bootstrap Method(--group_bootstrap_method):**
    - The bootstrap method for phylogeny of each hierarchical group.
    - Available options: `ultrafast`, `nonparametric`
    - `ultrafast`: Use ultrafast bootstrap method.
    - `nonparametric`: Use standard nonparametric bootstrap method.
    - Default is `ultrafast`
- **Group Bootstrap Number(--group_bootstrap_number):**
    - The number of bootstrap replicates for phylogeny of each hierarchical group.
    - If method is `ultrafast`, default is 1000.
    - If method is `nonparametric`, default is 100.
- **CladeBreaker(--use_cladebreaker):**
    - Whether or not to use CladeBreaker to restrain the strain composition using the closely related genomes in the WhatsGNU database.
    - Enable when investigating putative novel or locally-restricted strains that should be genomically distinct from globally circulating strains. 
    - `True` or `False`, default is `True`.
- **Threshold Ceiling(--threshold_ceiling):**
    - The ceiling of the range tested to search for the optimal phylothreshold (default: 500).
    - This parameter sets the upper limit of SNP phylothresholds considered when determining the optimal phylothreshold for defining strains within hierarchical clustering groups.
    - Adjust this ceiling based on the expected genomic diversity within your dataset and the species being analyzed. A higher ceiling allows for the inclusion of more distantly related genomes in the strain identification process, which may be appropriate for highly diverse species.
- **Singleton Threshold(--singleton_threshold):**
    - The SNP distance threshold above which every genome in the hierarchical clustering group is considered a singleton (default: 100).
    - If the minimal SNP distance between any two genomes in a hierarchical clustering group is equal to or greater than this threshold, all genomes in that group will be classified as singletons, meaning each genome forms its own unique strain.
    - This parameter helps to avoid defining strains in groups where genomes are too distantly related to form meaningful strain. Adjust this threshold based on the expected genomic diversity within your dataset and the species being analyzed.
-  **Correction Bootstrap Threshold(--correction_bootstrap):**
    - Minimum bootstrap support required to apply phylogenetic corrections to SNP strain composition (default: 0). Nodes with bootstrap values below this threshold are excluded from correction, except for the root node which is always retained.
    - Higher thresholds enforce stricter corrections but may reduce the number of corrections applied, particularly in trees with lower overall bootstrap support. The default of 0 applies all possible corrections regardless of bootstrap support, which may be appropriate for small genome groups where high bootstrap values are difficult to achieve. Adjust based on your bootstrap method and resampling depth.
- **Endpoint Method(--endpoint):**
    - The endpoint method to use for determing clusters and making plots.
    - Available Options: `plateau`, `peak`, `discrepancy`, `global`
        - `plateau` : Phylothreshold set at a plateau where further increases no longer change the number or composition of strains within the group
        - `peak`: Phylothreshold set at the peak number of clones defined within the group.
        - `discrepancy`: Phylothreshold set at the point where the discrepancy is minimized within the group.
        - `global`: Phylothreshold set at the first time a global genome is included in any strain within the group.
    - Default is `plateau`.
- **Plateau Length(--plateau_length):**
    - The plateau length for the plateau endpoint method.
    - Only used when endpoint method is `plateau`. Default is 15
- **Threads(--threads / -t):**
    - Thread number. Default is 1. It is highly recommended to increase the thread count using this option to improve performance, regardless of dataset size.
    - Bakta genome annotation runs `{threads}` parallel jobs, each requiring approximately 10 GB of RAM. Ensure your system has sufficient memory to support the requested thread count.
- **Prefix(--prefix):**
    - Prefix for config file. If not provided, defaults to timestamp: `YYYY_MM_DD_HHMMSS`
- **Conda Prefix(--conda_prefix):** 
    - Directory for conda environments needed for this analysis. If not provided, defaults to `<OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>`
- **Force(--force):**
    - Bypass system compatibility checks (operating system and available RAM) and force execution of the pipeline.
    - This may cause instability or failures.
## Output
The output of the New Full mode is in the same format of that of the Full Pipeline mode, in order to fully update all results including phylothresholds, strains, and clusters and be compatible with subsequent new full analyses added to the existing results.
The only new file added is updated cluster summaries that reflect the changes in strain compositions after adding new genomes and recalculating phylothresholds.
1. **Updated Clusters:**

    Summary table of strain/cluster composition with original and newly added genomes: `thresher/output/genomes_summary_new_full.csv` (columns identical to New SNPs mode output).