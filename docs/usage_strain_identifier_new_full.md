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
                        Thread number. Default is the maximum available.
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
    - Thread number. Default is the maximum available.
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