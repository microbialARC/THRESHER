# Strain Identifier - Full Pipeline: Run the complete analysis pipeline from scratch.
```
thresher strain_identifier full-pipeline -h

options:
  -h, --help            show this help message and exit
  --metadata METADATA   Path to the input metadata file.
                        The input metadata file should be a tab-delimited file with 3 or 5 columns.
                        3 columns: genome_name, genome_accession, genome_path
                        5 columns: genome_name, genome_accession, genome_path, patient_id, collection_date
                        At least 4 genomes should be provided to perform the analysis.
  -o OUTPUT, --output OUTPUT
                        Path to output directory. If not provided, defaults to thresher_strain_identifier_output_<YYYY_MM_DD_HHMMSS> under the current working directory.
  --species {sau,sepi,cdiff,kp}
                        Bacteria species.
                        Available options: [sau, sepi, cdiff, kp]
                        sau: Staphylococcus aureus
                        sepi: Staphylococcus epidermidis
                        cdiff: Clostridium difficile
                        kp: Klebsiella pneumoniae
  --analysis_mode ANALYSIS_MODE
                        Whether to make cluster plots and persistence plot.
                        Available Options: [full, lite]
                        full: Determine strains, clusters and make plots.
                        lite: Determine strains without determing clusters and making plots.
                        Default is full
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
                        Default is ultrafast.
  --core_bootstrap_number CORE_BOOTSTRAP_NUMBER
                        The number of bootstrap replicates for core genome phylogeny.
                        If method is ultrafast, default is 1000.
                        If method is nonparametric, default is 100.
  --group_bootstrap_method {ultrafast,nonparametric}
                        The bootstrap method for phylogeny of each hierarchical group.
                        Available options: [ultrafast, nonparametric]
                        ultrafast: Use ultrafast bootstrap method.
                        nonparametric: Use standard nonparametric bootstrap method.
                        Default is ultrafast.
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
                        The plateau length for the plateau endpoint method. Default is 15.
                        Only used when endpoint method is 'plateau'.
  -t THREADS, --threads THREADS
                        Thread number. Default is the maximum available.
  --prefix PREFIX       Prefix for config file, output files, and analysis naming. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS
  --conda_prefix CONDA_PREFIX
                        Directory for conda environments needed for this analysis. If not provided, defaults to <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>
  --force               Bypass system compatibility checks (operating system and available RAM) and force execution of the pipeline.
                        This may cause instability or failures.
```
## Required Input
**Note:** At least 4 genomes are required to perform the analysis. 
1. **Genome Assemblies:**
   - Genome assemblies for all samples to be analyzed.
2. **Input Metadata(--metadata):**
   - Path to a tab-delimited text file containing at least three columns (for lite mode) with no header. For full mode, additional columns are required.
   - Columns:
     - **Column 1:** Genome name (required for both lite and full modes).
     - **Column 2:** GenBank accession number (required for both lite and full modes). If unavailable, use "new" (all lowercase).
     - **Column 3:** Path to the genome (required for both lite and full modes).
     - **Column 4:** Patient ID (required for full mode).
     - **Column 5:** Collection date in the format `yyyy-mm-dd` (required for full mode).
    - Example metadata file for full mode: 
      
      [Example Input Metadata File](example/example_metadata.txt)

3. **Species(--species):**

   - Indicate the species being analyzed. Currently supported species include:
    - `sau` (*Staphylococcus aureus*)
    - `sepi` (*Staphylococcus epidermidis*)
    - `cdiff` (*Clostridium difficile*)
    - `kp` (*Klebsiella pneumoniae*)

## Optional Input
1. **Output Directory(--output):**
   - Path to the output directory. If not provided, defaults to `thresher_strain_identifier_output_<YYYY_MM_DD_HHMMSS>` under the current working directory.
2. **Analysis Mode(--analysis_mode):**
   - Choose between `full` (default) and `lite` modes.
   - `full`: Determines strains, clusters, and generates plots.
   - `lite`: Determines strains without determining clusters or generating plots.
3. **WhatsGNU Database Path(--whatsgnu_db_path):**
    - Path to an existing WhatsGNU database. If not provided, the database will be downloaded
      to `<OUTPUT>/whatsgnu/db`.

4. **Bakta Database Type(--bakta_db_type) and Path(--bakta_db_path):**
    - Specify the type of Bakta database (`full` or `light`, default is `full`).
    - Path to an existing Bakta database. If not provided, the database will be downloaded to `<OUTPUT>/bakta/db`.

5. **Bootstrap Methods and Numbers(--core/group_bootstrap_method, --core/group_bootstrap_number):**
    - Specify bootstrap methods (`ultrafast` or `nonparametric`) and the number of replicates for both core genome phylogeny and group phylogeny.
    - Default methods are `ultrafast`, with default replicate numbers of 1000 for ultrafast and 100 for nonparametric.

6. **Endpoint Method(--endpoint):**
    - Choose the endpoint method for determining clusters and generating plots. Options include `plateau` (default), `peak`, `discrepancy`, and `global`.
    - `plateau`: Phylothreshold set at a plateau where further increases no longer change the number or composition of strains within the group.
    - `peak`: Phylothreshold set at the peak number of clones defined within the group.
    - `discrepancy`: Phylothreshold set at the point where the discrepancy is minimized within the group.
    - `global`: Phylothreshold set at the first time a global genome is included in any strain within the group.
7. **Plateau Length(--plateau_length):**
    - Specify the plateau length for the plateau endpoint method (default is 15). Only used when the endpoint method is `plateau`.
8. **Threads(--threads / -t):**
    - Number of threads to use (default is the maximum available).
9. **Prefix(--prefix):**
    - Prefix for config file, output files, and analysis naming. If not provided, defaults to a timestamp in the format `YYYY_MM_DD_HHMMSS`.
10. **Conda Prefix(--conda_prefix):**
    - Directory for conda environments needed for this analysis. If not provided, defaults to `<OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>`.
11. **Force Execution(--force):**
    - Bypass system compatibility checks (operating system and available RAM) and force execution of the pipeline. This may cause instability or failures.
    
## Output
1. **Config Files:**
  - Config file used for the analysis: `config/config_{prefix}.yaml`
2. **MLST Analysis:**
  - Summary:
    - `mlst/summary/mlst_summary.csv`
    - Columns included in the MLST summary csv file:
      - genome: Genome name
      - ST: Sequence Type
      - MLST: Broader MLST grouping when available.
        - e.g., for *Staphylococcus aureus*, the Clonal Complex (CC).

  - Raw MLST results:
    - `mlst/raw/{genome_file_name}_mlst.csv`

3. **BlastX results for MRSA: (species: `sau`):**
  - Summary:
    - `blastx/mrsa/output/summary/blastx_MRSA_strains.csv`
    Columns included in the MRSA blastx strains summary csv files:
      - strain: Strain ID identified by Thresher
      - MRSA: Whether the strain is identified as MRSA (MSSA or MRSA)
    - `blastx/mrsa/output/summary/blastx_MRSA_genomes.csv`
    Columns included in the MRSA blastx genomes summary csv files:
      - genome: Genome name
      - MRSA: Whether the genome is identified as MRSA (MSSA or MRSA)
  - Raw blastx results:
    - `blastx/mrsa/output/raw/{genome_name}_blastx_mrsa.tsv`
4. **Bakta Annotation:**
  - `bakta_annotation/{genome_name}/`

5. **Assembly Scan:**
  - `assembly_scan/{genome_name}_assembly_scan.txt`

6. **Pair-wise comparison using Mummer4:**
  - Study concatenated reports: `mummer4_study/{genome_name}_concatenated.report`
  - Global concatenated reports: `mummer4_global/{genome_name}_concatenated.report`
  - Global SNP matrix: `mummer4_global/global_snp_matrix.RDS`
  - Study SNP matrix: `mummer4_study/study_snp_matrix.RDS`

7. **WhatsGNU:**
  - `whatsgnu/{genome_name}/{genome_name}_WhatsGNU_topgenomes.txt`

8. **Global genomes downloaded using Datasets:**
  - `datasets_topgenomes/`

9. **Pan-genome analysis using Panaroo:**
  - `panaroo/`

10. **Snippy:**
  - Snippy analysis of groups: `snippy/output/`
  - Snippy core-genome alignments: `snippy/output/cleaned_aln/`

11. **IQTree:**
  - Comprehensive tree:
    - `iqtree/comprehensive_tree/`
  - Group trees:
    - `iqtree/group_tree/`

12. **Hierarchical Clustering:**
  - Hierarchical clustering groups(RDS): `thresher/input/hierarchical_clustering_groups.RDS`
  - Hierarchical clustering groups(csv): `thresher/input/hierarchical_clustering_groups_simplified.csv`
  Columns included in the hierarchical clustering groups simplified csv file:
    - genome: Genome name
    - group: Hierarchical clustering group ID
    - overlimit: Whether the minimal SNP distance of the genome to any other genome in the same group exceeds 1000 SNPs,
    If TRUE, the genome will be considered an outlier and excluded from group analysis even though it is assigned to a hierarchical clustering group.

  - Thresher input: `thresher/input/thresher_input.RDS`

13. **Strain Composition:**

The strain composition files generated by all 4 endpoint methods contain the following columns:
1. strain_id : Strain ID assigned by Thresher using the specified endpoint method.
2. genome : Genome name.
  - Plateau strains:
    - `thresher/output/plateau_strains.RDS`
    - `thresher/output/plateau_strains.csv`
  - Peak strains:
    - `thresher/output/peak_strains.RDS`
    - `thresher/output/peak_strains.csv`
  - Discrepancy strains:
    - `thresher/output/discrepancy_strains.RDS`
    - `thresher/output/discrepancy_strains.csv`
  - Global strains:
    - `thresher/output/global_strains.RDS`
    - `thresher/output/global_strains.csv`
  - Group-specific thresholds:
    - Plateau thresholds: `thresher/output/group_plateau.csv`
      - Columns included in the group-specific thresholds csv files:
        - group: Hierarchical clustering group ID
        - plateau: Phylothreshold determined for the group using the plateau endpoint method.
        - plateau_length: Plateau length used for determining the phylothreshold.
    - Peak thresholds: `thresher/output/group_peak.csv`
      - Columns included in the group-specific thresholds csv files:
        - group: Hierarchical clustering group ID
        - peak: Phylothreshold determined for the group using the peak endpoint method.
    - Discrepancy thresholds: `thresher/output/group_discrepancy.csv`
      - Columns included in the group-specific thresholds csv files:
        - group: Hierarchical clustering group ID
        - discrepancy: Phylothreshold determined for the group using the discrepancy endpoint method.
    - Global thresholds: `thresher/output/group_global.csv`
      - Columns included in the group-specific thresholds csv files:
        - group: Hierarchical clustering group ID
        - global: Phylothreshold determined for the group using the global endpoint method.

14. **Sanity Check Plots:**
  - Plateau: `thresher/output/QC/plateau_qc_plot.pdf`
  - Peak: `thresher/output/QC/peak_qc_plot.pdf`
  - Global: `thresher/output/QC/global_qc_plot.pdf`
  - Discrepancy: `thresher/output/QC/discrepancy_qc_plot.pdf`
  - Figure descriptions:
    - X-axis: Phylogenetic Average Distance

      The average phylogenetic distance between two strains, calculated by averaging the pairwise branch-length distances between all genomes in strain A and all genomes in strain B on the core genome phylogeny (comprehensive tree).

    - Y-axis: gSNP Average Distance
    
      The average SNP distance between two strains, calculated by averaging the pairwise SNP distances between all genomes in strain A and all genomes in strain B. Pairwise SNP distances are computed using MUMmer4.
    
    - Colored dots: Each dot represents a comparison between two strains. Blue dots indicate strain pairs within the same hierarchical clustering group. Red dots indicate strain pairs from different hierarchical clustering groups. 

    - Doted line: A vertical dotted line at gSNP distance of 100.

15. **Multi-SNP-Threshold Plot:**
  - `thresher/output/MSTP/{Group_ID}_MSTP.pdf`
  - Figure descriptions:
    - Shared X-axis: SNP Threshold tested

      The SNP thresholds tested to define strains.

    - Upper panel: The median and mean bootstrap values of the clades containing the strains defined at the SNP threshold.

    - Lower panel: 
      - Lines: The number of singletons, clones, and the discrepant genomes when the strain compositions determined by pairwise SNP distances mapped to the corresponding reference-based phylogenetic tree of the group.
      - Annotations: The final phylothresholds determined by each endpoint method (plateau, peak, discrepancy, global) are indicated on the plot.

16. **Sanity Check Tables:**

  All csv files contain the following columns:

  1. subject: subject strain ID
  2. query: query strain ID
  3. snp_average_distance: average SNP distance between the subject and query strains
  4. phylogeny_average_distance: average phylogenetic distance between the subject and query strains
  5. same_group: whether the subject and query strains belong to the same hierarchical clustering group

  - Plateau: `thresher/output/QC/plateau_qc_table.csv`
  - Peak: `thresher/output/QC/peak_qc_table.csv`
  - Global: `thresher/output/QC/global_qc_table.csv`
  - Discrepancy: `thresher/output/QC/discrepancy_qc_table.csv`

17. **Visualization:**
  - Core-gene comprehensive tree annotated with hierarchical clustering groups: 
    - `plots/comprehensive_tree_group.pdf`
    - `plots/comprehensive_tree_group.RDS`(R object)
  - Core-gene comprehensive tree annotated with MLST:
    - `plots/comprehensive_tree_mlst.pdf`
    - `plots/comprehensive_tree_mlst.RDS`(R object)
  - SNP distance: `plots/SNP_Distance.pdf`

18. **Cluster Plots:** (if full analysis mode is enabled)
  - Cluster plots: `plots/ClusterPlots/`
  - Persistence plot: `plots/PersistencePlot.pdf`
  - Clusters summary:
    - RDS: `thresher/output/clusters_summary.RDS`
    - CSV: `thresher/output/clusters_summary.csv`
    Columns included in the clusters summary csv file:
      - cluster: Cluster ID assigned by Thresher.
      - strain: Strain ID assigned by Thresher.
      - MLST: Broader MLST grouping when available.
        - e.g., for *Staphylococcus aureus*, the Clonal Complex (CC).
      - AMR: If species is `sau`, whether the strain is identified as MRSA (MSSA or MRSA).
      - genomes: genome names belonging to the cluster, separated by "|".
      - patients: patient IDs involved in the cluster, separated by "|".
      - first_seen: The earliest collection date among genomes in the cluster.
      - last_seen: The latest collection date among genomes in the cluster.
      - persistence: The persistence of the cluster in days (last_seen - first_seen).
      