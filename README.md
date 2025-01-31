# **THRESHER**  
THRESHER is a pipeline designed to implement dynamic and unbiased SNP thresholds, corrected for phylogenetic structures. With the SNP thresholds, this tool determines clonality and identifies transmission clusters.

---

## **Workflow**  
The THRESHER workflow is illustrated below:  
![Thresher workflow](data/workflow_schematic.png)

---

## **Installation**  
THRESHER utilizes a Snakemake pipeline complemented by a Python script for input validation, configuration file generation, and execution of the workflow.

### **Bioconda Installation**  
*Coming soon.*

### **Manual Installation via Git**  
1. Clone the repository:  
   ```bash
   git clone https://github.com/microbialARC/THRESHER
   cd THRESHER
2. Install dependencies:
    ```bash
    conda env create -f thresher.yml
3. Activate THRESHER environment and install:
    ```bash
    conda activate thresher
    bash install.sh
## **Usage**  
### System Requirements
- **Minimum RAM:** 64 Gb (required)
- **Recommended CPUs:** 16 or more cores
### Command line options
```
usage: thresher [-h] -i INPUT -o OUTPUT --species SPECIES [--whatsgnu_db_path WHATSGNU_DB_PATH] [-db BAKTA_DB] [--bakta_db_path BAKTA_DB_PATH] [-t THREADS] [--memory MEMORY] [-m MODE] [-e ENDPOINT]
                [--plateau_length PLATEAU_LENGTH]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input file, which is a CSV file with at least 3 columns(lite mode), separated by tab, no header.
                        The first column is the genome name(lite & full).
                        The second column is the GenBank accession number(lite & full). If the accession number is not available, put "new"(all lowercase) in the column.
                        The third column is the path to the genome(lite & full).
                        The forth column is the patient ID (full).
                        The fifth column is the collection date. Format: yyyy-mm-dd (full).
  -o OUTPUT, --output OUTPUT
                        Path to the output directory.
  --species SPECIES     Bacteria species.
                        Available options: [sau, sepi, cdiff, kp]
                        sau: Staphylococcus aureus
                        sepi: Staphylococcus epidermidis
                        cdiff: Clostridium difficile
                        kp: Klebsiella pneumoniae
  --whatsgnu_db_path WHATSGNU_DB_PATH
                        The path to the existing WhatsGNU database.
                        If provided, the WhatsGNU database will not be downloaded
  -db BAKTA_DB, --bakta_db BAKTA_DB
                        Bakta database.
                        Available options: [full, light]
                        Default is full
  --bakta_db_path BAKTA_DB_PATH
                        The path of the directory where the existing Bakta database locates.
                        If provided, the Bakta database will not be downloaded
  -t THREADS, --threads THREADS
                        Number of threads to use.
                        Default is 16
  --memory MEMORY       RAM(Gb) to use.
                        WhatsGNU database requires loading the database into memory.
                        All WhatsGNU databases require at least 64Gb to run.
                        Default is 64Gb
  -m MODE, --mode MODE  Whether to make cluster plots and persistence plot.
                        Availabel Options: [full, lite] 
                        full: Determine strains, clusters and make plots.
                        lite: Determine strains without determing clusters and making plots.
                        Default is lite
  -e ENDPOINT, --endpoint ENDPOINT
                        The endpoint method to use for determing clusters and making plots.
                        Availabel Options: [plateau, peak, discrepancy, global]
                        plateau : Threshold set at the plateau where strain compositions are stable in the group.
                        peak: Threshold set at the peak number of clones defined in the group.
                        discrepancy: Threshold set at the point where the discrepancy is minimized in the group. 
                        global: Threshold set at the first time a global genome is included in any strain in the group.
                        Default is plateau
  --plateau_length PLATEAU_LENGTH
                        The plateau length for the plateau endpoint method.
                        Default is 15
```
### Required Inputs

1. **Input Metadata**

   - Path to a tab-delimited CSV file containing at least three columns (for lite mode) with no header. For full mode, additional columns are required.
   - Columns:
     - **Column 1:** Genome name (required for both lite and full modes).
     - **Column 2:** GenBank accession number (required for both lite and full modes). If unavailable, use "new" (all lowercase).
     - **Column 3:** Path to the genome (required for both lite and full modes).
     - **Column 4:** Patient ID (required for full mode).
     - **Column 5:** Collection date in the format `yyyy-mm-dd` (required for full mode).

2. **Output Directory**

   - Specify the directory where output files will be saved.

3. **Species**

   - Indicate the species being analyzed. Currently supported species include:
     - `sau` (Staphylococcus aureus)
     - `sepi` (Staphylococcus epidermidis)
     - `cdiff` (Clostridium difficile)
     - `kp` (Klebsiella pneumoniae)

### Outputs

The pipeline generates the following outputs, organized into corresponding directories based on the input configurations:

- **MLST Analysis:**
  - `mlst/mlst_results.csv`

- **BlastX results for MRSA: (species: `sau`):**
    - `blastx/mrsa/output/blastx_MRSA_strains.csv`
    - `blastx/mrsa/output/blastx_MRSA_genomes.csv`

- **Bakta Annotation:**
  - `bakta_annotation/`

- **Assembly Scan:**
  - `assembly_scan/{genome_name}_assembly_scan.txt`

- **Pair-wise comparison using Mummer4:**
  - Study concatenated reports: `mummer4_study/{genome_name}_concatenated.report`
  - Global concatenated reports: `mummer4_global/{genome_name}_concatenated.report`
  - Global SNP matrix: `mummer4_global/global_snp_matrix.RDS`
  - Study SNP matrix: `mummer4_study/study_snp_matrix.RDS`

- **WhatsGNU:**
  - `whatsgnu/{genome_name}/{genome_name}_WhatsGNU_topgenomes.txt`

- **Global genomes downloaded using Datasets:**
  - `datasets_topgenomes/`

- **Pan-genome analysis using Panaroo:**
  - `panaroo/`

- **Snippy:**
  - Snippy analysis of groups: `snippy/output/`
  - Snippy core-genome alignments: `snippy/output/cleaned_aln/`

- **IQTree:**
  - Everything tree:
    - `iqtree/everything_tree/`
  - Group trees:
    - `iqtree/group_tree`

#### Thresher Outputs
- **Input Data:**
  - Hierarchical clustering groups(RDS): `thresher/input/hierarchical_clustering_groups.RDS`
  - Hierarchical clustering groups(csv): `thresher/input/hierarchical_clustering_groups_simplified.csv`
  - Thresher input: `thresher/input/thresher_input.RDS`

- **Output Data:**
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
    - Peak thresholds: `thresher/output/group_peak.csv`
    - Discrepancy thresholds: `thresher/output/group_discrepancy.csv`
    - Global thresholds: `thresher/output/group_global.csv`

- **Sanity Check Plots:**
  - Plateau: `thresher/output/QC/plateau_qc_plot.pdf`
  - Peak: `thresher/output/QC/peak_qc_plot.pdf`
  - Global: `thresher/output/QC/global_qc_plot.pdf`
  - Discrepancy: `thresher/output/QC/discrepancy_qc_plot.pdf`

- **Sanity Check Tables:**
  - Plateau: `thresher/output/QC/plateau_qc_table.csv`
  - Peak: `thresher/output/QC/peak_qc_table.csv`
  - Global: `thresher/output/QC/global_qc_table.csv`
  - Discrepancy: `thresher/output/QC/discrepancy_qc_table.csv`

#### Plots
- **Visualization Outputs:**
  - Core-gene tree: `plots/everything_tree.pdf`
  - SNP distance: `plots/SNP_Distance.pdf`

#### Full Mode Outputs (if enabled)
- **Cluster Plots:**
  - Cluster plots: `plots/ClusterPlots/`
  - Persistence plot: `plots/PersistencePlot.pdf`
  - Clusters summary:
    - RDS: `thresher/output/clusters_summary.RDS`
    - CSV: `thresher/output/clusters_summary.csv`

## **Author**
Qianxuan(Sean) She


[![PennMedicine](data/PennMedicine.png)](https://www.pennmedicine.org/)  
[![CHOP_Research](data/CHOP_Research.png)](https://www.research.chop.edu/)  
[![PennCHOP](data/PennCHOP.png)](https://www.research.chop.edu/pennchop-microbiome-program)
