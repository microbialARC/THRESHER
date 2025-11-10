# **THRESHER**  
THRESHER is a bacterial genomics toolkit with three core functionalities:
- Strain Identifier: 
Determines bacterial clonality and identifies strains/transmission clusters using phylothresholds (phylogenetically-corrected SNP thresholds). 

- Genome Profiler (*Coming soon*): 
Infer substitution probabilities and mobile genetic element (MGE) dynamics leveraging publicly available complete genomes.

- Evolution Simulator (*Coming soon*):
Models bacterial evolution through substitution, gene gain/loss, and recombination events.

## **Development Status**
Please be advised that this pipeline is still in its early development stage. It is subject to significant changes in terms of options, outputs, and other functionalities. Users should be prepared for potential modifications and updates in future releases.

This page is temporary and will be replaced with a proper documentation site in the future. 

## **Workflow**  
### Strain Identifier
- Full Pipeline
  ![Strain Identifier Full-Pipeline Workflow](data/workflow/Thresher_strain_identifier_full_pipeline_workflow.png)
- Redo Endpoint
  ![Strain Identifier Redo Endpoint Workflow](data/workflow/Thresher_strain_identifier_redo_endpoint_workflow.png)
  <span style="color:#91a01e;"><strong>Elements highlighted in green (#91a01e) indicate updates applied to existing results.</strong></span>
- New SNPs 
  ![Strain Identifier New SNPs Workflow](data/workflow/Thresher_strain_identifier_new_snps_workflow.png)
  <span style="color:#91a01e;"><strong>Elements highlighted in green (#91a01e) indicate updates applied to existing results.</strong></span>
- New Full 
  ![Strain Identifier New Full Workflow](data/workflow/Thresher_strain_identifier_new_full_workflow.png)
  <span style="color:#91a01e;"><strong>Elements highlighted in green (#91a01e) indicate updates applied to existing results.</strong></span>

### Genome Profiler (*Coming soon.*)

### Evolution Simulator (*Coming soon.*)
- Preset 
- Custom 

## **Installation**
THRESHER utilizes a Snakemake pipeline complemented by a Python script for input validation, configuration file generation, and execution of the workflow.

### **Bioconda Installation**  
*Coming soon.*

### **Manual Installation via Git**  
1. Clone the GitHub repository:  
   ```bash
   git clone https://github.com/microbialARC/THRESHER
   cd THRESHER
   ```

2. Install THRESHER using the installation script:
    ```bash
    bash install.sh
    ```
3. Activate Thresher environment and display help message:
    ```bash
    conda activate thresher
    thresher -h
    ```

## **Usage**  
### Strain Identifier
Strain Identifier comprises four modes:
#### **Full Pipeline:** Run the complete analysis pipeline from scratch.
```
thresher strain_identifier full_pipeline -h
```

#### **Redo Endpoint:** Only rerun the final endpoint analysis using existing intermediate files.
```
thresher strain_identifier redo_endpoint -h
```

  - Current supported endpoints. For each hierarchical clustering group defined in the core gene comphrehensive tree:
    - Plateau: Phylothreshold set at a plateau where further increases no longer change the number or composition of strains within the group.
    - Peak: Phylothreshold set at the peak number of clones defined within the group.
    - Discrepancy:  Phylothreshold set at the point where the discrepancy is minimized within the group.
      - Global:  Phylothreshold set at the first time a global genome is included in any strain within the group.

#### **New SNPs:** Update existing strain/transmission compositions with new genomes using predefined phylothresholds.
```
thresher strain_identifier new_snps -h
```

#### **New Full:** Rerun the full pipeline to update the phylothresholds, and update strain/transmission compositions with new genomes.
```
thresher strain_identifier new_full -h
```

## Requirements
- **Operating System:** Linux
- **Minimum RAM:** 60 GB — required for Full Pipeline and New Full modes
(required for WhatsGNU database. This will be reduced in future versions after we optimize the WhatsGNU database usage.)
- Internet connection is required.

## Conda Environment Management
THRESHER uses Conda and Snakemake to automatically manage dependencies and create required environments. All four Strain Identifier modes share the same Conda environments. After the initial installation, subsequent runs can reuse existing environments by specifying the `--conda_prefix` flag in the THRESHER command line, eliminating the need to recreate environments.

## Resuming Interrupted Runs
Given the pipeline's complexity and potentially long runtime, THRESHER utilizes Snakemake's `--rerun-incomplete` functionality to handle interruptions. If the Strain Identifier pipeline is interrupted, Snakemake can resume from the last successful step without restarting the entire analysis.

To resume:
1. First run the unlock script: `bash <output_directory>/resume/unlock_snakemake_<output_prefix>.sh`
2. Then run the resume script: `bash <output_directory>/resume/resume_snakemake_<output_prefix>.sh`

This approach preserves completed work and saves computational resources.

### Required Inputs for Full Pipeline
**Note:** At least 4 genomes are required to perform the analysis. 
1. **Genome Assemblies:**
   - Genome assemblies for all samples to be analyzed.
2. **Input Metadata**
   - Path to a tab-delimited text file containing at least three columns (for lite mode) with no header. For full mode, additional columns are required.
   - Columns:
     - **Column 1:** Genome name (required for both lite and full modes).
     - **Column 2:** GenBank accession number (required for both lite and full modes). If unavailable, use "new" (all lowercase).
     - **Column 3:** Path to the genome (required for both lite and full modes).
     - **Column 4:** Patient ID (required for full mode).
     - **Column 5:** Collection date in the format `yyyy-mm-dd` (required for full mode).

3. **Output Directory**

   - Specify the directory where output files will be saved.

4. **Species**

   - Indicate the species being analyzed. Currently supported species include:
    - `sau` (*Staphylococcus aureus*)
    - `sepi` (*Staphylococcus epidermidis*)
    - `cdiff` (*Clostridium difficile*)
    - `kp` (*Klebsiella pneumoniae*)

### Outputs of Strain Identifier (Full Pipeline mode)

- **MLST Analysis:**
  - Summary:
    - `mlst/summary/mlst_summary.csv`
  - Raw MLST results:
    - `mlst/raw/{genome_file_name}_mlst.csv`

- **BlastX results for MRSA: (species: `sau`):**
  - Summary:
    - `blastx/mrsa/output/summary/blastx_MRSA_strains.csv`
    - `blastx/mrsa/output/summary/blastx_MRSA_genomes.csv`
  - Raw blastx results:
    - `blastx/mrsa/output/raw/{genome_name}_blastx_mrsa.tsv`
- **Bakta Annotation:**
  - `bakta_annotation/{genome_name}/`

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
  - Comprehensive tree:
    - `iqtree/comprehensive_tree/`
  - Group trees:
    - `iqtree/group_tree/`

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
  - Core-gene comprehensive tree annotated with hierarchical clustering groups: 
    - `plots/comprehensive_tree_group.pdf`
    - `plots/comprehensive_tree_group.RDS`(R object)
  - Core-gene comprehensive tree annotated with MLST:
    - `plots/comprehensive_tree_mlst.pdf`
    - `plots/comprehensive_tree_mlst.RDS`(R object)
  - SNP distance: `plots/SNP_Distance.pdf`

#### Full Mode Outputs (if enabled)
- **Cluster Plots:**
  - Cluster plots: `plots/ClusterPlots/`
  - Persistence plot: `plots/PersistencePlot.pdf`
  - Clusters summary:
    - RDS: `thresher/output/clusters_summary.RDS`
    - CSV: `thresher/output/clusters_summary.csv`

### Main Outputs of Strain Identifier (Redo Endpoint mode)
After rerunning the endpoint analysis, updated cluster summaries and plots will be generated based on the predefined phylothresholds.
#### Updated Clusters:
- RDS: `files/clusters_summary_redo_endpoint.RDS`
- CSV: `files/clusters_summary_redo_endpoint.csv`
#### Updated Plots:
- Cluster plots: `plots/Cluster{Cluster ID}.pdf` (for each cluster)
- Persistence plot PDF: `plots/PersistencePlot_redo_endpoint.pdf`
- Persistence plot R object: `plots/PersistencePlot_redo_endpoint.RDS`
### Main Outputs of Strain Identifier (New SNPs mode)
#### Updated Genomes:
- Plateau: `thresher/output/new_plateau_genomes.csv`
- Peak: `thresher/output/new_peak_genomes.csv`
- Discrepancy: `thresher/output/new_discrepancy_genomes.csv`
- Global: `thresher/output/new_global_genomes.csv`
#### Updated Strains:
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
#### Updated Clusters:
- RDS: `thresher/output/clusters_summary_new_snps.RDS`
- CSV: `thresher/output/clusters_summary_new_snps.csv`
- Summary table listing each strain/cluster and its composition, including original genomes and newly added genomes: `thresher/output/genomes_summary_new_snps.csv`
  - Column 1 `genome_name`: Name of the genome
  - Column 2 `genome_category`: "original" or "new" indicating whether the genome was part of the original analysis or newly added
  - Column 3 `original_strain`: ID of the strain the genome belongs to in the original analysis ("New-Genome" for new genomes)
  - Column 4 `new_strain`: ID of the strain the genome belongs to in the updated analysis
  - Column 5 `original_cluster`: ID of the cluster the genome belongs to in the original analysis ("New-Cluster" for new genomes, "Non-Cluster" if not belonging to any cluster)
  - Column 6 `new_cluster`: ID of the cluster the genome belongs to in the updated analysis ("Non-Cluster" if not belonging to any cluster)

#### Updated Plots:
- Cluster plots: `plots/Cluster{Cluster ID}.pdf` (for each cluster)
- Persistence plot PDF: `plots/PersistencePlot_new_snps.pdf`
- Persistence plot R object: `plots/PersistencePlot_new_snps.RDS`
### Main Outputs of Strain Identifier (New Full mode)
The output of the New Full mode is in the same format of that of the Full Pipeline mode, in order to fully update all results including phylothresholds, strains, and clusters and be compatible with subsequent new full analyses added to the existing results.
The only new file added is updated cluster summaries that reflect the changes in strain compositions after adding new genomes and recalculating phylothresholds.
#### Updated Clusters:
- Summary table of strain/cluster composition with original and newly added genomes: `thresher/output/genomes_summary_new_full.csv` (columns identical to New SNPs mode output).

## **Author**
Qianxuan (Sean) She


[![PennMedicine](data/logo/PennMedicine.png)](https://www.pennmedicine.org/) [![CHOP_Research](data/logo/CHOP_Research.png)](https://www.research.chop.edu/) [![PennCHOP](data/logo/PennCHOP.png)](https://www.research.chop.edu/pennchop-microbiome-program)
