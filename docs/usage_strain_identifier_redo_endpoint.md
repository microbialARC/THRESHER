# Strain Identifier - Redo Endpoint: Only rerun the final endpoint analysis using existing intermediate files.
```
thresher strain_identifier redo-endpoint -h

options:
  -h, --help            show this help message and exit
  --original_metadata ORIGINAL_METADATA
                        Path to the original input file used for the THRESHER full-pipeline.
                        The file must be tab-delimited and contain 5 columns becuase the full-pipeline mode requires patient ID and collection date information.
                        5 columns: genome_name, genome_accession, genome_path, patient_id, collection_date.
  --thresher_output THRESHER_OUTPUT
                        Path to the existing THRESHER directory.
                        The existing analysis directory should contain the previous analysis results.
  --output OUTPUT       Path to output directory.
                        If not provided, defaults to thresher_strain_identifier_redo_endpoint_<YYYY_MM_DD_HHMMSS> under the current working directory.
  --endpoint ENDPOINT   The endpoint method to use for determing clusters and making plots.
                        Available Options: [plateau, peak, discrepancy, global]
                        plateau : Phylothreshold set at a plateau where further increases no longer change the number or composition of strains within the group
                        peak: Phylothreshold set at the peak number of clones defined within the group.
                        discrepancy: Phylothreshold set at the point where the discrepancy is minimized within the group.
                        global: Phylothreshold set at the first time a global genome is included in any strain within the group.
                        Default is plateau.
  --prefix PREFIX       Prefix for config files. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS
  --conda_prefix CONDA_PREFIX
                        Directory for conda environments needed for this analysis. If not provided, defaults to OUTPUT/conda_envs
```
## Required Input
1. **Original Metadata File:**

   Path to the original input file used for the THRESHER full-pipeline.
2. **Existing THRESHER Output Directory:**
   
   Path to the existing THRESHER directory containing previous analysis results.

3. **Endpoint Method:**
   
   The endpoint method to use for determining clusters and making plots. Available options are:
   - plateau
   - peak
   - discrepancy
   - global

## Optional Input
1. **Output Directory:**
   
   Path to the output directory. If not provided, defaults to `thresher_strain_identifier_redo_endpoint_<YYYY_MM_DD_HHMMSS>` under the current working directory.

2. **Prefix:**

    Prefix for config files. If not provided, defaults to timestamp: `YYYY_MM_DD_HHMMSS`.

3. **Conda Environment Directory:**
    Directory for conda environments needed for this analysis. If not provided, defaults to `OUTPUT/conda_envs`. You can reuse the conda environments from previous THRESHER runs to save time and disk space.

## Output
1. **Config Files:**
  - Config file used for the analysis: `config/config_{prefix}.yaml`
2. **Updated Clusters and Plots:**
  - Updated clusters: 
    - `clusters_summary_redo_endpoint.csv`
    - `clusters_details_redo_endpoint.RDS`(R object)

  - Updated plots:
    - Cluster plots: `plots/Cluster{Cluster ID}.pdf` (for each cluster)
    - Persistence plot PDF: `plots/PersistencePlot_redo_endpoint.pdf`
    - Persistence plot R object: `plots/PersistencePlot_redo_endpoint.RDS`(R object)