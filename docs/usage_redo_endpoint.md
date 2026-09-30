# THRESHER Redo Endpoint: Only rerun the final endpoint analysis using existing intermediate files.
```
thresher redo-endpoint -h

options:
  -h, --help            show this help message and exit
   --original_metadata ORIGINAL_METADATA
                        Path to the original input file used for the THRESHER full run.
                        The file must be tab-delimited and contain 5 columns because redo-endpoint requires patient ID and collection date information.
  --thresher_output THRESHER_OUTPUT
                        Path to the existing THRESHER directory.
                        The existing analysis directory should contain the previous analysis results.
  --output OUTPUT       Path to output directory.
                        If not provided, defaults to thresher_strain_identifier_redo_endpoint_output_<YYYY_MM_DD_HHMMSS> under the current working directory.
  --endpoint ENDPOINT   The endpoint method to use for determing clusters and making plots.
                        Available Options: [plateau, peak, discrepancy, public]
                        plateau : Phylothreshold set at a plateau where further increases no longer change the number or composition of strains within the group
                        peak: Phylothreshold set at the peak number of clones defined within the group.
                        discrepancy: Phylothreshold set at the point where the discrepancy is minimized within the group.
                        public: Phylothreshold set at the first time a public genome is included in any strain within the group.
                        Default is plateau.
  --prefix PREFIX       Prefix for config files. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS
  --conda_prefix CONDA_PREFIX
                        Directory for conda environments needed for this analysis. If not provided, defaults to <THRESHER_OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>
```
## Required Input
1. **Original Metadata File(--original_metadata):**

   Path to the original input metadata file used for the THRESHER full run.
2. **Existing THRESHER Output Directory(--thresher_output):**
   
   Path to the existing THRESHER output directory containing the previous analysis results.
3. **Endpoint Method(--endpoint):**
   
   The endpoint method to use for determining clusters and making plots. Available options are:
   - plateau
   - peak
   - discrepancy
   - public

## Optional Input
1. **Output Directory(--output):**
   
   Path to the output directory. If not provided, defaults to `thresher_redo_endpoint_output_<YYYY_MM_DD_HHMMSS>` under the current working directory.

2. **Prefix(--prefix):**

    Prefix for config files. If not provided, defaults to timestamp: `YYYY_MM_DD_HHMMSS`.

3. **Conda Environment Directory(--conda_prefix):**

    Directory for conda environments needed for this analysis. If not provided, defaults to `<THRESHER_OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>`. You can reuse the conda environments from previous THRESHER runs to save time and disk space.

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