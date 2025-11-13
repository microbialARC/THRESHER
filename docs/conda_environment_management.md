# Conda Environment Management
THRESHER uses Conda and Snakemake to automatically manage dependencies and create required environments. All four Strain Identifier modes share the same Conda environments. After the initial installation, subsequent runs can reuse existing environments by specifying the `--conda_prefix` flag in the THRESHER command line, eliminating the need to recreate environments.

For example, run the Strain Identifier in Full Pipeline mode with previously created conda environments:
```
thresher strain_identifier full-pipeline --metadata /path/to/metadata.txt --species <species> --conda_prefix /path/to/conda/environments/
```
