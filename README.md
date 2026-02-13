# **THRESHER**  
THRESHER is a bacterial genomics toolkit with three core functionalities:
- Strain Identifier: 
Identifies strains/transmission clusters using phylothresholds (phylogenetically-corrected SNP thresholds). 

- Genome Profiler:
Infers substitution probabilities and mobile genetic element (MGE) dynamics leveraging publicly available genomes.

- Evolution Simulator:
Models bacterial evolution through substitution, gene gain/loss, and recombination events.

## **Workflow**  
### Strain Identifier
- Full Pipeline
  ![Strain Identifier Full-Pipeline Workflow](data/workflow/Thresher_strain_identifier_full_pipeline_workflow.png)
- Cladebreaker OFF
  ![Strain Identifier Cladebreaker OFF Workflow](data/workflow/Thresher_strain_identifier_cladebreaker_off_workflow.png)
- Redo Endpoint
  ![Strain Identifier Redo Endpoint Workflow](data/workflow/Thresher_strain_identifier_redo_endpoint_workflow.png)
  Elements highlighted in green (#91a01e) indicate updates applied to existing results.
- New SNPs 
  ![Strain Identifier New SNPs Workflow](data/workflow/Thresher_strain_identifier_new_snps_workflow.png)
  Elements highlighted in green (#91a01e) indicate updates applied to existing results.
- New Full 
  ![Strain Identifier New Full Workflow](data/workflow/Thresher_strain_identifier_new_full_workflow.png)
  Elements highlighted in green (#91a01e) indicate updates applied to existing results.

### Genome Profiler
![Genome Profiler Workflow](data/workflow/Thresher_genome_profiler_workflow.png)

### Evolution Simulator
![Evolution Simulator Workflow](data/workflow/Thresher_evolution_simulator_workflow.png)

## Documentation
### Getting Started
- [System Requirements](docs/requirements.md)
- [Installation](docs/installation.md)
- [Quick Start](docs/quick_start.md)
### Core Functionalities
- [Strain Identifier](docs/usage_strain_identifier.md)
  - [Full Pipeline Mode](docs/usage_strain_identifier_full_pipeline.md)
  - [Cladebreaker OFF Mode](docs/usage_strain_identifier_cladebreaker_off.md)
  - [Redo Endpoint Mode](docs/usage_strain_identifier_redo_endpoint.md)
  - [New SNPs Mode](docs/usage_strain_identifier_new_snps.md)
  - [New Full Mode](docs/usage_strain_identifier_new_full.md)
- [Genome Profiler](docs/usage_genome_profiler.md)
- [Evolution Simulator](docs/usage_evolution_simulator.md)
### Advanced Usage
- [Conda Environment Management](docs/conda_environment_management.md)
- [Resuming Interrupted Runs](docs/resuming_interrupted_runs.md)


## **Development Status**
Please be advised that this pipeline is still in its early development stage. It is subject to significant changes in terms of options, outputs, and other functionalities. Users should be prepared for potential modifications and updates in future releases.

## Authors

**Qianxuan (Sean) She**

## Acknowledgments
Developed under the mentorship of:

**Dr. Paul J. Planet**  
Email: planetp@chop.edu

**Dr. Ahmed M. Moustafa**  
Email: moustafaam@chop.edu

**Dr. Joseph P. Zackular**  
Email: joseph.zackular@pennmedicine.upenn.edu
##
<a href="https://www.research.chop.edu/center-for-microbial-medicine"><img src="data/logo/CMM.png" height="70"></a><a href="https://www.pennmedicine.org/"><img src="data/logo/PennMedicine.png" height="70"></a>