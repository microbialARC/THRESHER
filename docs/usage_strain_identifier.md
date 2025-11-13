# Strain Identifier
Determines bacterial clonality and identifies strains/transmission clusters using phylothresholds (phylogenetically-corrected SNP thresholds).

Strain Identifier comprises four modes:
## **[Full Pipeline](docs/Usage_Strain_Identifier_Full_Pipeline.md):** Run the complete analysis pipeline from scratch.
```
thresher strain_identifier full-pipeline -h
```

## **[Redo Endpoint](docs/Usage_Strain_Identifier_Redo_Endpoint.md):** Only rerun the final endpoint analysis using existing intermediate files.
```
thresher strain_identifier redo-endpoint -h
```

  - Current supported endpoints. For each hierarchical clustering group defined in the core gene comphrehensive tree:
    - Plateau: Phylothreshold set at a plateau where further increases no longer change the number or composition of strains within the group.
    - Peak: Phylothreshold set at the peak number of clones defined within the group.
    - Discrepancy:  Phylothreshold set at the point where the discrepancy is minimized within the group.
      - Global:  Phylothreshold set at the first time a global genome is included in any strain within the group.

## **[New SNPs](docs/Usage_Strain_Identifier_New_SNPs.md):** Update existing strain/transmission compositions with new genomes using predefined phylothresholds.
```
thresher strain_identifier new-snps -h
```

## **[New Full](docs/Usage_Strain_Identifier_New_Full.md):** Rerun the full pipeline to update the phylothresholds, and update strain/transmission compositions with new genomes.
```
thresher strain_identifier new-full -h
```