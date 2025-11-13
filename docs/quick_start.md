# Quick Start

Run THRESHER's core functions using default settings. Each command requires minimal input to get started.

## Strain Identifier
To run the full pipeline to identify strains and detect transmission clusters:
```
thresher strain_identifier full-pipeline --metadata /path/to/metadata.txt --species <species>
```

[Example metadata file](example/example_metadata.txt)

## Genome Profiler
```
thresher genome_profiler --input_genome /path/to/genome.fasta --species <species>
```

## Evolution Simulator
*Coming soon.*

---

For detailed usage, configuration options, and advanced features, see the [Strain Identifier](usage_strain_identifier.md), [Genome Profiler](usage_genome_profiler.md), and [Evolution Simulator (*Coming Soon.*)](usage_evolution_simulator.md) documentation pages.