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
The evolution simulator uses genome profiler outputs to simulate bacterial evolution:

1. Run genome profiler on your reference genome
2. Configure simulation parameters (mutation rates, recombination, gene gain/loss)
3. Run the evolution simulator with the output from genome profiler


---

For detailed usage, configuration options, and advanced features, see the [Strain Identifier](usage_strain_identifier.md), [Genome Profiler](usage_genome_profiler.md), and [Evolution Simulator](usage_evolution_simulator.md) documentation pages.