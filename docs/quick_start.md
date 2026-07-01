# Quick Start

Run THRESHER's core functions using default settings. Each command requires minimal input to get started.

## Strain Identifier
To run the full pipeline to identify strains and detect transmission clusters:
```
thresher full --metadata /path/to/metadata.txt --species <species>
```

[Example metadata file](example/example_metadata.txt)

How do I know the Genbank accession for my genomes?

See the [Genbank Accession](genbank_accession.md) page.

Note:
- The default thread count is 1, which may result in lengthy runtimes. It is highly recommended to increase the thread count using the `-t/--threads` option to improve performance, regardless of dataset size.

- Bakta genome annotation runs `{threads}` parallel jobs, each requiring approximately 10 GB of RAM. Ensure your system has sufficient memory to support the requested thread count.

For detailed usage, configuration options, and advanced features, see the [THRESHER](docs/usage_thresher_overview.md) documentation.
