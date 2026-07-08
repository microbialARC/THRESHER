# Quick Start

This guide is written for a first-time run of THRESHER. Since new users have no prior THRESHER output to build on, **full** mode is the only mode to start with. It executes the complete pipeline from raw input through strain identification and transmission cluster detection. (The remaining modes support specialized and continued workflows)

Running **full** mode with default settings requires only two inputs: the path to your metadata file and the species under analysis.

```
thresher full --metadata /path/to/metadata.txt --species <species>
```

[Example metadata file](example/example_metadata.txt)

**How do I find the Genbank accession for my genomes?**

See the [Genbank Accession](genbank_accession.md) page.

**Notes**

- The default thread count is 1, which may result in lengthy runtimes. Increasing the thread count with the `-t/--threads` option is strongly recommended to improve performance, regardless of dataset size.

- Bakta genome annotation runs `{threads}` parallel jobs, each requiring approximately 10 GB of RAM. Ensure your system has sufficient memory to support the requested thread count.

For detailed usage, configuration options, and advanced features, see the [THRESHER](docs/usage_thresher_overview.md) documentation.