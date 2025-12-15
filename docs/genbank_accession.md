# GenBank Accessions

Before running THRESHER, check whether input genome assemblies already exist in public databases. If they do, THRESHER will exclude the public versions of the input genomes from the "global genomes" used in downstream analyses. This exclusion is required for two features: the global endpoint method and the `--use_cladebreaker` option. CladeBreaker assumes your strains are circulating locally rather than globally and is enabled by default; to disable it, use `--use_cladebreaker False`.

## Option 1: No GenBank Accession

If the sequencing reads or assemblies have not been submitted to public databases like NCBI or ENA, there are no GenBank accessions yet. Just enter `new` in the second column of your metadata file.

## Option 2: Look Up via NCBI Website

If the sequencing reads or assemblies have already been submitted to public databases, you may find your GenBank accessions through the web interface:

1. Go to [NCBI](https://www.ncbi.nlm.nih.gov/)
2. Search for your BioProject or BioSample ID
3. Navigate to the "Assembly details" section to find the accession numbers

## Option 3: Retrieve with Entrez Direct using BioSample IDs

If there's no "Assembly details" section on the BioProject/BioSample page, you can fetch accessions using NCBI's [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) command-line tools.

### Installing EDirect

The easiest method is via conda:
```bash
conda create -n entrez-direct -c bioconda entrez-direct
conda activate entrez-direct
```

### Fetching Accessions

To retrieve the GenBank accession for a BioSample ID (e.g., `SAMN04456594`):
```bash
esearch -db biosample -query "SAMN04456594" | elink -target assembly | esummary | xtract -pattern DocumentSummary -element AssemblyAccession
```

This returns `GCA_016534985.2`. Replace the BioSample ID with your own to get the corresponding accessions, and put the results in the second column of your input metadata file.

## Option 4: Retrieve with R package [Rentrez](https://github.com/ropensci/rentrez/)
You can also use the R package `rentrez` to fetch GenBank accessions in R.
A custom function to do this is provided below:
### Install Rentrez in R
```R
library(devtools)
install_github("ropensci/rentrez")
```
### R function to get GenBank accession from BioSample ID
```R
library(rentrez)
get_genbank_accession <- function(biosample_acc) {
  tryCatch({
    search <- entrez_search(db = "biosample", term = biosample_acc)
    if (length(search$ids) == 0) return(NA)
    
    links <- entrez_link(dbfrom = "biosample", id = search$ids, db = "assembly")
    assembly_id <- links$links$biosample_assembly
    
    if (is.null(assembly_id)) return(NA)
    
    summary <- entrez_summary(db = "assembly", id = assembly_id)
    return(summary$assemblyaccession)
  }, error = function(e) NA)
}
```

### Fetch accession for a BioSample ID
```R
get_genbank_accession("SAMN04456594")
# [1] "GCA_016534985.2"
```