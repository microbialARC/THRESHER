# Analysis of the output of blastx ----
library(dplyr)
# Function to determine MRSA
analyze_blastx_output <- function(raw_dir,
                                  summary_dir,
                                  strains){
  
  output_list <- list.files(path = raw_dir,
                            pattern = "blastx_*",
                            all.files = TRUE,
                            full.names = TRUE,
                            recursive = FALSE)
  
  mrsa_genome_df <- do.call(rbind,
                            lapply(output_list, 
                                   function(output){
                                     
                                     genome <- gsub("_blastx_mrsa.tsv","",basename(output))
                                     
                                     output_df <- read.table(output,
                                                             header = FALSE)
                                     
                                     colnames(output_df) <- c("qseqid",
                                                              "sseqid",
                                                              "pident",
                                                              "length",
                                                              "mismatch",
                                                              "gapopen",
                                                              "qstart",
                                                              "qend",
                                                              "sstart",
                                                              "send",
                                                              "evalue",
                                                              "bitscore")
                                     #only keep the hits with highest bitscore of different subject seq id
                                     output_df <- output_df %>% group_by(sseqid) %>% slice_max(order_by = bitscore,
                                                                                               n=1,
                                                                                               with_ties = FALSE) %>% ungroup()
                                     
                                     #length of subject id
                                     #MecA:
                                     #WP_057521704.1:677
                                     #WP_063852670.1: 677
                                     #WP_063852677.1: 677
                                     #WP_063852683.1: 677
                                     #WP_000721309.1: 677
                                     #WP_000721306.1: 677
                                     #WP_063852626.1: 677
                                     #WP_063852617.1: 677
                                     #WP_000721310.1: 677
                                     #WP_063851348.1: 677
                                     #MecB:
                                     #WP_012655867.1: 683
                                     #WP_063852710.1: 683
                                     #MecC:
                                     #WP_000725529.1: 674
                                     output_df$coverage <- sapply(X = 1:nrow(output_df),
                                                                  function(X){
                                                                    if(output_df$sseqid[X] == "WP_012655867.1" | output_df$sseqid[X] == "WP_063852710.1"){
                                                                      (output_df$length[X] / 683) * 100
                                                                    }else if(output_df$sseqid[X] == "WP_000725529.1"){
                                                                      (output_df$length[X] / 674) * 100
                                                                    }else{
                                                                      (output_df$length[X] / 677) * 100
                                                                    }
                                                                  })
                                     
                                     
                                     #determine MRSA/MSSA of the isolates
                                     if(any(output_df$coverage >= 70 & output_df$pident >= 70)){
                                       category <- "MRSA"
                                     }else{
                                       category <- "MSSA"
                                     }
                                     
                                     return(data.frame(
                                       genome = genome,
                                       MRSA = category
                                     ))
                                   }))
  
  
  mrsa_strain_df <- do.call(rbind,
                            lapply(sort(unique(strains$strain_id)),
                                   function(strain){
                                     genomes <- strains$genome[strains$strain_id == strain]
                                     strain_category <- names(which.max(table(mrsa_genome_df$MRSA[mrsa_genome_df$genome %in% genomes])))
                                     return(
                                       data.frame(
                                         strain = strain,
                                         MRSA = strain_category
                                       )
                                     )
                                   })) %>%
    arrange(strain)
  
  
  
  write.csv(mrsa_strain_df,
            file.path(summary_dir,"blastx_MRSA_strains.csv"),
            quote = FALSE,
            row.names = FALSE)
  
  write.csv(mrsa_genome_df,
            file.path(summary_dir,"blastx_MRSA_genomes.csv"),
            quote = FALSE,
            row.names = FALSE)
  
}

# Import from snakemake
raw_dir <- snakemake@params[["raw_dir"]]
summary_dir <- snakemake@params[["summary_dir"]]
plateau_strains_rds_path <- snakemake@input[["plateau_strains_rds"]]
peak_strains_rds_path <- snakemake@input[["peak_strains_rds"]]
global_strains_rds_path <- snakemake@input[["global_strains_rds"]]
discrepancy_strains_rds_path <- snakemake@input[["discrepancy_strains_rds"]]

endpoint <- snakemake@params[["endpoint"]]

if(endpoint == "plateau"){
  strains <- readRDS(plateau_strains_rds_path)[["strains"]]
}else if(endpoint == "peak"){
  strains <- readRDS(peak_strains_rds_path)[["strains"]]
}else if(endpoint == "global"){
  strains <- readRDS(global_strains_rds_path)[["strains"]]
}else if(endpoint == "discrepancy"){
  strains <- readRDS(discrepancy_strains_rds_path)[["strains"]]
}

analyze_blastx_output(raw_dir = raw_dir,
                      summary_dir = summary_dir,
                      strains = strains)
