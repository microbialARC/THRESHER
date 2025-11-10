# Analysis of the output of blastx ----
library(dplyr)
# Function to determine MRSA
analyze_blastx_output_new_full <- function(plateau_strains_rds_path,
                                           peak_strains_rds_path,
                                           global_strains_rds_path,
                                           discrepancy_strains_rds_path,
                                           thresher_output_dir,
                                           raw_dir,
                                           summary_dir){
  # Make the summary directory if it does not exist
  if(!dir.exists(summary_dir)){
    dir.create(summary_dir,
               recursive = TRUE)
  }
  # Analyze blastx output for new genomes ----
  new_output_list <- list.files(path = raw_dir,
                                pattern = "blastx_*",
                                all.files = TRUE,
                                full.names = TRUE,
                                recursive = FALSE)
  
  new_mrsa_genome_df <- do.call(rbind,
                                lapply(new_output_list, 
                                       function(output_entry){
                                         
                                         genome <- gsub("_blastx_mrsa.tsv","",basename(output_entry))
                                         
                                         output_df <- read.table(output_entry,
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
                                         if(any(output_df$pident >= 70 & output_df$coverage >= 70)){
                                           category <- "MRSA"
                                         }else{
                                           category <- "MSSA"
                                         }
                                         
                                         return(data.frame(
                                           genome = genome,
                                           MRSA = category
                                         ))
                                       }))
  
  # Combine with original genome results
  
  original_genome_results <- read.csv(file.path(thresher_output_dir,"blastx","mrsa","output","summary","blastx_MRSA_genomes.csv"))
  
  updated_genome_results <- rbind(new_mrsa_genome_df,
                                  original_genome_results)
  
  # Export the updated genome results
  
  write.csv(updated_genome_results,
            file.path(summary_dir,"blastx_MRSA_genomes.csv"),
            quote = FALSE,
            row.names = FALSE)
  
  # Analyze MRSA status for strains using the updated genome results ----
  
  plateau_strains <- readRDS(plateau_strains_rds_path)[["strains"]]
  peak_strains <- readRDS(peak_strains_rds_path)[["strains"]]
  global_strains <- readRDS(global_strains_rds_path)[["strains"]]
  discrepancy_strains <- readRDS(discrepancy_strains_rds_path)[["strains"]]
  
  # Use the position to loop through the different strain compositions
  # and get endpoint method
  
  method_compositions <- list(
    list(endpoint = "plateau", composition = plateau_strains),
    list(endpoint = "peak", composition = peak_strains),
    list(endpoint = "global", composition = global_strains),
    list(endpoint = "discrepancy", composition = discrepancy_strains)
  )
  
  for(pos_idx in seq_along(method_compositions)){
    endpoint <- method_compositions[[pos_idx]]$endpoint
    updated_strain_composition <- method_compositions[[pos_idx]]$composition
    
    updated_mrsa_strain_df <- do.call(rbind,
                                      lapply(sort(unique(updated_strain_composition$strain_id)),
                                             function(strain_id_entry){
                                               genomes <- updated_strain_composition$genome[updated_strain_composition$strain_id == strain_id_entry]
                                               strain_category <- names(which.max(table(updated_genome_results$MRSA[updated_genome_results$genome %in% genomes])))
                                               return(
                                                 data.frame(
                                                   strain = strain_id_entry,
                                                   MRSA = strain_category
                                                 )
                                               )
                                             })) %>%
      arrange(strain)
    
    write.csv(updated_mrsa_strain_df,
              file.path(summary_dir,paste0("blastx_MRSA_",endpoint,"_strains.csv")),
              quote = FALSE,
              row.names = FALSE) 
  }
}

# Import from snakemake
raw_dir <- snakemake@params[["raw_dir"]]
summary_dir <- snakemake@params[["summary_dir"]]
plateau_strains_rds_path <- snakemake@input[["plateau_strains_rds"]]
peak_strains_rds_path <- snakemake@input[["peak_strains_rds"]]
global_strains_rds_path <- snakemake@input[["global_strains_rds"]]
discrepancy_strains_rds_path <- snakemake@input[["discrepancy_strains_rds"]]
thresher_output_dir <- snakemake@params[["thresher_output"]]


analyze_blastx_output_new_full(plateau_strains_rds_path = plateau_strains_rds_path,
                               peak_strains_rds_path = peak_strains_rds_path,
                               global_strains_rds_path = global_strains_rds_path,
                               discrepancy_strains_rds_path = discrepancy_strains_rds_path,
                               thresher_output_dir = thresher_output_dir,
                               raw_dir = raw_dir,
                               summary_dir = summary_dir)

