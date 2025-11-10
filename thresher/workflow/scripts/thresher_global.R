get_global_strains <- function(thresher_input,
                               hierarchical_clustering_groups,
                               global_snp_matrix,
                               output_dir,
                               ncores){
  
  group_output <- mclapply(thresher_input, function(group_input){
    
    group_id <- group_input[[1]]$HC_group
    
    # Find genomes in the group
    group_genomes <- unlist(sapply(hierarchical_clustering_groups, \(x) x$genomes[x$hc_group == group_id]))
    group_global_snp_matrix <- global_snp_matrix[global_snp_matrix$subject %in% group_genomes,]
    
    # The earliest cutoff a study genome hit a global genome
    # if the min SNP distance is above 500, threshold is set to 500
  
    global_cutoff <- if(round(min(group_global_snp_matrix$gsnp)) >= 500) 500 else round(min(group_global_snp_matrix$gsnp))
    
    # If the global_cutoff < min(sapply(group_input, `[[`, "cutoff"))
    # Use min(sapply(group_input, `[[`, "cutoff"))
    global_cutoff <- if(global_cutoff < min(sapply(group_input, `[[`, "cutoff"))) min(sapply(group_input, `[[`, "cutoff")) else global_cutoff
    # Pull the global strains
    global_strains <- group_input[which(sapply(group_input, `[[`, "cutoff") == global_cutoff)][[1]]
    
    return(list(
      group = group_id,
      global = global_cutoff,
      composition = global_strains
    ))
  },
  mc.cores = ncores)
  
  # A data frame describing the peaks for groups
  
  global_df <- do.call(rbind,
                     lapply(group_output,
                            function(group){
                              data.frame(
                                group = group$group,
                                global = group$global
                              )
                            })) %>%
    arrange(group)
  
  # Use group_output to generate the final strain data frame
  strain_df <- do.call(rbind,
                       lapply(group_output,
                              function(group){
                                
                                group_strain_composition <- lapply(group$composition$strain_composition,
                                                                   `[[`,
                                                                   "genome")
                                return(do.call(rbind,
                                               lapply(seq_along(group_strain_composition),
                                                      function(strain_entry){
                                                        return(
                                                          data.frame(
                                                            strain_id = paste0(group$group,
                                                                               "_",
                                                                               strain_entry),
                                                            genome = unlist(group_strain_composition[[strain_entry]])
                                                          )
                                                        )
                                                      })))
                                
                              }))
  
  write.csv(strain_df,
            quote = FALSE,
            row.names = FALSE,
            file = "global_strains.csv")
  
  write.csv(global_df,
            quote = FALSE,
            row.names = FALSE,
            file = "group_global.csv")
  
  return(list(
    strains = strain_df,
    global = global_df,
    composition_details = group_output
  ))
}

# Get the input from Snakemake
thresher_input_path <- snakemake@input[["thresher_input"]]
hierarchical_clustering_groups_path <- snakemake@input[["hc_groups"]]
global_snp_matrix_path <- snakemake@input[["global_snp_matrix"]]
output_dir <- snakemake@params[["output_dir"]]
ncores <- snakemake@threads
setwd(dir = output_dir)
#Libraries 
library(dplyr)
library(parallel)
hierarchical_clustering_groups <- readRDS(hierarchical_clustering_groups_path)
thresher_input <- readRDS(thresher_input_path)
global_snp_matrix <- readRDS(global_snp_matrix_path)

final_strains <- get_global_strains(thresher_input,
                                    hierarchical_clustering_groups,
                                    global_snp_matrix,
                                    output_dir,
                                    ncores)

saveRDS(final_strains,
        snakemake@output[["global_strains_rds"]])

