get_peak_strains <- function(determine_strains_input,
                             hierarchical_clustering_groups,
                             output_dir){
  
  group_output <- lapply(determine_strains_input, function(group_input){
    
    group_id <- group_input[[1]]$HC_group
    # Find the threshold that generates most clones in this group
    peak_cutoff <- min(group_input[[which.max(sapply(group_input, \(x) x$after_correction_clones))]]$cutoff)
    # Pull the peak strains
    peak_strains <- group_input[which(sapply(group_input, `[[`, "cutoff") == peak_cutoff)][[1]]
   
    return(list(
      group = group_id,
      peak = peak_cutoff,
      composition = peak_strains
    ))
  })
  
  # A data frame describing the peaks for groups
  
  peak_df <- do.call(rbind,
                        lapply(group_output,
                               function(group){
                                 data.frame(
                                   group = group$group,
                                   peak = group$peak
                                 )
                               })) %>%
    arrange(group)
  
  # Use group_output to generate the final strain data frame
  
  # Groups >= 2 genomes 
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
  
  # Groups with only 1 genome
  singleton_groups <-  hierarchical_clustering_groups[sapply(hierarchical_clustering_groups,
                                                             function(group){
                                                               length(group$genomes) == 1
                                                             })]
  strain_df <- rbind(strain_df,
                     do.call(rbind,
                             lapply(singleton_groups,
                                    function(singleton_group){
                                      return(
                                        data.frame(
                                          strain_id = paste0(singleton_group$hc_group,
                                                             "_1"),
                                          genome = unlist(singleton_group$genomes)
                                        )
                                      )
                                    })))
  
  write.csv(strain_df,
            quote = FALSE,
            row.names = FALSE,
            file = "peak_strains.csv")
  
  write.csv(peak_df,
            quote = FALSE,
            row.names = FALSE,
            file = "group_peak.csv")
  
  return(list(
    strains = strain_df,
    peaks = peak_df,
    composition_details = group_output
  ))
}

#Libraries
library(dplyr)
# Get the input from Snakemake
determine_strains_input_path <- snakemake@input[["thresher_input"]]
hierarchical_clustering_groups_path <- snakemake@input[["hc_groups"]]
output_dir <- snakemake@params[["output_dir"]]
system(paste0("mkdir -p ",output_dir))
setwd(dir = output_dir)
hierarchical_clustering_groups <- readRDS(hierarchical_clustering_groups_path)
determine_strains_input <- readRDS(determine_strains_input_path)

final_strains <- get_peak_strains(determine_strains_input,
                                  hierarchical_clustering_groups,
                                  output_dir)
saveRDS(final_strains,
        snakemake@output[["peak_strains_rds"]])

