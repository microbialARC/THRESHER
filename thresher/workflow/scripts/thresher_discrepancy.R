get_discrepancy_strains <- function(determine_strains_input,
                                    hierarchical_clustering_groups,
                                    output_dir){
  
  group_output <- lapply(determine_strains_input, function(group_input){
    
    group_id <- group_input[[1]]$HC_group
    # After the peak at the beginning of the threshold 
    # Find the threshold that have minimal discrepancy in this group
    # Find the peak <= 50 gSNP
    
    peak_discrepancy <- group_input[[which.max(sapply(group_input, \(x) x$discrepancy * (x$cutoff <= 50)))]]$cutoff
    
    discrepancy_cutoff <- group_input[[which.min(sapply(group_input, \(x) if(x$cutoff > peak_discrepancy) x$discrepancy else Inf))]]$cutoff
    
    # Pull the discrepancy strains
    discrepancy_strains <- group_input[which(sapply(group_input, `[[`, "cutoff") == discrepancy_cutoff)][[1]]
    
    return(list(
      group = group_id,
      discpreancy = discrepancy_cutoff,
      composition = discrepancy_strains
    ))
  })
  
  # A data frame describing the discrepancy for groups
  
  discrepancy_df <- do.call(rbind,
                     lapply(group_output,
                            function(group){
                              data.frame(
                                group = group$group,
                                discrepancy = group$discpreancy
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
            file = "discrepancy_strains.csv")
  
  write.csv(discrepancy_df,
            quote = FALSE,
            row.names = FALSE,
            file = "group_discrepancy.csv")
  
  return(list(
    strains = strain_df,
    discrepancy = discrepancy_df,
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

final_strains <- get_discrepancy_strains(determine_strains_input,
                                         hierarchical_clustering_groups,
                                         output_dir)

saveRDS(final_strains,
        snakemake@output[["discrepancy_strains_rds"]])