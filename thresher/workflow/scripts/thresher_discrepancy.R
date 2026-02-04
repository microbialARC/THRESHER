get_discrepancy_strains <- function(thresher_input,
                                    hierarchical_clustering_groups,
                                    output_dir,
                                    singleton_threshold,
                                    ncores){
  
  group_output <- mclapply(thresher_input, function(group_input){
    
    group_id <- group_input[[1]]$HC_group
    # After the peak at the beginning of the threshold 
    # Find the threshold that have minimal discrepancy in this group
    # Find the phylothreshold using discrepancy endpoint before singleton threshold
    
    if(length(group_input) == 1){
    # If length(group_input) == 1, skip because the discrepancy function does not work here
    group_result <- list(
        group = group_id,
        discpreancy = NA,
        composition = list(strain_composition = group_input[[1]]$strain_composition)
      )
    
    }else{
      # Only consider the cutoff range below or equal to singleton threshold
      peak_discrepancy <- group_input[[which.max(sapply(group_input, \(x) x$discrepancy * (x$cutoff <= singleton_threshold)))]]$cutoff
      discrepancy_cutoff <- group_input[[which.min(sapply(group_input, \(x) if(x$cutoff > peak_discrepancy) x$discrepancy else Inf))]]$cutoff
      # Pull the discrepancy strains
      discrepancy_strains <- group_input[which(sapply(group_input, `[[`, "cutoff") == discrepancy_cutoff)][[1]]
      group_result <- list(
        group = group_id,
        discpreancy = discrepancy_cutoff,
        composition = discrepancy_strains
      )
    }
    return(group_result)
  },
  mc.cores = ncores)
  
  # A data frame describing the discrepancy for groups
  
  discrepancy_df <- do.call(rbind,
                     lapply(group_output,
                            function(group){
                              data.frame(
                                group = group$group,
                                discrepancy = sapply(group$discpreancy,
                                                     function(disc_entry){
                                                       if(is.na(disc_entry)){
                                                         "Singletons"
                                                       }else{
                                                         as.integer(disc_entry)
                                                       }
                                                     })
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
library(parallel)
# Get the input from Snakemake
thresher_input_path <- snakemake@input[["thresher_input"]]
hierarchical_clustering_groups_path <- snakemake@input[["hc_groups"]]
output_dir <- snakemake@params[["output_dir"]]
ncores = snakemake@threads
singleton_threshold = as.integer(snakemake@params[["singleton_threshold"]])
system(paste0("mkdir -p ",output_dir))
setwd(dir = output_dir)
hierarchical_clustering_groups <- readRDS(hierarchical_clustering_groups_path)
thresher_input <- readRDS(thresher_input_path)

final_strains <- get_discrepancy_strains(thresher_input,
                                         hierarchical_clustering_groups,
                                         output_dir,
                                         singleton_threshold,
                                         ncores)

saveRDS(final_strains,
        snakemake@output[["discrepancy_strains_rds"]])