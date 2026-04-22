get_peak_strains <- function(thresher_input,
                             hierarchical_clustering_groups,
                             output_dir,
                             ncores){
  
  group_output <- mclapply(thresher_input, function(group_input){
    
    group_id <- group_input[[1]]$HC_group
    # Find the threshold that generates most clones in this group
    # First determine whether or not this group has only singletons. 
    # If it does then no peak cutoff and the value is NA

    # there are 3 conditions to tell whether this group has only singletons
    # 1. Only 1 entry in the group_input
    condition_1 <- length(group_input) == 1
    # 2. Cutoff is NA for the entry
    condition_2 <- all(sapply(group_input, function(x) is.na(x$cutoff)))
    # 3. All strain compositions are singletons (i.e. only one genome in each composition)
    condition_3 <- all(sapply(group_input, function(x) all(sapply(x$strain_composition, function(comp) length(comp$genome) == 1))))
    
    if(condition_1 || condition_2 || condition_3){
      peak_cutoff <- NA
      peak_strains <- group_input[[1]]
    } else {
      peak_cutoff <- min(group_input[[which.max(sapply(group_input, \(x) x$after_correction_clones))]]$cutoff)
      # Pull the peak strains
      peak_strains <- group_input[which(sapply(group_input, `[[`, "cutoff") == peak_cutoff)][[1]]
    }
   
    return(list(
      group = group_id,
      peak = peak_cutoff,
      composition = peak_strains
    ))
  },
  mc.cores = ncores)
  
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
library(parallel)
# Get the input from Snakemake
thresher_input_path <- snakemake@input[["thresher_input"]]
hierarchical_clustering_groups_path <- snakemake@input[["hc_groups"]]
output_dir <- snakemake@params[["output_dir"]]
ncores <- snakemake@threads
system(paste0("mkdir -p ",output_dir))
setwd(dir = output_dir)
hierarchical_clustering_groups <- readRDS(hierarchical_clustering_groups_path)
thresher_input <- readRDS(thresher_input_path)

final_strains <- get_peak_strains(thresher_input,
                                  hierarchical_clustering_groups,
                                  output_dir,
                                  ncores)
saveRDS(final_strains,
        snakemake@output[["peak_strains_rds"]])

