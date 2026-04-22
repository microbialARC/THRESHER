get_global_strains <- function(thresher_input,
                               hierarchical_clustering_groups,
                               global_snp_matrix,
                               threshold_ceiling,
                               output_dir,
                               ncores){
  
  group_output <- mclapply(thresher_input, function(group_input){
    
    group_id <- group_input[[1]]$HC_group
    # there are 3 conditions to tell whether this group has only singletons
    # 1. Only 1 entry in the group_input
    condition_1 <- length(group_input) == 1
    # 2. Cutoff is NA for the entry
    condition_2 <- all(sapply(group_input, function(x) is.na(x$cutoff)))
    # 3. All strain compositions are singletons (i.e. only one genome in each composition)
    condition_3 <- all(sapply(group_input, function(x) all(sapply(x$strain_composition, function(comp) length(comp$genome) == 1))))
    
    if(condition_1 & condition_2 & condition_3){
      return(list(
        group = group_id,
        global = NA,
        composition = group_input[[1]]
      ))
    }
    # Find genomes in the group
    group_genomes <- unlist(sapply(hierarchical_clustering_groups, \(x) x$genomes[x$hc_group == group_id]))
    group_global_snp_matrix <- global_snp_matrix[global_snp_matrix$subject %in% group_genomes,]
    
    # Steps to determine the global cutoff for the group:
    # For each study genome, find the minimum SNP distance to any global genome as the global threshold for the group.
    # Since the definition of global is the smallest SNP distance of any study genome in the group to a global genome
    global_cutoff <- if(round(min(group_global_snp_matrix$gsnp)) >= threshold_ceiling) threshold_ceiling else round(min(group_global_snp_matrix$gsnp))
    # and then we want to ensure that the global cutoff is not below the largest cutoff among study genomes in the group
    global_cutoff <- if(global_cutoff < min(sapply(group_input, `[[`, "cutoff"))) min(sapply(group_input, `[[`, "cutoff")) else global_cutoff
    # However, we also want to ensure that the global cutoff does not exceed the threshold ceiling
    global_cutoff <- if(global_cutoff > threshold_ceiling) threshold_ceiling else global_cutoff
    
    # (Optional): If there is only 1 composition in the group 
    # That is, the smallest SNP distance among study genomes in the group is above the singleton threshold
    # In this case, we set the global cutoff to be the cutoff of that single composition
    if(length(group_input) == 1){
      global_cutoff <- group_input[[1]]$cutoff
    } 
    # Find the index of the composition in the group that has the same cutoff as the global cutoff after the above adjustments.
    global_index <- which(sapply(group_input, `[[`, "cutoff") == global_cutoff)
    # However, if there is still no global_index found after the above adjustments, use the maximum cutoff among the compositions in the group that is less than the global cutoff
    if(length(global_index) == 0){
      global_index <- which.max(sapply(group_input, `[[`, "cutoff")[sapply(group_input, `[[`, "cutoff") < global_cutoff])
    }
    # Use the global cutoff to determine the global strains for the group
    global_strains <- group_input[[global_index]]
    
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
threshold_ceiling <- as.integer(snakemake@params[["threshold_ceiling"]])
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
                                    threshold_ceiling,
                                    output_dir,
                                    ncores)

saveRDS(final_strains,
        snakemake@output[["global_strains_rds"]])

