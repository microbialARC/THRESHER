
get_plateau_strains <- function(determine_strains_input,
                                hierarchical_clustering_groups,
                                plateau_length,
                                output_dir){
  
  group_output <- lapply(determine_strains_input, function(group_input){
    
    group_id <- group_input[[1]]$HC_group
    
    # Compare the strain compositions and get an array showing if identical with previous cutoff
    
    rle_array <- sapply(group_input[-1], function(cutoff_composition) {
      
      cutoff <- cutoff_composition$cutoff
      
      get_sorted_composition <- function(target_cutoff) {
        strain_data <- group_input[sapply(group_input, `[[`, "cutoff") == target_cutoff]
        composition <- lapply(strain_data[[1]]$strain_composition, `[[`, "genome")
        composition[order(sapply(composition, length))]
      }
      
      current_composition <- get_sorted_composition(cutoff)
      previous_composition <- get_sorted_composition(cutoff - 1)
      
      if (length(current_composition) != length(previous_composition)) {
        return(0)
      }
      
      current_sorted <- lapply(current_composition, sort)
      previous_sorted <- lapply(previous_composition, sort)
      
      return(if(identical(current_sorted, previous_sorted)) 1 else 0)
    })
    
    names(rle_array) <- sapply(group_input[-1], `[[`, "cutoff")
    
    
    # Find the earliest plateau
    
    plateau <- rle(rle_array)
    plateau_pos <- which(plateau$values == 1 & plateau$lengths >= plateau_length)[1]
    
    plateau_cutoff <- if(!is.na(plateau_pos)){
      if (plateau_pos == 1) as.integer(names(rle_array[1]))
      else as.integer(names(rle_array)[sum(plateau$lengths[1:(plateau_pos-1)]) + 1])
    }else "No Plateau Found"
    
    # Pull out the plateau strains
    
    plateau_strains <- if(plateau_cutoff != "No Plateau Found"){
      group_input[which(sapply(group_input, `[[`, "cutoff") == plateau_cutoff)][[1]]
    }else{
      group_input[1][[1]]
    }
    
    return(list(
      group = group_id,
      plateau = plateau_cutoff,
      plateau_length = plateau_length,
      composition = plateau_strains
    ))
  })
  
  # A data frame describing the plateaus for groups
  
  plateau_df <- do.call(rbind,
                        lapply(group_output,
                               function(group){
                                 data.frame(
                                   group = group$group,
                                   plateau = group$plateau,
                                   plateau_length = plateau_length
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
            file = "plateau_strains.csv")
  
  write.csv(plateau_df,
            quote = FALSE,
            row.names = FALSE,
            file = "group_plateau.csv")
  
  return(list(
    strains = strain_df,
    plateaus = plateau_df,
    plateau_length = plateau_length,
    composition_details = group_output
  ))
  
}
# Libraries
library(dplyr)
# Get the input from Snakemake
determine_strains_input_path <- snakemake@input[["thresher_input"]]
hierarchical_clustering_groups_path <- snakemake@input[["hc_groups"]]
plateau_length <- as.integer(snakemake@params[["plateau_length"]])
output_dir <- snakemake@params[["output_dir"]]
system(paste0("mkdir -p ",output_dir))
setwd(dir = output_dir)
hierarchical_clustering_groups <- readRDS(hierarchical_clustering_groups_path)
determine_strains_input <- readRDS(determine_strains_input_path)

final_strains <- get_plateau_strains(determine_strains_input,
                                     hierarchical_clustering_groups,
                                     plateau_length,
                                     output_dir)
saveRDS(final_strains,
        snakemake@output[["plateau_strains_rds"]])

