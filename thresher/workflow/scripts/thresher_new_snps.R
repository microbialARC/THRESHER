# Use only SNP to perform single-linkage clustering to update strain/transmission cluster compositions

# Libraries ----
library(parallel)
library(dplyr)
library(purrr)

# Helper function to export the strain compositions ----
export_new_strain <- function(endpoint, new_strain) {
  
  # Define output key patterns based on endpoint
  output_keys <- list(
    rds = paste0("new_", endpoint, "_rds"),
    strains = paste0("new_", endpoint, "_strains_csv"),
    genomes = paste0("new_", endpoint, "_genomes_csv")
  )
  
  # Save RDS file
  saveRDS(new_strain, file = snakemake@output[[output_keys$rds]])
  
  # Helper function to write CSV
  write_clean_csv <- function(data, key) {
    write.csv(data, 
              file = snakemake@output[[key]], 
              row.names = FALSE, 
              quote = FALSE)
  }
  
  # Write CSV files
  write_clean_csv(new_strain$new_strain_compositions, output_keys$strains)
  write_clean_csv(new_strain$new_genomes_df, output_keys$genomes)
}

# Helper function to update strain compositions ----
update_strains <- function(strain_compositions,
                           endpoint,
                           group_thresholds,
                           new_genomes,
                           new_snp_matrix){
  
  ## Get a character to keep track of unassigned new genomes
  new_genomes_left <- new_genomes
  ## Copy a data frame to keep track of updated and complete strain compositions
  new_strain_compositions <- strain_compositions
  ## Create a data frame to keep track of only new genomes
  new_genomes_df <- data.frame()
  
  # Update strains with the match within the threshold ----
  ## sort() to make sure we handle the strains within the same hierarchical cluster in a consistent order
  for(strain_entry in sort(unique(strain_compositions$strain_id))){
    # Find the genomes in the current strain
    strain_genomes <- strain_compositions$genome[strain_compositions$strain_id == strain_entry]
    
    # Subset the SNP matrix to only include comparisons involving the strain genomes
    strain_matrix <- new_snp_matrix[new_snp_matrix$subject %in% strain_genomes |
                                      new_snp_matrix$query %in% strain_genomes,]
    
    # Get the threshold for the current strain group
    strain_group <- strsplit(strain_entry,
                             split = "\\_")[[1]][1]
    
    strain_threshold <- group_thresholds[[endpoint]][group_thresholds$group == strain_group]
    
    
    strain_threshold <- if(is.na(strain_threshold) || 
                           strain_threshold == "Singletons" || 
                           strain_threshold == Inf) {
      NA
    }else{
      as.integer(strain_threshold)
    }
    
    # The additional step for those strains whose groups have no defined threshold
    if(is.na(strain_threshold)){
      
      # If the group containing this strain has no defined threshold
      # which means this group is a group with only this strain and this strain is singleton
      # Then use the SNP distance to iteratively find the closest group to this strain until the closest group has the defined threshold
      
      # Get the order of closest original genomes to this strain
      # The codes below are not the most elegant since I split the logic into multiple sapply functions
      # but it self-explains the logic
      # which is important for the future maintenance by other graduate students or postdocs
      
      # Find the closest genomes to this strain genomes
      closest_genomes <- unlist(sapply(order(strain_matrix$gsnp),
                                       function(row_index){
                                         
                                         setdiff(c(strain_matrix$subject[row_index],
                                                   strain_matrix$query[row_index]),
                                                 c(strain_genomes,
                                                   new_genomes))
                                         
                                       }))
      
      closest_groups <- sapply(closest_genomes,
                               function(genome_entry){
                                 genome_strain <- strain_compositions$strain_id[strain_compositions$genome == genome_entry]
                                 genome_group <- strsplit(genome_strain,
                                                          split = "\\_")[[1]][1]
                               })
      
      closest_threshold <- na.omit(sapply(closest_groups,
                                          function(group_entry){
                                            group_threshold <- group_thresholds[[endpoint]][group_thresholds$group == group_entry]
                                            
                                            if(is.na(group_threshold) || 
                                               group_threshold == "Singletons" || 
                                               group_threshold == Inf) {
                                              NA
                                            }else{
                                              as.integer(group_threshold)
                                            }
                                          }))
      
      # The final threshold is the first non-NA threshold in the closest_threshold vector
      strain_threshold <- as.integer(closest_threshold[1])
      
      # Remove the closest_genomes, closest_groups, closest_threshold variables to avoid these variables being used in the next iteration
      rm(closest_genomes,
         closest_groups,
         closest_threshold)
    }
    
    # Now go back to the strain_matrix to find new genomes within the threshold
    
    new_strain_genomes <- strain_matrix %>% 
      filter(strain_matrix$gsnp <= strain_threshold) %>%
      select(subject,query) %>%
      unlist() %>%
      unique()
    
    # Remove strain_threshold variable to avoid this variable being used in the next iteration
    rm(strain_threshold)
    # If there are new genomes within the threshold, update the strain compositions
    # Otherwise if there is no new genome within the threshold, skip to the next strain
    if(length(new_strain_genomes) > 0){
      if(!(all(new_strain_genomes %in% strain_genomes) & 
           all(strain_genomes %in% new_strain_genomes))){
        # Also make sure the new genomes are in the new_genomes_left
        # otherwise, they have already been assigned to other strains
        new_strain_genomes <- setdiff(new_strain_genomes,
                                      strain_genomes)
        
        new_genomes_left <- setdiff(new_genomes_left, new_strain_genomes)
        
        # Update the strain compositions data frame
        # Only do this if there are new genomes assigned to this strain
        if(length(new_strain_genomes) > 0){
          new_strain_compositions <- rbind(new_strain_compositions,
                                           data.frame(
                                             strain_id = strain_entry,
                                             genome = new_strain_genomes
                                           )) 
          
          new_genomes_df <- rbind(new_genomes_df,
                                  data.frame(
                                    strain_id = strain_entry,
                                    category = "existing_strain",
                                    genome = new_strain_genomes
                                  ))
        }else{
          next
        }
      }
    }else{
      next
    }
  }
  # Update strains with NO defined threshold ----
  # Handle the remaining new genomes that are not assigned to any existing strains
  # the rationale is that if there is any new strain formed by the new genomes only,
  # the new strain will be from new genomes left or the new genomes left would be the missing intermediate to link the genomes in the previous analysis to form the strain
  
  if(length(new_genomes_left) > 0){
    # If there are still new genomes left unassigned to any existing strains
    # Get a list to keep track of redundant strains formed by new_genomes_left only to perform single linkage clustering
    redundant_strain_list <- list()
    # Assign each remaining new genome to a new strain ID
    # Again, without a core genome tree, we use SNP distance to find the closest hierarchical group with the defined threshold
    # But this time, we find the closest genome using the new_snp_matrix only involving the new genomes left
    for(genome_entry in new_genomes_left){
      
      closest_genomes <- new_snp_matrix %>%
        filter((subject == genome_entry | query == genome_entry) &
                 !(subject %in% new_genomes[-which(new_genomes == genome_entry)]) &
                 !(query %in% new_genomes[-which(new_genomes == genome_entry)])) %>%
        arrange(gsnp) %>%
        select(subject, query) %>%
        unlist() %>%
        setdiff(genome_entry)
      
      
      
      closest_groups <- sapply(closest_genomes,
                               function(genome_entry){
                                 genome_strain <- strain_compositions$strain_id[strain_compositions$genome == genome_entry]
                                 genome_group <- strsplit(genome_strain,
                                                          split = "\\_")[[1]][1]
                               })
      
      closest_threshold <- na.omit(sapply(closest_groups,
                                          function(group_entry){
                                            group_threshold <- group_thresholds[[endpoint]][group_thresholds$group == group_entry]
                                            
                                            if(is.na(group_threshold) || 
                                               group_threshold == "Singletons" || 
                                               group_threshold == Inf) {
                                              NA
                                            }else{
                                              as.integer(group_threshold)
                                            }
                                          }))
      
      # The final threshold is the first non-NA threshold in the closest_threshold vector
      genome_threshold <- as.integer(closest_threshold[1])
      
      # Find all single links using the genome_threshold
      genome_cutoff <- new_snp_matrix %>% 
        filter((subject == genome_entry | query == genome_entry) &
                 gsnp <= genome_threshold) %>%
        select(subject, query) %>%
        unlist() %>%
        unique()
      
      if(length(genome_cutoff) == 0){
        genome_cutoff <- genome_entry
      }
      
      redundant_strain_list <- c(redundant_strain_list,
                                 list(genome_cutoff))
      
    }
    
    # Now perform single-linkage clustering on the redundant_strain_list
    # to form new strains from new_genomes_left only
    # Initialize a vector to keep track of assigned genomes
    new_strain_id <- 1
    unique_new_strain_df <- data.frame(strain_id = NA,
                                       genome = new_genomes_left)
    
    
    for(genome_entry in new_genomes_left){
      
      # Find all redundant strains involving the genome in redundant_strain_list
      redundant_entry <- which(sapply(redundant_strain_list,function(redundant_list) genome_entry %in% redundant_list))
      
      # the genomes in redundant strains with redundant_entry are now considered the same strain 
      unique_strain_genomes <- unique(as.character(unlist(sapply(redundant_entry,
                                                                 function(entry) redundant_strain_list[[entry]]))))
      
      # Check if any of unique_strain_genomes is already labeled with strain_id
      # if all are labeled, no action needed 
      # first look for the previous strain label
      
      existing_strain_id <- unique(na.omit(new_strain_compositions$strain[new_strain_compositions$genome %in% unique_strain_genomes]))
      
      if(length(existing_strain_id) == 0){
        
        not_labeled_unique_new_strain_genomes <- unique_new_strain_df$genome[unique_new_strain_df$genome %in% unique_strain_genomes & 
                                                                               is.na(unique_new_strain_df$strain_id)]
        
        if(length(not_labeled_unique_new_strain_genomes) > 0){
          if(identical(sort(unique_strain_genomes),
                       sort(not_labeled_unique_new_strain_genomes))){
            # if none are labeled, label all genomes in unique_strain_genomes with new_strain_id, and add 1 to new_strain_id
            unique_new_strain_df$strain_id[unique_new_strain_df$genome %in% unique_strain_genomes] <- new_strain_id
            new_strain_id <- new_strain_id + 1
          }else{
            # If any one of the genome is already labelled a with strain_id
            # the rest of the unlabeled genomes will be labeled the same strain_id
            # If there are multiple strain_id, use the smallest id
            smallest_strain_id <- as.integer(min(na.omit(unique_new_strain_df$strain_id[unique_new_strain_df$genome %in% unique_strain_genomes])))
            unique_new_strain_df$strain_id[unique_new_strain_df$genome %in% unique_strain_genomes] <- smallest_strain_id
          }
        }
        
      }else{
        
        unique_new_strain_df$strain_id[unique_new_strain_df$genome %in% unique_strain_genomes] <- existing_strain_id
      }
    }
    
    
    # Now append the unique_new_strain_df to new_strain_compositions
    new_strain_compositions <- rbind(new_strain_compositions,
                                     unique_new_strain_df)
    # Also append to new_genomes_df
    new_genomes_df <- rbind(new_genomes_df,
                            unique_new_strain_df %>%
                              transmute(
                                strain_id = strain_id,
                                category = sapply(strain_id,
                                                  function(id_entry){
                                                    if(grepl("_",id_entry)){
                                                      "existing_strain"
                                                    }else{
                                                      "new_strain"
                                                    }
                                                  }),
                                genome = genome
                              )
    )
  }
  # Return the endpoint method, new_strain_compositions, and new_genomes_df
  return(list(endpoint = endpoint,
              new_strain_compositions = new_strain_compositions,
              new_genomes_df = new_genomes_df))
}

# Function to get the new strain compositions  ----
update_thresher_new_snps <- function(new_metadata_path,
                                     original_metadata_path,
                                     output_dir,
                                     thresher_output_dir,
                                     ncores){
  ## Input ----
  # Read new_metadata
  new_metadata <- read.delim(new_metadata_path,
                             header = FALSE,
                             sep = "\t")
  
  # Read original_metadata
  original_metadata <- read.delim(original_metadata_path,
                                  header = FALSE,
                                  sep = "\t")
  
  # Read the SNP matrix
  new_snp_matrix <- readRDS(new_snp_matrix_path)
  ## Plateau strains ----
  plateau_strains <- readRDS(file.path(thresher_output_dir,"thresher","output","plateau_strains.RDS"))
  ## Peak strains ----
  peak_strains <- readRDS(file.path(thresher_output_dir,"thresher","output","peak_strains.RDS"))
  ## Global strains ----
  global_strains <- readRDS(file.path(thresher_output_dir,"thresher","output","global_strains.RDS"))
  ## Discrepancy strains ----
  discrepancy_strains <- readRDS(file.path(thresher_output_dir,"thresher","output","discrepancy_strains.RDS"))
  
  # Update strain compositions ----
  ## Plateau strains ----
  new_plateau_strains <- update_strains(strain_compositions = plateau_strains$strains,
                                        endpoint = "plateau",
                                        group_thresholds = plateau_strains$plateaus,
                                        new_genomes = new_metadata$V1,
                                        new_snp_matrix = new_snp_matrix)
  ## Peak strains ----
  new_peak_strains <- update_strains(strain_compositions = peak_strains$strains,
                                     endpoint = "peak",
                                     group_thresholds = peak_strains$peaks,
                                     new_genomes = new_metadata$V1,
                                     new_snp_matrix = new_snp_matrix)
  ## Global strains ----
  new_global_strains <- update_strains(strain_compositions = global_strains$strains,
                                       endpoint = "global",
                                       group_thresholds = global_strains$global,
                                       new_genomes = new_metadata$V1,
                                       new_snp_matrix = new_snp_matrix)
  
  ## Discrepancy strains ----
  new_discrepancy_strains <- update_strains(strain_compositions = discrepancy_strains$strains,
                                            endpoint = "discrepancy",
                                            group_thresholds = discrepancy_strains$discrepancy,
                                            new_genomes = new_metadata$V1,
                                            new_snp_matrix = new_snp_matrix)
  
  return(
    list(
      new_plateau_strains = new_plateau_strains,
      new_peak_strains = new_peak_strains,
      new_global_strains = new_global_strains,
      new_discrepancy_strains = new_discrepancy_strains
    )
  )
  
}

# Import input from Snakemake
new_metadata_path <- snakemake@params[["new_metadata"]]
original_metadata_path <- snakemake@params[["original_metadata"]]
output_dir <- snakemake@params[["output_dir"]]
thresher_output_dir <- snakemake@params[["thresher_output"]]
new_snp_matrix_path <- snakemake@params[["new_snp_matrix"]]
ncores = snakemake@threads

# Create the output directories if not exist
if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# Execute the function  ----
new_thresher_strains <- update_thresher_new_snps(new_metadata_path = new_metadata_path ,
                                                 original_metadata_path = original_metadata_path,
                                                 output_dir = output_dir,
                                                 thresher_output_dir = thresher_output_dir,
                                                 ncores = ncores)

# Save the updated strain compositions and export results using the export_new_strain function ----
export_new_strain(endpoint = "plateau",
                  new_strain = new_thresher_strains$new_plateau_strains)
export_new_strain(endpoint = "peak",
                  new_strain = new_thresher_strains$new_peak_strains)
export_new_strain(endpoint = "global",
                  new_strain = new_thresher_strains$new_global_strains)
export_new_strain(endpoint = "discrepancy",
                  new_strain = new_thresher_strains$new_discrepancy_strains)
