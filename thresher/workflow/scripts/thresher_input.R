# Prepare thresher input for strain determination ----
# Libraries ----
library(parallel)
library(TreeTools)
library(ape)
library(cluster)
library(dplyr)
library(purrr)
# Function to get the thresher input ----
get_thresher_input <- function(hierarchical_clustering_groups_path,
                               group_tree_dir,
                               study_snp_matrix_path,
                               study_global_snp_matrix_path,
                               ncores){
  ## Input for the function ----
  ###  SNP-distance matrix
  study_snp_matrix <- do.call(rbind,
                              lapply(study_snp_matrix_path,
                                     function(matrix_path){
                                       readRDS(matrix_path)
                                     }))

  ### Read HC groups RDS
  hierarchical_clustering_groups <- readRDS(hierarchical_clustering_groups_path)

  ### Get the list of trees from group_tree_path
  group_trees <- list.files(path = group_tree_dir,
                            pattern = "\\.contree$",
                            all.files = TRUE,
                            full.names = TRUE,
                            recursive = FALSE)

  # Process each hierarchical clustering group in hierarchical_clustering_groups
  # and return the corresponding thresher input list

  thresher_input <- mclapply(hierarchical_clustering_groups,
                             function(group_entry){
                               
                               group_id <- group_entry$hc_group
                               group_genomes <- group_entry$genomes
                               # Exclude genomes over gsnp limit 
                               group_genomes <- setdiff(group_genomes,
                                                        group_entry$genomes_overlimit)
                               
                               # theoretically the group tree will have path:
                               group_tree_path <- file.path(group_tree_dir,paste0("Group",group_id,".contree"))
                               # IQtree: It makes no sense to perform bootstrap with less than 4 sequences.
                               # Thus group with less than 4 genomes will not have the tree
                               # In this case no phylogenetic correction will be applied
                               if_phylogenetic_correction <- file.exists(group_tree_path)
                               
                               group_snp_matrix <- study_snp_matrix[study_snp_matrix$subject %in% group_genomes &
                                                                      study_snp_matrix$query %in% group_genomes,]
                               
                               # Summarize group tree if group with no less than 4 genomes
                               if(if_phylogenetic_correction){
                                 group_tree_newick <- read.tree(group_tree_path)
                                 group_tree_newick <- phytools::midpoint_root(group_tree_newick)
                                 group_tree_newick <- Preorder(group_tree_newick)
                                 group_tree_newick_all_nodes <- (length(group_tree_newick$tip.label) + 1):(length(group_tree_newick$tip.label) + group_tree_newick$Nnode)
                                 
                                 group_tree_newick_sum <- lapply(seq_along(group_tree_newick_all_nodes),
                                                                 function(node_entry){
                                                                   sub_group_tree_newick <- Subtree(group_tree_newick,
                                                                                                    group_tree_newick_all_nodes[node_entry])
                                                                   
                                                                   return(list(
                                                                     node = group_tree_newick_all_nodes[node_entry],
                                                                     bootstrap_support = ifelse(
                                                                       node_entry == 1,
                                                                       "Root",
                                                                       as.numeric(group_tree_newick$node.label[node_entry])
                                                                     ),
                                                                     genomes = sub_group_tree_newick$tip.label
                                                                   ))
                                                                 }) 
                               }
                               
                               # Test gsnp from min to 500 ----
                               # If the minimal snp distance in the group is larger than 500, every genome in the group is singleton
                               
                               if(ceiling(min(group_snp_matrix$gsnp)) >= 500){
                                 
                                 
                                 strain_composition <- lapply(seq_along(group_genomes),
                                                              function(genome_entry){
                                                                list(
                                                                  strain_id = paste0(group_id,"_",genome_entry),
                                                                  category = 'singleton',
                                                                  genome = group_genomes[genome_entry],
                                                                  correction = FALSE,
                                                                  bootstrap_support = 0
                                                                )
                                                              })
                                 
                                 
                                 hc_group_input <- list(
                                   list(
                                     HC_group = group_id,
                                     cutoff = ceiling(min(group_snp_matrix$gsnp)),
                                     discrepancy = NA,
                                     before_correction_singletons = length(group_genomes),
                                     before_correction_clones = 0,
                                     clones_corrected = 0,
                                     after_correction_singletons = length(group_genomes),
                                     after_correction_clones = 0,
                                     median_strain_bootstrap_support = 0,
                                     mean_strain_bootstrap_support = 0,
                                     strain_composition = strain_composition
                                   )
                                 )
                                 
                                 
                               }else{
                                 
                                 hc_group_input <- lapply((ceiling(min(group_snp_matrix$gsnp)):500),
                                                          function(gsnp_cutoff){
                                                            
                                                            ### SNP only ----
                                                            #Find all single links among the study genomes at this gsnp cutoff
                                                            group_snp_matrix_cutoff <- group_snp_matrix[group_snp_matrix$gsnp <= gsnp_cutoff,]
                                                            # Within those links, the genomes in group_snp_matrix_cutoff would be linked to at least one another genome and considered part of a strain
                                                            # Which means that those genomes will not be singletons
                                                            strain_genomes <- unique(c(group_snp_matrix_cutoff$subject,
                                                                                       group_snp_matrix_cutoff$query))
                                                            
                                                            # Find singleton genomes
                                                            singleton_genomes <- c(setdiff(group_genomes,strain_genomes),
                                                                                   group_entry$genomes_overlimit)
                                                            
                                                            # Summrize the single-linkage 
                                                            redundant_strain_list <- lapply(seq_along(strain_genomes),
                                                                                            function(genome_entry){
                                                                                              genome <- strain_genomes[genome_entry]
                                                                                              # Find the links associated with this genome
                                                                                              # Every genomes involved in the links associated with this genome will be considered from the same strain
                                                                                              genome_link <- group_snp_matrix_cutoff[group_snp_matrix_cutoff$subject == genome | 
                                                                                                                                       group_snp_matrix_cutoff$query == genome,]
                                                                                              strain_linked_genomes <- unique(c(genome_link$subject,
                                                                                                                                genome_link$query))
                                                                                              return(
                                                                                                list(id = genome_entry,
                                                                                                     genomes = strain_linked_genomes)
                                                                                              )
                                                                                            })
                                                            
                                                            
                                                            
                                                            # Perform single-linkage clustering with redundant_strain_list
                                                            # Strain_id_index is used to track the strain ID, starting with 1
                                                            
                                                            strain_id_index <- 1
                                                            
                                                            # Use both unique_strain_df to perform single-linkage clustering
                                                            
                                                            unique_strain_df <- data.frame(genome = strain_genomes,
                                                                                           category = "clone",
                                                                                           strain_id = NA)
                                                            
                                                            for(genome in strain_genomes){
                                                              # Find all redundant strains involving the genome in redundant_strain_list
                                                              redundant_entry <- which(sapply(redundant_strain_list,function(redundant_list) genome %in% redundant_list$genomes))
                                                              # the genomes in redundant strains with redundant_entry are now considered the same strain 
                                                              unique_strain_genomes <- unique(as.character(unlist(sapply(redundant_entry,
                                                                                                                         function(entry) redundant_strain_list[[entry]]$genomes))))
                                                              # Check if any of unique_strain_genomes is already labeled with strain_id
                                                              # if all are labeled, no action needed 
                                                              not_labeled_unique_strain_genomes <- unique_strain_df$genome[unique_strain_df$genome %in% unique_strain_genomes & 
                                                                                                                             is.na(unique_strain_df$strain_id)]
                                                              
                                                              if(length(not_labeled_unique_strain_genomes) > 0){
                                                                if(identical(sort(unique_strain_genomes),
                                                                             sort(not_labeled_unique_strain_genomes))){
                                                                  # if none are labeled, label all genomes in unique_strain_genomes with strain_id_index, and add 1 to strain_id_index
                                                                  unique_strain_df$strain_id[unique_strain_df$genome %in% unique_strain_genomes] <- strain_id_index
                                                                  strain_id_index <- strain_id_index + 1
                                                                }else{
                                                                  # If any one of the genome is already labelled a with strain_id
                                                                  # the rest of the unlabeled genomes will be labeled the same strain_id
                                                                  # If there are multiple strain_id, use the smallest id
                                                                  smallest_strain_id <- min(na.omit(unique_strain_df$strain_id[unique_strain_df$genome %in% unique_strain_genomes]))
                                                                  unique_strain_df$strain_id[unique_strain_df$genome %in% unique_strain_genomes] <- smallest_strain_id
                                                                }
                                                              }
                                                            }
                                                            
                                                            if(length(singleton_genomes) > 0){
                                                              # Add singletons to unique_strain_df if there is any
                                                              unique_strain_df <- rbind(unique_strain_df, 
                                                                                        data.frame(genome = singleton_genomes, 
                                                                                                   category = "singleton",
                                                                                                   strain_id = seq_along(singleton_genomes) + max(unique_strain_df$strain_id)))
                                                            }
                                                            
                                                            # Count before-correction singleton/clones
                                                            before_correction_clones <- length(unique(unique_strain_df$strain_id[unique_strain_df$category == "clone"]))
                                                            before_correction_singletons <- length(unique(unique_strain_df$strain_id[unique_strain_df$category == "singleton"]))
                                                            
                                                            
                                                            if(if_phylogenetic_correction){
                                                              ### Phylogenetic Tree ----
                                                              # Only process clones (strains that have no less than 2 genomes, not singleton strains)
                                                              # Find the strain number of the strain with no less than 2 genomes (clones)
                                                              # Add bootstrap support for the strains
                                                              unique_strain_df <- unique_strain_df %>% mutate(correction = NA,
                                                                                                              bootstrap_support = NA)
                                                              
                                                              unique_strain_df$bootstrap_support[unique_strain_df$category == "singleton"] <- 0
                                                              
                                                              unique_strain_tree <- lapply(sort(unique(unique_strain_df$strain_id[unique_strain_df$category == "clone"])),
                                                                                           function(strain_id){
                                                                                             # strain genomes determined only with gsnp
                                                                                             strain_snp_genomes <- unique_strain_df$genome[unique_strain_df$strain_id == strain_id]
                                                                                             #iterate group_tree_newick_sum to find the clades covering all strain_snp_genomes
                                                                                             strain_clade <- group_tree_newick_sum[sapply(group_tree_newick_sum, function(clade) 
                                                                                               all(strain_snp_genomes %in% clade$genomes))] %>%
                                                                                               #the best final clade is the one with least number of genomes 
                                                                                               {.[which.min(sapply(., function(clade) length(clade$genomes)))]}
                                                                                             
                                                                                             #return result to the unique_strain_tree
                                                                                             list(
                                                                                               strain_id = strain_id,
                                                                                               snp_only_genomes = sort(strain_snp_genomes),
                                                                                               snp_tree_genomes = sort(strain_clade[[1]]$genomes),
                                                                                               bootstrap_support = strain_clade[[1]]$bootstrap_support,
                                                                                               node = strain_clade[[1]]$node
                                                                                             )
                                                                                           })
                                                              
                                                              ### Corrected for phylogenetic structure -----
                                                              
                                                              correction_strain_tree <- unique_strain_tree[sapply(unique_strain_tree, function(strain) !identical(strain$snp_only_genomes, strain$snp_tree_genomes))]
                                                              # Count discrepancy
                                                              discrepancy <- sum(unlist(sapply(correction_strain_tree,
                                                                                               function(strain){
                                                                                                 length(setdiff(strain$snp_tree_genomes[!grepl("GCA_|GCF_",strain$snp_tree_genomes)],
                                                                                                                strain$snp_only_genomes))
                                                                                               })))
                                                              # Only perform correction when there is discrepancy
                                                              if(length(correction_strain_tree) > 0){
                                                                # Filter the list to keep only the strains with discrepancy between tree and gsnp
                                                                # And only accept the correction only if the bootstrap support >= 70 (https://doi.org/10.1093/sysbio/42.2.182)
                                                                # Or this is the root of the tree
                                                                
                                                                correction_strain_tree <- correction_strain_tree[sapply(correction_strain_tree, function(strain) ((strain$bootstrap_support >= 70) | (strain$bootstrap_support == "Root")))]
                                                                clones_corrected <- length(correction_strain_tree)
                                                                # The function to perform Cladebreaker to return the list of Cladebreaker-divided strains
                                                                correction_strain_list <- lapply(correction_strain_tree,
                                                                                                 function(strain_clade){
                                                                                                   # 2 situations are expected:
                                                                                                   if(!any(grepl("GCA_|GCF_",strain_clade$snp_tree_genomes))){
                                                                                                     # 1st situation: only discrepancy genomes with no breaker genomes
                                                                                                     # use the minimal ID for all snp_tree_genomes in unique_strain_df
                                                                                                     return(
                                                                                                       list(
                                                                                                         list(original_strain = strain_clade$strain_id,
                                                                                                              corrected_genomes =  strain_clade$snp_tree_genomes,
                                                                                                              category = "clone",
                                                                                                              bootstrap_support = as.integer(strain_clade$bootstrap_support)
                                                                                                         ))
                                                                                                     )
                                                                                                   }else{
                                                                                                     # 2nd situation: breaker genomes were introduced (with or without discrepancy genomes)
                                                                                                     
                                                                                                     strain_clade_newick <- Subtree(group_tree_newick,
                                                                                                                                    strain_clade$node)
                                                                                                     
                                                                                                     strain_clade_newick_all_nodes <- (length(strain_clade_newick$tip.label) + 1):(length(strain_clade_newick$tip.label) + strain_clade_newick$Nnode)
                                                                                                     # From strain_clade_newick_all_nodes
                                                                                                     # We need to find the biggest clades in strain_clade_newick without being broken by breaker genomes
                                                                                                     # And the genomes in those biggest pure-study-genomes clades are from the same strains
                                                                                                     tmp_list <- lapply(seq_along(strain_clade_newick_all_nodes),
                                                                                                                        function(node_entry){
                                                                                                                          
                                                                                                                          node <- strain_clade_newick_all_nodes[node_entry]
                                                                                                                          sub_strain_clade_newick <- Subtree(strain_clade_newick,
                                                                                                                                                             node)
                                                                                                                          #Check if this clade has only study genomes
                                                                                                                          if(!any(grepl("GCA_|GCF_",sub_strain_clade_newick$tip.label))){
                                                                                                                            # If there are only study genomes
                                                                                                                            # proceed to check if this is the biggest pure-study-genome clade (There are breaker genomes in the parent clade)
                                                                                                                            # Parent node:
                                                                                                                            parent_node <- strain_clade_newick$edge[which(strain_clade_newick$edge[,2] == node), 1]
                                                                                                                            parent_sub_strain_clade_newick <- Subtree(strain_clade_newick,
                                                                                                                                                                      parent_node)
                                                                                                                            if(any(grepl("GCA_|GCF_",parent_sub_strain_clade_newick$tip.label))){
                                                                                                                              
                                                                                                                              return(
                                                                                                                                list(
                                                                                                                                  original_strain = strain_clade$strain_id,
                                                                                                                                  corrected_genomes = sub_strain_clade_newick$tip.label,
                                                                                                                                  category = "clone",
                                                                                                                                  bootstrap_support = as.integer(strain_clade_newick$node.label[node_entry])
                                                                                                                                )
                                                                                                                              )
                                                                                                                            }
                                                                                                                          }
                                                                                                                        }) %>% purrr::compact()
                                                                                                     
                                                                                                     #Add the singletons to the list
                                                                                                     singletons <- setdiff(strain_clade$snp_tree_genomes[!grepl("GCA_|GCF_", strain_clade$snp_tree_genomes)],
                                                                                                                           unlist(sapply(tmp_list, `[[`, "corrected_genomes")))
                                                                                                     
                                                                                                     #Add the singleton list to the tmp_list
                                                                                                     if(length(singletons) > 0){
                                                                                                       tmp_list <- c(tmp_list,lapply(singletons, function(singleton) list(
                                                                                                         original_strain = strain_clade$strain_id,
                                                                                                         corrected_genomes = singleton,
                                                                                                         category = "singleton",
                                                                                                         bootstrap_support = 0
                                                                                                       )))
                                                                                                     }
                                                                                                     return(tmp_list)
                                                                                                     
                                                                                                   }
                                                                                                 })
                                                                
                                                                
                                                                # Use correction_strain_list to update info in unique_strain_df
                                                                for(correction_strain in correction_strain_list){
                                                                  for(i in seq_along(correction_strain)){
                                                                    if(length(correction_strain) == 1){
                                                                      unique_strain_df$strain_id[unique_strain_df$genome %in% correction_strain[[i]]$corrected_genomes] <- correction_strain[[i]]$original_strain
                                                                    }else{
                                                                      unique_strain_df$strain_id[unique_strain_df$genome %in% correction_strain[[i]]$corrected_genomes] <- paste0(correction_strain[[i]]$original_strain,
                                                                                                                                                                                  "_",
                                                                                                                                                                                  i)
                                                                    }
                                                                    
                                                                    unique_strain_df$category[unique_strain_df$genome %in% correction_strain[[i]]$corrected_genomes] <- correction_strain[[i]]$category
                                                                    unique_strain_df$correction[unique_strain_df$genome %in% correction_strain[[i]]$corrected_genomes] <- TRUE
                                                                    unique_strain_df$bootstrap_support[unique_strain_df$genome %in% correction_strain[[i]]$corrected_genomes] <- correction_strain[[i]]$bootstrap_support
                                                                  }
                                                                }
                                                                
                                                              }else{
                                                                clones_corrected <- 0
                                                              }
                                                              
                                                              
                                                              unique_strain_df$correction[is.na(unique_strain_df$correction)] <- FALSE
                                                              
                                                              # Fill in the bootstrap support
                                                              for(id in unique(unique_strain_df$strain_id[is.na(unique_strain_df$bootstrap_support)])){
                                                                
                                                                if(suppressWarnings((is.null(unique_strain_tree[[as.integer(id)]]$bootstrap_support)))){
                                                                  unique_strain_df$bootstrap_support[unique_strain_df$strain_id == id] <- 0
                                                                }else{
                                                                  unique_strain_df$bootstrap_support[unique_strain_df$strain_id == id] <- unique_strain_tree[[as.integer(id)]]$bootstrap_support 
                                                                }
                                                              }
                                                              
                                                              # Count before-correction singleton/clones
                                                              after_correction_clones <- length(unique(unique_strain_df$strain_id[unique_strain_df$category == "clone"]))
                                                              after_correction_singletons <- length(unique(unique_strain_df$strain_id[unique_strain_df$category == "singleton"]))
                                                              
                                                              
                                                              ### Return result ----
                                                              # Sort the output using unique_strain_df
                                                              strain_composition <- lapply(sort(unique(unique_strain_df$strain_id)),
                                                                                           function(id){
                                                                                             list(
                                                                                               strain_id = id,
                                                                                               category = unique(unique_strain_df$category[unique_strain_df$strain_id == id]),
                                                                                               genome = unique_strain_df$genome[unique_strain_df$strain_id == id],
                                                                                               correction = unique(unique_strain_df$correction[unique_strain_df$strain_id == id]),
                                                                                               bootstrap_support = unique(unique_strain_df$bootstrap_support[unique_strain_df$strain_id == id])
                                                                                             )
                                                                                           })
                                                              
                                                              # Calculate the mean/median bootstrap support for strains
                                                              
                                                              
                                                              median_strain_bootstrap_support <- unlist(median(sapply(strain_composition[sapply(strain_composition, function(x) x$category == "clone")], `[[`, "bootstrap_support")))
                                                              
                                                              median_strain_bootstrap_support <- ifelse(is.null(median_strain_bootstrap_support), 
                                                                                                        0,
                                                                                                        ifelse(is.na(as.numeric(median_strain_bootstrap_support)),
                                                                                                               0,
                                                                                                               as.numeric(median_strain_bootstrap_support)))
                                                              
                                                              mean_strain_bootstrap_support <- mean(sapply(strain_composition[sapply(strain_composition, function(x) x$category == "clone")], `[[`, "bootstrap_support"))
                                                              
                                                              mean_strain_bootstrap_support <- ifelse(is.null(mean_strain_bootstrap_support), 
                                                                                                      0,
                                                                                                      ifelse(is.na(as.numeric(mean_strain_bootstrap_support)),
                                                                                                             0,
                                                                                                             as.numeric(mean_strain_bootstrap_support)))
                                                            }else{
                                                              
                                                              # Sort the output using unique_strain_df but with bootstrap_support or correction as NA
                                                              strain_composition <- lapply(sort(unique(unique_strain_df$strain_id)),
                                                                                           function(id){
                                                                                             list(
                                                                                               strain_id = id,
                                                                                               category = unique(unique_strain_df$category[unique_strain_df$strain_id == id]),
                                                                                               genome = unique_strain_df$genome[unique_strain_df$strain_id == id],
                                                                                               correction = NA,
                                                                                               bootstrap_support = NA
                                                                                             )
                                                                                           })
                                                              clones_corrected <- NA
                                                              after_correction_singletons <- NA
                                                              after_correction_clones <- NA
                                                              median_strain_bootstrap_support <- NA 
                                                              mean_strain_bootstrap_support <- NA
                                                              discrepancy <- NA
                                                            }
                                                            
                                                            
                                                            return(
                                                              list(
                                                                HC_group = group_id,
                                                                cutoff = gsnp_cutoff,
                                                                discrepancy = discrepancy,
                                                                before_correction_singletons = before_correction_singletons,
                                                                before_correction_clones = before_correction_clones,
                                                                clones_corrected = clones_corrected,
                                                                after_correction_singletons = after_correction_singletons,
                                                                after_correction_clones = after_correction_clones,
                                                                median_strain_bootstrap_support = median_strain_bootstrap_support,
                                                                mean_strain_bootstrap_support = mean_strain_bootstrap_support,
                                                                strain_composition = strain_composition
                                                              )
                                                            )
                                                          })
                               }
                             },
                             mc.cores = ncores)

                             return(thresher_input)
}

# Get the input from Snakemake
ncores <- snakemake@threads
hierarchical_clustering_groups_path <- snakemake@input[["hc_groups"]]
group_tree_dir <- snakemake@params[["group_tree_dir"]]
study_snp_matrix_path <- snakemake@input[["study_snp_matrix"]]
study_global_snp_matrix_path <- snakemake@input[["global_snp_matrix"]]
output_dir <- snakemake@params[["thresher_input_dir"]]

# Execute the function to get thresher input
setwd(dir = output_dir)

thresher_input <- get_thresher_input(hierarchical_clustering_groups_path = hierarchical_clustering_groups_path,
                                     group_tree_dir = group_tree_dir,
                                     study_snp_matrix_path = study_snp_matrix_path,
                                     study_global_snp_matrix_path = study_global_snp_matrix_path,
                                     ncores = ncores)

# Export the output back to Snakemake
saveRDS(thresher_input, snakemake@output[["thresher_input"]])