# Hierarchical clustering of the everything tree to perform the clustering of the genomes
# Libraries ----
#packages to analyze the trees
library(TreeTools)
library(ape)
#clustering tool
library(cluster)
#others 
library(dplyr)
library(purrr)

# Get the input from Snakemake
comprehensive_tree_path <- snakemake@input[["comprehensive_tree_path"]]
study_snp_matrix_path <- snakemake@input[["study_snp_matrix_path"]]
mlst_results_path <- snakemake@input[["mlst_results_path"]]
output_dir <- snakemake@params[["thresher_input_dir"]]
system(paste0("mkdir -p ",output_dir))
# Function for hierarchical clustering ----
# Hierarchical Clustering ----
#A CC-free way to group the study genomes in the everything-tree using the pairwise distances
hierarchical_clustering <- function(comprehensive_tree_path,
                                    study_snp_matrix_path,
                                    mlst_results_path){
  
  # read the SNP matrix of study genomes comparing to study genomes to perform the QC
  study_snp_matrix <- readRDS(study_snp_matrix_path)
  # read mlst results
  mlst_results <- read.csv(mlst_results_path,sep = "\t")
  # read the tree
  comprehensive_tree <- ape::read.tree(comprehensive_tree_path)
  # mid-point rooting of the tree
  comprehensive_tree <- phytools::midpoint_root(comprehensive_tree)
  # pre-order the tree
  comprehensive_tree <- Preorder(comprehensive_tree)
  # summarize the tree by iterating through nodes in the tree
  all_nodes <- (length(comprehensive_tree$tip.label) + 1):(length(comprehensive_tree$tip.label) + comprehensive_tree$Nnode)
  
  comprehensive_tree_sum <- lapply(seq_along(all_nodes),
                                function(node_entry){
                                  sub_comprehensive_tree <- Subtree(comprehensive_tree,all_nodes[node_entry])
                                  
                                  if(comprehensive_tree$node.label[node_entry] == "Root"){
                                    
                                    return(list(
                                      node = all_nodes[node_entry],
                                      SHaLRT_support = "Root",
                                      bootstrap_support = "Root",
                                      genomes = sub_comprehensive_tree$tip.label
                                      ))
                                    
                                  }else{
                                    
                                    return(list(
                                      node = all_nodes[node_entry],
                                      SHaLRT_support = as.numeric(strsplit(comprehensive_tree$node.label[node_entry],
                                                                split = "/")[[1]][1]),
                                      bootstrap_support = as.numeric(strsplit(comprehensive_tree$node.label[node_entry],
                                                                   split = "/")[[1]][2]),
                                      genomes = sub_comprehensive_tree$tip.label
                                    ))
                                  }
                                })
  
  # The typing groups(CC, Clade...) mapped to everything tree
  # After excluding the "Unassigned" group:
  # Condition1: if there is only one typing group
  # and this only one typing group covers all genomes in the tree
  # (At smallest clade covers all genomes of this group is the root)
  condition1 <- if(sum(unique(unlist(mlst_results[[3]])) != "Unassigned", na.rm=TRUE) == 1){
    # the genomes of identified groups
    condition1_genomes <- mlst_results$genome[mlst_results[[3]] != "Unassigned"]
    #iterate comprehensive_tree_sum to find the clades covering all genomes
    condition1_clade <- comprehensive_tree_sum[sapply(comprehensive_tree_sum, function(clade) 
      all(condition1_genomes %in% clade$genomes))] %>%
      #the best final clade is the one with least number of genomes 
      {.[which.min(sapply(., function(clade) length(clade$genomes)))]}
    
    if(condition1_clade[[1]]$bootstrap_support == "Root"){
      TRUE
    }else{
      FALSE
    }
    
  }else{
    FALSE
  }

  # Condition2: the pair-wise snp distances of genomes in the tree are all within 500 gsnp
  # If condition1 and 2 are met, the Hierarchical Clustering will NOT be performed
  # thus all genomes are in the same group proceeding to the subsequent fine-grained analysis
  condition2 <- max(study_snp_matrix$gsnp) <= 500
  
  
  if(condition1 & condition2){
    
    # Hierarchical Clustering will NOT be performed
    final_group <- setNames(rep(1,length(comprehensive_tree$tip.label)),
                            comprehensive_tree$tip.label)
    
  }else{
    # Hierarchical Clustering will be performed
    # matrix of pairwise distances from tree
    distance_matrix <-  cophenetic.phylo(comprehensive_tree)
    # the various methods can become options in the future
    hclust_result <- hclust(as.dist(distance_matrix),
                            method = "complete")
    
    # we test the range from 2 to (number of genomes in the everything - 1) clustering patterns 
    hclust_result_silhouette_scores <- as.numeric(unlist(sapply((2:(length(comprehensive_tree$tip.label)-1)),
                                                                function(clustering_pattern_entry){
                                                                  clusters <- cutree(hclust_result, k = clustering_pattern_entry)
                                                                  return(mean(silhouette(clusters,distance_matrix)[, "sil_width"]))
                                                                })))
    
    final_group <- cutree(hclust_result,
                          k = (which.max(hclust_result_silhouette_scores) + 1))
    
  }
  
  #Look at each group to check the pair-wise SNP distance
  #If the genome in the group has no link with other genomes within 1000 gsnps (can be changed)
  #The genome is a singletons
  
  hc_groups_qc <- lapply(sort(unique(final_group)),
                         function(hc_group){
                           
                           group_genomes <- names(final_group)[final_group == hc_group]
                           
                           if(length(group_genomes) > 1){
                             
                             # Find the SHaLRT/bootstrap support of this HC group
                             sum_entry <- which(sapply(comprehensive_tree_sum,
                                                       function(summary){identical(sort(summary$genomes),
                                                                                   sort(group_genomes))}))
                             
                             group_SHaLRT_support <- comprehensive_tree_sum[[sum_entry]]$SHaLRT_support
                             group_bootstrap_support <- comprehensive_tree_sum[[sum_entry]]$bootstrap_support
                             
                             # If all genomes over 1000gsnp limit? 
                             group_study_snp_matrix <- study_snp_matrix[(study_snp_matrix$subject %in% group_genomes) &
                                                                          (study_snp_matrix$query %in% group_genomes),]
                             
                             if(min(group_study_snp_matrix$gsnp) >= 1000){
                               
                               group_all_overlimit <- TRUE
                               genomes_overlimit <- group_genomes
                               
                             }else{
                               
                               group_all_overlimit <- FALSE
                               
                               # Iterate through group_genomes to find which genomes has >= 1000 gsnp
                               genomes_overlimit <- as.character(unlist(sapply(group_genomes,
                                                                  function(genome){
                                                                    group_study_snp_matrix_genome <- group_study_snp_matrix[group_study_snp_matrix$subject == genome |
                                                                                                                              group_study_snp_matrix$query == genome,]
                                                                    if(min(group_study_snp_matrix_genome$gsnp) >= 1000){
                                                                      return(genome)
                                                                    }
                                                                  })))
                             }
                             
                             
                             return(list(
                               hc_group = hc_group,
                               genomes = group_genomes,
                               all_overlimit = group_all_overlimit,
                               genomes_overlimit = genomes_overlimit,
                               SHaLRT_support = group_SHaLRT_support,
                               bootstrap_support = group_bootstrap_support
                             ))
                             
                           }else{
                             # If there is only 1 group genome: singleton
                             return(list(
                               hc_group = hc_group,
                               genomes = group_genomes,
                               all_overlimit = FALSE,
                               genomes_overlimit = NULL,
                               SHaLRT_support = "N/A",
                               bootstrap_support = "N/A"
                             ))
                           }
                         })
  
  
  hc_groups_simplified <- do.call(rbind,
                                  lapply(hc_groups_qc,
                                         function(group){
                                           data.frame(
                                             genome = group$genomes,
                                             group = group$hc_group,
                                             overlimit = sapply(group$genomes, 
                                                                function(genome) if(genome %in% group$genomes_overlimit) TRUE else FALSE)
                                           )
                                         }))
  #return hc_groups_qc full and simplified csv
  return(list("full" = hc_groups_qc,
              "simplified " = hc_groups_simplified))
}

hierarchical_clustering_groups <- hierarchical_clustering(comprehensive_tree_path,
                                                          study_snp_matrix_path,
                                                          mlst_results_path)

# Export the output back to Snakemake
saveRDS(hierarchical_clustering_groups[["full"]],
        snakemake@output[["hc_groups"]])
# Export the simplified output as csv to Snakemake
write.csv(hierarchical_clustering_groups$`simplified `,
          snakemake@output[["hc_groups_csv"]],
          quote = FALSE,
          row.names = FALSE)
