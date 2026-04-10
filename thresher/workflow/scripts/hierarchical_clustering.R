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
# Helper function
# Replace any non-alphanumeric characters (except dash/underscore) with underscore in the first column
parse_genome_name <- function(original_name) {
  parsed_name <- gsub("[^a-zA-Z0-9._-]", "_", original_name)
  parsed_name <- gsub("_+", "_", parsed_name)
  if (parsed_name == "Reference") {
    parsed_name <- "Reference_Genome"
  }
  trimws(parsed_name, whitespace = "_")
} 
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
  mlst_results$genome <- sapply(mlst_results$genome,
                               parse_genome_name)
  # read the tree
  comprehensive_tree <- ape::read.tree(comprehensive_tree_path)
  # mid-point rooting of the tree
  comprehensive_tree <- phytools::midpoint_root(comprehensive_tree)
  # pre-order the tree
  comprehensive_tree <- Preorder(comprehensive_tree)
  # replace any non-alphanumeric characters in tip.label
  comprehensive_tree$tip.label <- sapply(comprehensive_tree$tip.label,
                                               parse_genome_name)
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
  
  # The typing groups(CC, Clade...) mapped to comprehensive tree
  # After excluding the "Unassigned" group:
  # Condition1: if there is only one typing group
  # and this only one typing group covers all genomes in the tree
  # (At smallest clade covers all genomes of this group is the root)
  condition1 <- if(sum(unique(unlist(mlst_results[[3]])) != "Unassigned", na.rm=TRUE) == 1){
    # the genomes of identified groups
    condition1_genomes <- mlst_results$genome[mlst_results[[3]] != "Unassigned"]
    #iterate comprehensive_tree_sum to find the clade covering all genomes
    condition1_clade <- comprehensive_tree_sum[sapply(comprehensive_tree_sum, function(clade) 
      all(condition1_genomes %in% clade$genomes))] %>%
      #the best final clade is the one with least number of genomes 
      {.[which.min(sapply(., function(clade) length(clade$genomes)))]}
    # The bootstrap support of this clade should be "Root" because the clade covers all genomes in the tree, which means the clade is the root of the tree
    if(!is.na(condition1_clade[[1]]$bootstrap_support) && condition1_clade[[1]]$bootstrap_support == "Root"){
      TRUE
    }else{
      FALSE
    }
    
  }else{
    FALSE
  }

  # Condition2: Check if there is multiple cluster of pair-wise snp distances within the genomes in the tree if Condition1 is met
  # The hypothesis here is that:
  # If there is just one cluster of pair-wise snp distances determined by local minima of the density estimation of pair-wise snp distances
  # The genomes are closely related enough in terms of single-linkage clustering to not do hierarchical clustering
  # If there are multiple groups
  # This means that there is natural discontinuity among the genomes and hierarchical clustering should be performed to separate them

  # 1000bp is a good balance between sensitivity and specificity for separating closely related genomes

  condition2 <- if(condition1){
    
    # Use 1000 bp bandwidth to estimate the density
    snp_density <- density(study_snp_matrix$gsnp,bw = 1000)
    snp_local_min <- which(diff(sign(diff(snp_density$y))) == 2) + 1
    snp_density_cluster <- length(snp_local_min) + 1
    
    if(snp_density_cluster == 1){
      # If snp_density_cluster is 1, all genomes are clustered together in one group
      # thus hierarchical clustering will NOT be performed
      TRUE
    }else{
      # Otherwise, genomes in the tree are not closely related enough
      # thus hierarchical clustering will be performed
      FALSE
    }
  }else{
    # If condition1 is not met, condition2 is automatically FALSE
    FALSE
  }
  
  if(condition1 & condition2){
    # If condition1 and 2 are met, the Hierarchical Clustering will NOT be performed
    # thus all genomes are in the same group proceeding to the subsequent fine-grained analysis
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
  
  # Correct the hirerarchical clustering groups using tree topology if necessary
  # First step: fix final_group to be topology-consistent
  for(hc_group in sort(unique(final_group))) {
    group_genomes <- names(final_group)[final_group == hc_group]
    if(length(group_genomes) > 1) {
      sum_entry <- which(sapply(comprehensive_tree_sum, function(sum_entry)
        identical(sort(as.character(sum_entry$genomes)), sort(group_genomes))))
      if(length(sum_entry) == 0) {
        fix_sum_entry <- which.min(sapply(comprehensive_tree_sum, function(sum_entry) {
          if(all(group_genomes %in% sum_entry$genomes)) length(sum_entry$genomes) else Inf
        }))
        discrepant_genomes <- setdiff(comprehensive_tree_sum[[fix_sum_entry]]$genomes, group_genomes)
        final_group[discrepant_genomes] <- hc_group
      }
    }
  }

  # Second step: build hc_groups_qc from the corrected final_group
  
  hc_groups_qc <- lapply(sort(unique(final_group)), function(hc_group) {
    
    group_genomes <- names(final_group)[final_group == hc_group]
    
    if (length(group_genomes) > 1) {
      
      # Find the SHaLRT/bootstrap support of this HC group
      sum_entry <- which(sapply(comprehensive_tree_sum, function(summary) {
        identical(sort(as.character(summary$genomes)), sort(group_genomes))
      }))
      
      if (length(sum_entry) == 0) {
        # Find the smallest clade covering all group genomes
        sum_entry <- which.min(sapply(comprehensive_tree_sum, function(summary) {
          if (all(group_genomes %in% summary$genomes)) length(summary$genomes) else Inf
        }))
      }
      
      group_SHaLRT_support <- comprehensive_tree_sum[[sum_entry]]$SHaLRT_support
      group_bootstrap_support <- comprehensive_tree_sum[[sum_entry]]$bootstrap_support
      
      # Check pair-wise SNP distances within group
      group_study_snp_matrix <- study_snp_matrix[
        (study_snp_matrix$subject %in% group_genomes) &
          (study_snp_matrix$query %in% group_genomes), ]
      
      # Look at each group to check the pair-wise SNP distance
      # If the genome in the group has no link with other genomes within 1000 gsnps (can be changed)
      # The genome is a singleton and will be separated from the group
      
      if (min(group_study_snp_matrix$gsnp) >= 1000) {
        group_all_overlimit <- TRUE
        genomes_overlimit <- group_genomes
      } else {
        group_all_overlimit <- FALSE
        genomes_overlimit <- as.character(unlist(sapply(group_genomes, function(genome) {
          genome_snps <- group_study_snp_matrix[
            group_study_snp_matrix$subject == genome |
              group_study_snp_matrix$query == genome, ]
          if (min(genome_snps$gsnp) >= 1000) return(genome)
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
      
    } else {
      # Singleton
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
  # Remove those entries with no genomes (if any) because of the correction/fix of the hierarchical clustering groups
  hc_groups_qc <- hc_groups_qc[sapply(hc_groups_qc, function(group) length(group$genomes) > 0)]
  
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
  # Return hc_groups_qc full and simplified csv
  return(list("full" = hc_groups_qc,
              "simplified" = hc_groups_simplified))
}

hierarchical_clustering_groups <- hierarchical_clustering(comprehensive_tree_path,
                                                          study_snp_matrix_path,
                                                          mlst_results_path)

# Export the output back to Snakemake
saveRDS(hierarchical_clustering_groups[["full"]],
        snakemake@output[["hc_groups"]])
# Export the simplified output as csv to Snakemake
write.csv(hierarchical_clustering_groups[["simplified"]],
          snakemake@output[["hc_groups_csv"]],
          quote = FALSE,
          row.names = FALSE)
