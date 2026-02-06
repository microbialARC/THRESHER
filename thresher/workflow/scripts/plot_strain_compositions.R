# Visualize the strain compositions determined by Thresher using ML tree of each group and SNP distances
# Libraries ----
library(dplyr)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(phytools)
library(ggExtra)
library(TreeTools)
library(aplot)
library(purrr)
# Function ----

# Helper function to calculate PDF dimensions based on number of genomes
get_pdf_dimensions <- function(n_genomes) {
  # Define scaling factors and minimum dimensions
  # Height and width increase with number of genomes to ensure readability
  height_per_genome <- 0.6
  min_height <- 6

  width_per_genome <- 0.6
  min_width <- 12
  
  height <- max(min_height, n_genomes * height_per_genome)
  width <- max(min_width, n_genomes * width_per_genome + 10)
  
  return(list(width = round(width, 1), height = round(height, 1)))
}

get_strain_compositions_plot <- function(endpoint_method,
                                         plateau_strains_rds_path,
                                         peak_strains_rds_path,
                                         discrepancy_strains_rds_path,
                                         global_strains_rds_path,
                                         study_snp_matrix_path,
                                         global_snp_matrix_path,
                                         group_tree_path,
                                         use_cladebreaker,
                                         output_dir){
  
  ## Input ----
  # Set output dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  plot_export_dir <- file.path(output_dir,endpoint_method)
  dir.create(plot_export_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(dir = plot_export_dir)
  
  
  # Combining Study and Global SNP Matrix 
  combined_snp_matrix <- rbind(readRDS(study_snp_matrix_path),
                               readRDS(global_snp_matrix_path))
  
  
  # thresher strain depending on the endpoint method
  
  thresher_strain_path <- switch(
    endpoint_method,
    "plateau" = plateau_strains_rds_path,
    "peak" = peak_strains_rds_path,
    "discrepancy" = discrepancy_strains_rds_path,
    "global" = global_strains_rds_path
  )
  
  thresher_strain <- readRDS(thresher_strain_path)
  thresher_strain_df <- thresher_strain$strains
  
  # Get the optimal phylothreshold for each group depending on the endpoint method
  thresher_phylothresholds_df <- switch(
    endpoint_method,
    "plateau" = thresher_strain$plateaus %>% 
      select(group,plateau) %>%
      rename(group = group,
             phylothreshold = plateau),
    "peak" = thresher_strain$peaks %>% 
      select(group,peak) %>%
      rename(group = group,
             phylothreshold = peak),
    "discrepancy" = thresher_strain$discrepancy %>% 
      select(group,discrepancy) %>%
      rename(group = group,
             phylothreshold = discrepancy),
    "global" = thresher_strain$global %>% 
      select(group,global) %>%
      rename(group = group,
             phylothreshold = global)
  )

  # Get the paths to the tree of each HC group
  hc_tree_list <- readLines(group_tree_path)
  
  # A list to store the plots of each group for current endpoint method
  plot_rds <- list()
  
  for(tree_entry in hc_tree_list){
    
    group_id <- gsub("\\.contree|Group","",basename(tree_entry))
    group_tree_newick <- ggtree::read.tree(tree_entry)
    group_total_genomes <- group_tree_newick$tip.label
    group_tree_newick <- phytools::midpoint_root(group_tree_newick)
    group_tree_newick <- TreeTools::Preorder(group_tree_newick)
    n_tips <- length(group_total_genomes)
    group_tree_newick_all_nodes <- (n_tips + 1):(n_tips + group_tree_newick$Nnode)
    # Get the information of each node in the group_tree_newick
    # The code is reused from thresher_input.R but we don't take the bootstrap support value here
    group_tree_newick_sum <- lapply(group_tree_newick_all_nodes,
                                    function(node_entry){
                                      sub_group_tree_newick <- Subtree(group_tree_newick, node_entry)
                                      return(list(
                                        node = node_entry,
                                        genomes = sub_group_tree_newick$tip.label
                                      ))
                                    })
    # Get the data frame recording the node of the clones (strain with no less than 2 genomes) in this group 
    # As the input for geom_hilight() and geom_cladelab()
    
    thresher_clones <- names(which(table(thresher_strain_df$strain_id[thresher_strain_df$genome %in% group_total_genomes]) > 1))
    
    # Generate annotation data frame mapping phylogenetic tree nodes to THRESHER-identified strains
    # When use_cladebreaker is TRUE, strain compositions exclude global genomes (cladebreakers)
    # reflecting the hypothesis of hyper-local transmission where global genomes should not be considered the same strain
    # with local strains. When FALSE, global genomes are included in strain composition,
    # allowing for the situation where the genomes deposited at NCBI GenBank could also be the same strain with study genomes
    
    group_tree_annotation_df <- do.call(rbind,
                                        lapply(thresher_clones,
                                               function(clone_entry){
                                                 
                                                 clone_genomes <- thresher_strain_df$genome[thresher_strain_df$strain_id == clone_entry]
                                                 
                                                 clone_node <- sapply(group_tree_newick_sum,
                                                                      function(node_entry){
                                                                        
                                                                        node_genomes <- if(use_cladebreaker){
                                                                          node_entry$genomes
                                                                        }else if(!use_cladebreaker){
                                                                          # Exclude global genomes from node only when use_cladebreaker is FALSE
                                                                          node_entry$genomes[!grepl("GCA_|GCF_", node_entry$genomes)]
                                                                        }
                                                                        
                                                                        if(identical(sort(node_genomes),
                                                                                     sort(clone_genomes))){
                                                                          
                                                                          node_entry$node
                                                                        }else{
                                                                          NULL
                                                                        }
                                                                      }) %>%
                                                   compact() %>%
                                                   as.character()
                                                 # If no node found, no data frame is returned
                                                 if(length(clone_node) == 0){
                                                   return(NULL)
                                                 }else{
                                                    return(data.frame(
                                                      node = clone_node,
                                                      strain_id = clone_entry,
                                                      method = tools::toTitleCase(endpoint_method),
                                                      stringsAsFactors = FALSE
                                                    ))
                                                 }
                                               }))
    
    if(!is.null(group_tree_annotation_df)){
      
      if(!use_cladebreaker){
        # If cladebreaker is not used 
        # The same strain could be mapped to multiple nodes due to global genomes inclusion
        # So we pick the node with least genomes
        group_tree_annotation_df <- group_tree_annotation_df %>%
          group_by(strain_id) %>%
          filter(node == max(as.numeric(node))) %>%
          ungroup()
      }
      group_tree_annotation_df$node <- as.numeric(group_tree_annotation_df$node)
      
      # Offset the label 
      group_tree_annotation_df <- group_tree_annotation_df %>%
        group_by(node) %>%
        mutate(
          label_offset = (row_number() - 1) * 0.5,
        ) %>%
        ungroup()
    }
    
    ## Visualization ----
    # Tree
    group_tree <- ggtree(group_tree_newick,
                         layout = "rectangular") +
      geom_treescale(x = 0,
                     y = length(group_tree_newick$tip.label) - 1,
                     width = signif(max(group_tree_newick$edge.length)/10, 1),
                     fontsize = 3.5,
                     linesize = 0.5,
                     offset = -0.5) +
      geom_tiplab(size = 0,
                  as_ylab = TRUE,
                  align = TRUE) +
      geom_tippoint(aes(color = ifelse(grepl("GCA_", label),
                                       "Global",
                                       "Study")),
                    size = 3.5,
                    alpha = 1) +
      scale_color_manual(values = c("Global" = "#786452", "Study" = "#41b6e6")) +
      labs(color = "Genome Category") +
      # Add phylothreshold annotation
      annotate("text",
               x = 0,
               y = length(group_tree_newick$tip.label),
               label = paste0("Optimal Phylothreshold: ", thresher_phylothresholds_df$phylothreshold[thresher_phylothresholds_df$group == group_id]),
               hjust = 0,
               vjust = 1,
               size = 3.5,
               fontface = "bold") +
      {if (!is.null(group_tree_annotation_df) && nrow(group_tree_annotation_df) > 0) {
        c(
          list(
            ggnewscale::new_scale_color(),
            ggnewscale::new_scale_fill()
          ),
          list(
            geom_hilight(
              data = group_tree_annotation_df,
              aes(node = node, fill = method, color = method),
              alpha = 0.35,
              type = "rect",
              extend = 1
            ),
            geom_cladelab(
              data = group_tree_annotation_df,
              mapping = aes(node = node, label = strain_id, color = method),
              hjust = 0,
              vjust = 2,
              geom = "text",
              angle = 0,
              barsize = 0,
              fontface = 2,
              size = 10
            )
          ),
          list(
            guides(color = "none"),
            labs(fill = "Method"),
            scale_color_manual(values = "#91a01e"),
            scale_fill_manual(values = "#91a01e",
                              guide = guide_legend(order = 2))
          )
        )
      }}
    
    
    # Create the SNP distance heatmap
    
    # Get the dataframe for making the heatmap 
    # In the heatmap, the snp distance between study genomes and the global genomes are not shown
    # because in the hypothesis we are not considering study and global genomes are not the same strain in hyper local transmission 
    
    group_tree_snp_matrix <- combined_snp_matrix %>%
      filter(
        subject %in% group_total_genomes & query %in% group_total_genomes
      )
    
    all_tip_order <- group_tree_newick$tip.label[group_tree_newick$edge[group_tree_newick$edge[,2] <= length(group_tree_newick$tip.label), 2]]
    study_tip_order <- all_tip_order[!grepl("GCA_|GCF_",all_tip_order)]
    
    
    heatmap_df <- do.call(rbind,
                          lapply(study_tip_order,
                                 function(genome_entry){
                                   
                                   data.frame(
                                     column_genome = genome_entry,
                                     row_genome = all_tip_order,
                                     gsnp = sapply(all_tip_order,
                                                   function(row_genome_entry){
                                                     if (row_genome_entry == genome_entry) {
                                                       0
                                                     } else if (grepl("GCA_|GCF_", row_genome_entry)) {
                                                       NA_real_
                                                     } else {
                                                       snp_value <- group_tree_snp_matrix$gsnp[
                                                         (group_tree_snp_matrix$subject == genome_entry &
                                                            group_tree_snp_matrix$query == row_genome_entry) |
                                                           (group_tree_snp_matrix$subject == row_genome_entry &
                                                              group_tree_snp_matrix$query == genome_entry)
                                                       ]
                                                       if (length(snp_value) == 0) NA_real_ else snp_value[1]
                                                     }
                                                   }),
                                     snp_quality = sapply(all_tip_order,
                                                          function(row_genome_entry){
                                                            if (row_genome_entry == genome_entry) {
                                                              "good"
                                                            } else if (grepl("GCA_|GCF_", row_genome_entry)) {
                                                              NA_character_
                                                            } else {
                                                              snp_quality <- group_tree_snp_matrix$snp_quality[
                                                                (group_tree_snp_matrix$subject == genome_entry &
                                                                   group_tree_snp_matrix$query == row_genome_entry) |
                                                                  (group_tree_snp_matrix$subject == row_genome_entry &
                                                                     group_tree_snp_matrix$query == genome_entry)
                                                              ]
                                                              if (length(snp_quality) == 0) NA_character_ else snp_quality[1]
                                                            }
                                                          })
                                   )
                                 }))
    
    
    heatmap_df <- heatmap_df %>%
      mutate(
        column_genome = factor(column_genome, levels = study_tip_order)
      )
    
    # Another SNP matrix containing study genomes compared to each of their top 3 global genomes
    global_heatmap_df <- do.call(rbind,
                                 lapply(study_tip_order,
                                        function(genome_entry){
                                          
                                          global_genome_entry <- grep("GCA_|GCF",group_tree_snp_matrix$query[
                                            group_tree_snp_matrix$subject == genome_entry
                                          ],
                                          value = TRUE)

                                          # If there is no global genome entry, return NULL
                                          if(length(global_genome_entry) == 0){
                                            return(NULL)
                                          }else{
                                            return(
                                              data.frame(
                                            column_genome = genome_entry,
                                            row_genome = global_genome_entry,
                                            gsnp = sapply(global_genome_entry,
                                                          function(row_genome_entry){
                                                            group_tree_snp_matrix$gsnp[
                                                              group_tree_snp_matrix$subject == genome_entry &
                                                                group_tree_snp_matrix$query == row_genome_entry
                                                            ]
                                                          }),
                                            snp_quality = sapply(global_genome_entry,
                                                                 function(row_genome_entry){
                                                                   group_tree_snp_matrix$snp_quality[
                                                                     group_tree_snp_matrix$subject == genome_entry &
                                                                       group_tree_snp_matrix$query == row_genome_entry
                                                                   ]
                                                                 })
                                            
                                          )
                                            )
                                          }
                                          
                                          
                                        }))
    
    global_heatmap_df <- global_heatmap_df %>%
      mutate(
        column_genome = factor(column_genome, levels = study_tip_order)
      )
    
    snp_heatmap <- ggplot(heatmap_df,
                          aes(x = column_genome,
                              y = row_genome,
                              fill = gsnp)) +
      geom_tile(color = "black", linewidth = 0.5) +
      # Text of study genomes compared to other study genomes
      geom_text(aes(label = ifelse(is.na(gsnp), "", as.character(gsnp)),
                    color = snp_quality),
                size = 3.5,
                fontface = "bold",
                show.legend = TRUE) +
      # Text of study genomes compared to their top global genomes (cladebreaker)
      geom_text(data = global_heatmap_df,
                aes(x = column_genome,
                    y = row_genome,
                    label = gsnp,
                    color = snp_quality),
                inherit.aes = FALSE,
                size = 3.5,
                fontface = "bold",
                show.legend = TRUE) +
      scale_fill_gradient(
        low = "white",
        high = "#005587",
        limits = c(min(heatmap_df$gsnp, na.rm = TRUE), max(heatmap_df$gsnp, na.rm = TRUE)),
        na.value = "grey",
        guide = guide_colorbar(order = 4)
      ) +
      scale_color_manual(
        values = c("good" = "green", "poor" = "red"),
        labels = c("good" = "Good", "poor" = "Poor"),
        na.value = "transparent",
        na.translate = FALSE,
        guide = guide_legend(order = 3)
      ) +
      labs(fill = "SNP Distance",
           color = "SNP Quality") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = "transparent"),
        panel.grid = element_blank()
      ) + 
      coord_fixed()
    
    
    tree_snp_matrix_plot <- snp_heatmap %>% insert_left(group_tree)
    
    
    # Export the plot as pdf file
    
    n_genomes <- length(group_tree_newick$tip.label)
    dims <- get_pdf_dimensions(n_genomes)
    
    pdf(file = file.path(plot_export_dir, paste0("Group", group_id, "_", endpoint_method,"_strain_composition.pdf")),
        width = dims$width,
        height = dims$height)
    print(tree_snp_matrix_plot)
    dev.off()
    
    plot_rds[[paste0("Group", group_id)]] <- tree_snp_matrix_plot
  }
  return(plot_rds)
}

# Import input from snakemake
# Strain compositions for four endpoint methods
plateau_strains_rds_path <- snakemake@input[["plateau_strains_rds"]]
peak_strains_rds_path <- snakemake@input[["peak_strains_rds"]]
discrepancy_strains_rds_path <- snakemake@input[["discrepancy_strains_rds"]]
global_strains_rds_path <- snakemake@input[["global_strains_rds"]]
# Paths to ML trees of different groups
group_tree_path <- snakemake@input[["iqtree_group_path"]]
# Paths to study and global SNP matrices
study_snp_matrix_path <- snakemake@input[["study_snp_matrix"]]
global_snp_matrix_path <- snakemake@input[["global_snp_matrix"]]
# Whether or not to perfrom cladebreaker
use_cladebreaker <- snakemake@params[["use_cladebreaker"]]
use_cladebreaker <- as.logical(use_cladebreaker)
# Output dir for strain composition plots
output_dir <- snakemake@params[["output_dir"]]

# Run the function for four endpoint methods
# Plateau
plateau_plots_rds <- get_strain_compositions_plot(
  endpoint_method = "plateau",
  plateau_strains_rds_path = plateau_strains_rds_path,
  peak_strains_rds_path = peak_strains_rds_path,
  discrepancy_strains_rds_path = discrepancy_strains_rds_path,
  global_strains_rds_path = global_strains_rds_path,
  study_snp_matrix_path = study_snp_matrix_path,
  global_snp_matrix_path = global_snp_matrix_path,
  group_tree_path = group_tree_path,
  use_cladebreaker = use_cladebreaker,
  output_dir = output_dir
)

saveRDS(plateau_plots_rds,
        file = file.path(output_dir,"plateau","plateau_strain_tree_snp.RDS"))

# Peak
peak_plots_rds <- get_strain_compositions_plot(
  endpoint_method = "peak",
  plateau_strains_rds_path = plateau_strains_rds_path,
  peak_strains_rds_path = peak_strains_rds_path,
  discrepancy_strains_rds_path = discrepancy_strains_rds_path,
  global_strains_rds_path = global_strains_rds_path,
  study_snp_matrix_path = study_snp_matrix_path,
  global_snp_matrix_path = global_snp_matrix_path,
  group_tree_path = group_tree_path,
  use_cladebreaker = use_cladebreaker,
  output_dir = output_dir
)

saveRDS(peak_plots_rds,
        file = file.path(output_dir,"peak","peak_strain_tree_snp.RDS"))

# Discrepancy
discrepancy_plots_rds <- get_strain_compositions_plot(
  endpoint_method = "discrepancy",
  plateau_strains_rds_path = plateau_strains_rds_path,
  peak_strains_rds_path = peak_strains_rds_path,
  discrepancy_strains_rds_path = discrepancy_strains_rds_path,
  global_strains_rds_path = global_strains_rds_path,
  study_snp_matrix_path = study_snp_matrix_path,
  global_snp_matrix_path = global_snp_matrix_path,
  group_tree_path = group_tree_path,
  use_cladebreaker = use_cladebreaker,
  output_dir = output_dir
)

saveRDS(discrepancy_plots_rds,
        file = file.path(output_dir,"discrepancy","discrepancy_strain_tree_snp.RDS"))

# Global
global_plots_rds <- get_strain_compositions_plot(
  endpoint_method = "global",
  plateau_strains_rds_path = plateau_strains_rds_path,
  peak_strains_rds_path = peak_strains_rds_path,
  discrepancy_strains_rds_path = discrepancy_strains_rds_path,
  global_strains_rds_path = global_strains_rds_path,
  study_snp_matrix_path = study_snp_matrix_path,
  global_snp_matrix_path = global_snp_matrix_path,
  group_tree_path = group_tree_path,
  use_cladebreaker = use_cladebreaker,
  output_dir = output_dir
)

saveRDS(global_plots_rds,
        file = file.path(output_dir,"global","global_strain_tree_snp.RDS"))
