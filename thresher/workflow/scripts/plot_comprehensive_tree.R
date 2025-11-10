# Funciton to visualize the comprehensive-tree annotated with HC groups and MLST ----
comprehensive_tree_visual <- function(comprehensive_tree_path,
                                      hc_groups_path,
                                      mlst_path,
                                      genome_species,
                                      output_dir,
                                      ncores){
  setwd(dir = output_dir)
  # Libraries
  library(ggplot2)
  library(ggnewscale)
  library(ape)
  library(phytools)
  library(colorspace)
  library(ggtree)
  library(ggtreeExtra)
  library(ggstar)
  library(TreeTools)
  library(parallel)
  
  # read the tree
  comprehensive_tree <- ape::read.tree(comprehensive_tree_path)
  # mid-point rooting of the tree
  comprehensive_tree <- phytools::midpoint_root(comprehensive_tree)
  # pre-order the tree
  comprehensive_tree <- Preorder(comprehensive_tree)
  # Read HC groups 
  hc_group_df <- read.csv(hc_groups_path,
                        sep = ",",
                        header = TRUE)
  # Read MLST
  mlst_df <- read.csv(mlst_path,sep = "\t",header = TRUE)
  # Find the node of the HC groups
  all_nodes <- (length(comprehensive_tree$tip.label) + 1):(length(comprehensive_tree$tip.label) + comprehensive_tree$Nnode)
  
  hc_groups_node <- do.call(rbind,
                            mclapply(sort(unique(hc_group_df$group)),
                                     function(group_id){
                                       
                                       group_genome <- hc_group_df$genome[hc_group_df$group == group_id]
                                       
                                       # Find the node for the group
                                       if(length(group_genome) > 1){
                                         
                                         group_node <- unlist(sapply(all_nodes,
                                                                     function(node){
                                                                       
                                                                       sub_comprehensive_tree <- Subtree(comprehensive_tree,node)
                                                                       
                                                                       if(identical(sort(sub_comprehensive_tree$tip.label),
                                                                                    sort(group_genome))){
                                                                         return(node)
                                                                       }
                                                                     }))
                                       }else{
                                         
                                         group_node <- which(comprehensive_tree$tip.label == group_genome)
                                         
                                       }
                                       
                                       # MLST for this group
                                       
                                       group_mlst <- unique(mlst_df$MLST[mlst_df$genome %in% group_genome])
                                       
                                       group_mlst <- if(length(group_mlst[group_mlst!= "Unassigned"])==1){
                                         
                                         group_mlst[group_mlst!= "Unassigned"]
                                         
                                       }else if(length(group_mlst[group_mlst!= "Unassigned"])==0){
                                         
                                         "Unassigned"
                                         
                                       }else if(length(group_mlst[group_mlst!= "Unassigned"])>1){
                                         
                                         names(which.max(table(mlst_df$MLST[mlst_df$genome %in% group_genome])))
                                       }
                                       
                                       return(data.frame(group = group_id,
                                                         node = group_node,
                                                         MLST = group_mlst))
                                     },
                                     mc.cores = ncores))

  
  hc_groups_node$group <- factor(as.character(hc_groups_node$group),
                                 levels = as.character(sort(as.integer(hc_groups_node$group))))
  # HC_groups color 
  hc_groups_color <- setNames(qualitative_hcl(length(unique(hc_groups_node$group)), palette = "Dark 3"),
                              as.character(sort(unique(hc_groups_node$group))))
  
  mlst_color <- setNames(c(qualitative_hcl(length(sort(unique(mlst_df$MLST[mlst_df$MLST!="Unassigned"]))), palette = "Dark 3"),
                           "#808080"),
                         c(sort(unique(mlst_df$MLST[mlst_df$MLST!="Unassigned"])),
                           "Unassigned"))
  
  genome_mlst_name <- switch(genome_species,
                             "sau" = "Clonal Complex",
                             "cdiff" = "MLST Clade",
                             "sepi" = "Sequence Type",
                             "kp" = "Sequence Type")
  
  mlst_df$MLST <- switch(genome_species,
                      "sau" = factor(mlst_df$MLST,
                                     levels = c(paste0("CC",
                                                       sort(as.integer(gsub("CC","",unique(mlst_df$MLST)[unique(mlst_df$MLST) != "Unassigned"])))),
                                                "Unassigned")),
                      "cdiff" = factor(ifelse(mlst_df$MLST != "Unassigned", paste0("Clade", mlst$MLST), mlst_df$MLST),
                                       levels = c(paste0("Clade",
                                                         sort(as.integer(unique(mlst_df$MLST)[unique(mlst_df$MLST) != "Unassigned"]))),
                                                  "Unassigned")),
                      "sepi" = factor(as.character(mlst_df$MLST),
                                      levels = c(as.character(sort(as.integer(unique(mlst_df$MLST)[unique(mlst_df$MLST) != "Unassigned"]))),
                                                 "Unassigned")),
                      "kp" = factor(as.character(mlst_df$MLST),
                                    levels = c(as.character(sort(as.integer(unique(mlst_df$MLST)[unique(mlst_df$MLST) != "Unassigned"]))),
                                               "Unassigned")))
  
  # Visualize the tree with ggtree and ggtreeextra
  
  circular_tree <- ggtree(comprehensive_tree,
                          aes(color = as.numeric(label),
                              subset = !isTip & !is.na(as.numeric(label))),
                          layout = "fan",
                          size=0.5,
                          open.angle=15) +
    scale_color_continuous(name = "Bootstrap Support",
                           low = "#9E1A1A",
                           high = "#2E9B57") + 
    theme(
      legend.background = element_rect(fill = NA),
      legend.position = c(0, 0.5),
      legend.direction = "vertical"  
    )
  
  
  circular_tree_group <- circular_tree +
    ggnewscale::new_scale_color() + 
    geom_cladelab(
      data = hc_groups_node,
      mapping = aes(
        node = node,
        label = group,
        color = group
      ),
      geom = "text",
      angle = "auto",
      barsize = 2) + 
    geom_treescale(x = 0,
                   y = 0) +
    scale_fill_manual(values = hc_groups_color,
                      guide = "none") +
    scale_color_manual(values = hc_groups_color,
                       guide = "none") +
    new_scale_fill() + 
    theme(
      legend.background = element_rect(fill = NA),
      legend.position = "inside",
      legend.direction = "vertical"
    )
  
  circular_tree_MLST <- circular_tree +
    ggnewscale::new_scale_color() + 
    geom_cladelab(
      data = hc_groups_node,
      mapping = aes(
        node = node,
        label = MLST,
        color = MLST
      ),
      geom = "text",
      angle = "auto",
      barsize = 2) + 
    geom_treescale(x = 0,
                   y = 0) +
    scale_fill_manual(values = mlst_color,
                      guide = "none") +
    scale_color_manual(values = mlst_color,
                       guide = "none") +
    new_scale_fill() + 
    theme(
      legend.background = element_rect(fill = NA),
      legend.position = "inside",
      legend.direction = "vertical"
    )
    
  # Export circular tree annotated with MLST
  pdf(file=snakemake@output[["comprehensive_tree_mlst_pdf"]],
      width=7.5,
      height=5)
  print(circular_tree_MLST)
  dev.off()
  
  saveRDS(circular_tree_MLST,
          snakemake@output[["comprehensive_tree_mlst_rds"]])
 
  # Export circular tree annotated with HC group
  pdf(file=snakemake@output[["comprehensive_tree_group_pdf"]],
      width=7.5,
      height=5)
  print(circular_tree_group)
  dev.off()
  
  saveRDS(circular_tree_group,
          snakemake@output[["comprehensive_tree_group_rds"]])
  
}

# Import from Snakemake
comprehensive_tree_path <- snakemake@input[["contree"]]
hc_groups_path <- snakemake@input[["hc_groups"]]
mlst_path <- snakemake@input[["mlst_results"]]
genome_species = snakemake@params[["genome_species"]]
output_dir <- snakemake@params[["output_dir"]]
ncores <- snakemake@threads
# Create output directory if not exists
if(!dir.exists(output_dir)){
  dir.create(output_dir,recursive = TRUE)
  }
# Execute the function
comprehensive_tree_visual(comprehensive_tree_path,
                          hc_groups_path,
                          mlst_path,
                          genome_species,
                          output_dir,
                          ncores)
