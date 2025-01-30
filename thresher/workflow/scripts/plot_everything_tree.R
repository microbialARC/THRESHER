# Funciton to visualize the everything-tree annotated with HC groups and MLST ----
everything_tree_visual <- function(everything_tree_path,
                                   hc_groups_path,
                                   mlst_path,
                                   output_dir){
  setwd(dir = output_dir)
  # Libraries
  library(ggplot2)
  library(ape)
  library(phytools)
  library(colorspace)
  library(dplyr)
  library(ggtree)
  library(TreeTools)
  
  # read the tree
  everything_tree <- ape::read.tree(everything_tree_path)
  # mid-point rooting of the tree
  everything_tree <- phytools::midpoint_root(everything_tree)
  # pre-order the tree
  everything_tree <- Preorder(everything_tree)
  # Read HC groups 
  hc_groups <- read.csv(hc_groups_path,
                        sep = ",",
                        header = TRUE)
  # Find the node of the HC groups
  all_nodes <- (length(everything_tree$tip.label) + 1):(length(everything_tree$tip.label) + everything_tree$Nnode)
  hc_groups_node <- do.call(rbind,
                            lapply(all_nodes,
                                   function(node){
                                     sub_everything_tree <- Subtree(everything_tree,node)
                                     for(group in sort(unique(hc_groups$group))){
                                       if(identical(sort(sub_everything_tree$tip.label),sort(hc_groups$genome[hc_groups$group == group]))){
                                         
                                           return(data.frame(group = group,
                                                      node = node))
                                       }
                                     }
                                   }
                            ))
  
  hc_groups_node$group <- factor(as.character(hc_groups_node$group),
                                 levels = as.character(sort(as.integer(hc_groups_node$group))))
  # HC_groups color 
  hc_groups_color <- setNames(qualitative_hcl(length(unique(hc_groups$group)), palette = "Dark 3"),
                              as.character(sort(unique(hc_groups$group))))
  # Read MLST
  mlst <- read.csv(mlst_path,sep = "\t",header = TRUE)
  mlst_color <- setNames(c(qualitative_hcl(length(sort(unique(mlst$MLST[mlst$MLST!="Unassigned"]))), palette = "Dark 3"),
                           "#808080"),
                         c(sort(unique(mlst$MLST[mlst$MLST!="Unassigned"])),
                           "Unassigned"))
  # Visualize the tree with ggtree and ggtreeextra
  circular_tree <- ggtree(everything_tree,
                                layout = "fan",
                                size=1) + 
    geom_hilight(
      data = hc_groups_node,
      aes(node = node,
          fill = group),
      alpha = 0.35,
      type = "rect"
    ) +
  scale_fill_manual(name = "Group",
                    values = hc_groups_color)
  
  
  circular_tree <- circular_tree %<+% mlst + geom_tippoint(aes(color=MLST),
                                                           size = 3.5,
                                                           alpha = 0.8) +
    scale_color_manual(name = "MLST",
                       values = mlst_color)
  
  pdf(file=file.path(output_dir,"everything_tree.pdf"),
      width=7.5,
      height=5)
  print(circular_tree)
  dev.off()
}

# Import from snakemake
everything_tree_path <- snakemake@input[["contree"]]
hc_groups_path <- snakemake@input[["hc_groups"]]
mlst_path <- snakemake@input[["mlst_results"]]
output_dir <- snakemake@params[["output_dir"]]
# Execute the function
everything_tree_visual(everything_tree_path,
                       hc_groups_path,
                       mlst_path,
                       output_dir)
