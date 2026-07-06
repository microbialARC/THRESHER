# Function to visualize the core-gene tree annotated with HC groups and MLST ----
core_gene_tree_visual <- function(core_gene_tree_path,
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
  library(ggstar)
  library(TreeTools)
  library(parallel)
  
  # read the tree
  core_gene_tree <- ape::read.tree(core_gene_tree_path)
  # mid-point rooting of the tree
  core_gene_tree <- phytools::midpoint_root(core_gene_tree)
  # pre-order the tree
  core_gene_tree <- Preorder(core_gene_tree)
  # Read HC groups 
  hc_group_df <- read.csv(hc_groups_path,
                        sep = ",",
                        header = TRUE)
  # Read MLST
  mlst_df <- read.csv(mlst_path,sep = "\t",header = TRUE)
  # This should not happen; temporary patch to avoid glitched MLST clade
  # visualization for C. difficile (strips "nan"/".0" clade-label artifacts).
  if (genome_species == "cdiff") {
    mlst_df$MLST <- gsub(".0", "", gsub("Cladenan", "Unassigned", mlst_df$MLST), fixed = TRUE)
  }
  # Find the node of the HC groups
  all_nodes <- (length(core_gene_tree$tip.label) + 1):(length(core_gene_tree$tip.label) + core_gene_tree$Nnode)
  
  hc_groups_node <- do.call(rbind,
                            mclapply(sort(unique(hc_group_df$group)),
                                     function(group_id){
                                       
                                       group_genome <- hc_group_df$genome[hc_group_df$group == group_id]
                                       
                                       # Find the node for the group
                                       if(length(group_genome) > 1){
                                         
                                         group_node <- unlist(sapply(all_nodes,
                                                                     function(node){
                                                                       
                                                                       sub_core_gene_tree <- Subtree(core_gene_tree,node)
                                                                       
                                                                       if(identical(sort(sub_core_gene_tree$tip.label),
                                                                                    sort(group_genome))){
                                                                         return(node)
                                                                       }
                                                                     }))
                                       }else{
                                         
                                         group_node <- which(core_gene_tree$tip.label == group_genome)
                                         
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
  
  mlst_prefix <- switch(genome_species,
                        "sau" = "CC",
                        "cdiff" = "Clade",
                        "sepi" = "ST",
                        "kp" = "ST")
  
  mlst_levels <- unique(mlst_df$MLST[mlst_df$MLST != "Unassigned"])
  mlst_levels <- paste0(mlst_prefix,
                        sort(as.integer(gsub(mlst_prefix, "", mlst_levels))))
  
  mlst_df$MLST <- factor(as.character(mlst_df$MLST),
                         levels = c(mlst_levels, "Unassigned"))

  # Size-adaptive parameters based on number of tips
  n_tips <- length(core_gene_tree$tip.label)
  # PDF dimensions: scale with sqrt(n_tips) so size grows but doesn't explode.
  pdf_size <- max(7.5, min(40, 7.5 * sqrt(n_tips / 100)))
  pdf_width <- pdf_size
  pdf_height <- pdf_size * (5 / 7.5)
  # Branch linewidth: thinner for dense trees so branches don't merge visually.
  tree_linewidth <- max(0.15, min(0.6, 0.5 * (100 / n_tips)^0.5))
  # Clade label / bar size: scales similarly to linewidth.
  clade_barsize  <- max(0.6, min(2.5, 2 * (100 / n_tips)^0.4))
  # Font size for clade labels (group / MLST). Smaller as tips increase.
  clade_fontsize <- max(1.5, min(4, 3.88 * (100 / n_tips)^0.4))
  # Legend / theme text size
  base_fontsize  <- max(7, min(11, 10 * (100 / n_tips)^0.15))

  # Visualize the tree with ggtree
  
  circular_tree <- ggtree(core_gene_tree,
                          aes(color = as.numeric(label),
                              subset = !isTip & !is.na(as.numeric(label))),
                          layout = "fan",
                          linewidth=tree_linewidth,
                          open.angle=15) +
    scale_color_continuous(name = "Bootstrap Support",
                           low = "#9E1A1A",
                           high = "#2E9B57") + 
    theme(
      legend.background = element_rect(fill = NA),
      legend.position = c(0, 0.5),
      legend.direction = "vertical",
      legend.text  = element_text(size = base_fontsize),
      legend.title = element_text(size = base_fontsize)
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
      barsize = clade_barsize,
      fontsize = clade_fontsize) + 
    geom_treescale(x = 0,
                   y = 0,
                   fontsize = clade_fontsize) +
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
      barsize = clade_barsize,
      fontsize = clade_fontsize) + 
    geom_treescale(x = 0,
                   y = 0,
                   fontsize = clade_fontsize) +
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
  pdf(file=snakemake@output[["core_gene_tree_mlst_pdf"]],
      width=pdf_width,
      height=pdf_height)
  print(circular_tree_MLST)
  dev.off()
  
  saveRDS(circular_tree_MLST,
          snakemake@output[["core_gene_tree_mlst_rds"]])
 
  # Export circular tree annotated with HC group
  pdf(file=snakemake@output[["core_gene_tree_group_pdf"]],
      width=pdf_width,
      height=pdf_height)
  print(circular_tree_group)
  dev.off()
  
  saveRDS(circular_tree_group,
          snakemake@output[["core_gene_tree_group_rds"]])
  
}

# Import from Snakemake
core_gene_tree_path <- snakemake@input[["contree"]]
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
core_gene_tree_visual(core_gene_tree_path,
                          hc_groups_path,
                          mlst_path,
                          genome_species,
                          output_dir,
                          ncores)
