# Libraries ----
library(TreeTools)
library(ape)
library(ggplot2)
library(ggExtra)
library(dplyr)
library(showtext)
library(ggrepel)
library(parallel)
# Function to calculate the average snp / phylogenetic distance and plot----

thresher_qc <- function(method_strain_path,
                        comprehensive_tree_path,
                        snp_matrix_path) {
  # Method Name 
  method_name <- gsub("\\_strains.RDS","",basename(method_strain_path))
  # Read the strains summary 
  method_strain <- readRDS(method_strain_path)
  # Read comprehensive tree and generate the distance matrix
  comprehensive_tree <- read.tree(comprehensive_tree_path)
  comprehensive_tree <- phytools::midpoint_root(comprehensive_tree)
  comprehensive_tree <- Preorder(comprehensive_tree)
  tree_distance_matrix <-  cophenetic.phylo(comprehensive_tree)
  # Read the study snp matrix 
  snp_matrix <- readRDS(snp_matrix_path)
  #All strains identified with the method
  all_strain_id <- sort(unique(method_strain$strains$strain_id))
  #unique comparisons among clones in the group
  if(length(all_strain_id) > 1){
    # Only proceed if there are more than one strain identified
    strain_comparisons <- data.frame(t(combn(all_strain_id, 2))) %>%
      setNames(c("subject", "query")) %>%
      mutate(
        # label if the comparison is in the same group
        same_group = sapply(strsplit(subject, "_"), `[`, 1) == sapply(strsplit(query, "_"), `[`, 1)
      ) 
    
    # Parallel 
    # Setup parallel processing
    cl <- makeCluster(detectCores())
    clusterExport(cl,
                  c("strain_comparisons",
                    "snp_matrix",
                    "tree_distance_matrix",
                    "method_strain"),
                  envir = environment())
    
    strain_comparison_results <- do.call(rbind,
                                         parLapplyLB(cl,
                                                     seq_len(nrow(strain_comparisons)),
                                                     function(row_entry){
                                                       subject_strain_genome <- method_strain$strains$genome[method_strain$strains$strain_id == strain_comparisons$subject[row_entry]]
                                                       query_strain_genome <- method_strain$strains$genome[method_strain$strains$strain_id == strain_comparisons$query[row_entry]]
                                                       
                                                       strain_comparisons_tmp <- data.frame(subject = subject_strain_genome,
                                                                                            query = rep(query_strain_genome,length(subject_strain_genome)))
                                                       
                                                       strain_comparisons_tmp$gsnp <- sapply(seq_len(nrow(strain_comparisons_tmp)),
                                                                                             function(secondary_row){
                                                                                               snp_matrix$gsnp[snp_matrix$subject %in% c(strain_comparisons_tmp$subject[secondary_row],
                                                                                                                                         strain_comparisons_tmp$query[secondary_row]) &
                                                                                                                 snp_matrix$query %in% c(strain_comparisons_tmp$subject[secondary_row],
                                                                                                                                         strain_comparisons_tmp$query[secondary_row])]
                                                                                             })
                                                       
                                                       strain_comparisons_tmp$phylogeny_distance <- sapply(seq_len(nrow(strain_comparisons_tmp)),
                                                                                                           function(secondary_row){
                                                                                                             tree_distance_matrix[strain_comparisons_tmp$subject[secondary_row],
                                                                                                                                  strain_comparisons_tmp$query[secondary_row]]
                                                                                                           })
                                                       
                                                       
                                                       
                                                       return(data.frame(
                                                         subject = strain_comparisons$subject[row_entry],
                                                         query = strain_comparisons$query[row_entry],
                                                         snp_average_distance = sum(strain_comparisons_tmp$gsnp) / (length(subject_strain_genome) * length(query_strain_genome)),
                                                         phylogeny_average_distance = sum(strain_comparisons_tmp$phylogeny_distance) / (length(subject_strain_genome) * length(query_strain_genome)),
                                                         same_group = strain_comparisons$same_group[row_entry]
                                                       ))
                                                       
                                                     }))
    stopCluster(cl)
    rm(cl)
    #export csv file
    write.csv(strain_comparison_results,
              file = paste0(method_name,"_qc_table.csv"),
              quote = FALSE,
              row.names = FALSE)
    
    
    # Plot 
    strain_comparison_results$same_group <- factor(strain_comparison_results$same_group,
                                                   levels = c(TRUE,FALSE))
    
    distance_plot <- ggplot(strain_comparison_results,
                            aes(x = phylogeny_average_distance,
                                y = snp_average_distance)) + 
      geom_point(aes(color = same_group),
                 size = 2.5,
                 alpha = 0.75) + 
      scale_color_manual(name = "Same Hierarchical Clustering Group",
                         values = setNames(c("#67A5D2",
                                             "#B22069"),
                                           c(TRUE,FALSE))) +
      geom_text_repel(
        data = subset(strain_comparison_results,
                      snp_average_distance <= 100),
        aes(label = paste0(subject,
                           " ",
                           query,
                           ":",
                           round(snp_average_distance,1)),
            color = same_group),
        size = 5,
        alpha = 1,
        box.padding = 0.5,
        point.padding = 0.5,
        force = 2,
        show.legend = FALSE
      ) + 
      labs(x = "Phylogenetic Average Distance",
           y = "gSNP Average Distance") +
      geom_hline(yintercept = 100,
                 color = "#829724",
                 linetype='dotted',
                 linewidth = 1) + 
      theme(axis.text.x = element_text(face = "bold",
                                       size = 20),
            axis.text.y = element_text(face = "bold",
                                       size = 20),
            axis.title.x = element_text(face = "bold",
                                        size = 25),
            axis.title.y = element_text(face = "bold",
                                        size = 25),
            plot.background = element_rect(fill = "transparent"),
            panel.background = element_rect(fill = "transparent"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = "transparent",
                                      fill = "transparent"),
            legend.text = element_text(size = 12.5),
            legend.title = element_text(size = 17.5),
            axis.line = element_line(),
            axis.ticks.length=unit(0.25,"cm"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
            legend.position="bottom"
      ) + 
      scale_x_continuous() + 
      scale_y_continuous()
    
    distance_plot <- ggMarginal(distance_plot,
                                type = "density",
                                groupFill = TRUE,
                                groupColour = TRUE)
    
    pdf(file=paste0(method_name,"_qc_plot.pdf"),
        width=10,
        height=10)
    print(distance_plot)
    dev.off()
  }else{
    # If there is only one strain identified, output a message in the output file
    one_strain_message <- paste0("Only one strain identified by the method ",method_name,". QC plot and table are not generated.")
    pdf(file=paste0(method_name,"_qc_plot.pdf"),
        width=10,
        height=1)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    text(0.5, 0.5, one_strain_message, cex = 1.2)
    dev.off()

    write.csv(data.frame(message = one_strain_message),
              file = paste0(method_name,"_qc_table.csv"),
              quote = FALSE,
              row.names = FALSE)
  }
  
}

# Import from snakemake ----
output_path <- snakemake@params[["QC_dir"]]
system(paste0("mkdir -p ",output_path))
setwd(dir = output_path)

comprehensive_tree_path <- snakemake@input[["comprehensive_tree"]]
snp_matrix_path <- snakemake@input[["study_snp_matrix"]]

#method strains 
peak_strains_path <- snakemake@input[["peak_strains_rds"]]
plateau_strains_path <- snakemake@input[["plateau_strains_rds"]]
global_strains_path <- snakemake@input[["global_strains_rds"]]
discrepancy_strains_path <- snakemake@input[["discrepancy_strains_rds"]]

# run the function ----
#plateau
thresher_qc(plateau_strains_path,
            comprehensive_tree_path,
            snp_matrix_path)
#peak
thresher_qc(peak_strains_path,
            comprehensive_tree_path,
            snp_matrix_path)
#global
thresher_qc(global_strains_path,
            comprehensive_tree_path,
            snp_matrix_path)
#discrepancy
thresher_qc(discrepancy_strains_path,
            comprehensive_tree_path,
            snp_matrix_path)

# remove the Rplots.pdf file generated by ggplot2
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}