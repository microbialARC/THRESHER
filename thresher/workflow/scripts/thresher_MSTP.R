# Libraries ----
library(ggplot2)
library(scales)
library(reshape2)
library(patchwork)

# Function to generate the MSTPs ----
thresher_MSTP <- function(output_path,
                          determine_strains_input,
                          plateau_strains,
                          peak_strains,
                          global_strains,
                          discrepancy_strains){
  #Universal variables/setting
  line_color <- c("#C72373",
                  "#86C2E8",
                  "#8F8578")
  names(line_color) <- c("Discrepancy",
                         "Clones",
                         "Singletons")
  # Create a list to store all plots
  all_combined_plots <- list()
  
  # Iterate groups with no less than 2 genomes 
  for(i in seq_along(determine_strains_input)){
    group_id <- determine_strains_input[[i]][[1]]$HC_group
    group_input <- determine_strains_input[[i]]
    ## The main plot showing the numbers of clones/singletons/discrepancy at SNP thresholds ----
    ### Data frame for main plot 
    main_plot_df <- melt(do.call(rbind, 
                                 lapply(group_input, function(row) 
                                   data.frame(threshold = row$cutoff, 
                                              Discrepancy = row$discrepancy, 
                                              Singletons = row$after_correction_singletons, 
                                              Clones = row$after_correction_clones))), 
                         id.vars = "threshold", 
                         measure.vars = c("Discrepancy", "Singletons", "Clones"), 
                         variable.name = "category", 
                         value.name = "value")
    
    # Plateau position
    plateau_pos <- plateau_strains$plateaus$plateau[plateau_strains$plateaus$group == group_id]
    plateau_length <- unique(plateau_strains$plateaus$plateau_length)
    # Peak position
    peak_pos <- peak_strains$peaks$peak[peak_strains$peaks$group == group_id]
    # Discrepancy position
    discrepancy_pos <- discrepancy_strains$discrepancy$discrepancy[discrepancy_strains$discrepancy$group == group_id]
    # Global position
    global_pos <- global_strains$global$global[global_strains$global$group == group_id]
    
    
    ### Main plot
    main_plot <- ggplot(main_plot_df,
                        aes(x = threshold,
                            y = value,
                            color = category)) +
      scale_color_manual(values = line_color,
                         name = "Number") +
      geom_line(alpha = 0.75,
                linewidth = 2) +
      # The plateau
      annotate("rect",
               xmin=plateau_pos,
               xmax=plateau_pos + plateau_length,
               ymin=-Inf,
               ymax=Inf,
               alpha=0.35,
               fill="#B3C16D") + 
      annotate(x=plateau_pos + plateau_length,
               y=+Inf,
               label=paste0("Plateau: ",plateau_pos," - ",plateau_pos+plateau_length),
               vjust=5,
               geom="label",
               color="#829726",
               fill="#E5EACC") + 
      # The peak 
      annotate(x=peak_pos,
               y=+Inf,
               label=paste0("Peak:", peak_pos),
               vjust=4,
               geom="label",
               color="#315784",
               fill="#DFF0FA") + 
      # The discrepancy 
      annotate(x=discrepancy_pos,
               y=+Inf,
               label=paste0("Discrepancy: ",discrepancy_pos),
               vjust=3,
               geom="label",
               color="#C72575",
               fill="#F0D1E2") + 
      # The global
      annotate(x=global_pos,
               y=+Inf,
               label=paste0("Global: ",global_pos),
               vjust=2,
               geom="label",
               color="#F2C100",
               fill="#F3EAAF") + 
      scale_x_continuous(name = "SNP Threshold",
                         breaks = seq(0, max(main_plot_df$threshold),50),
                         limits = c(min(main_plot_df$threshold),
                                    max(main_plot_df$threshold))) + 
      scale_y_continuous(
        name = "Number",
        breaks = pretty_breaks(),
      ) + 
      theme(
        axis.title.y.left = element_text(colour = "black",
                                         size = 30,
                                         face = "bold"),
        axis.text.y.left = element_text(colour = "black",
                                        size = 20),
        axis.title.x = element_text(color = "black",
                                    size = 30,
                                    face = "bold"),
        axis.text.x = element_text(colour = "black",
                                   size = 20),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.key.size = unit(1.5, 'cm'),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 2)
      )
    
    
    ## The secondary plot showing the mean/median bootstrap supports for strains ----
    
    secondary_plot_df <- melt(do.call(rbind, 
                                      lapply(group_input, function(row) 
                                        data.frame(threshold = row$cutoff, 
                                                   Mean = row$mean_strain_bootstrap_support,
                                                   Median = row$median_strain_bootstrap_support))), 
                              id.vars = "threshold", 
                              measure.vars = c("Mean", "Median"), 
                              variable.name = "Bootstrap", 
                              value.name = "value")
    
    
    secondary_plot <- ggplot(secondary_plot_df,
                             aes(x = threshold,
                                 y = Bootstrap,
                                 fill = value)) + 
      geom_tile() + 
      scale_fill_viridis_c(name = "Value") + 
      scale_x_continuous(name = "SNP Threshold",
                         breaks = seq(0, max(secondary_plot_df$threshold),50),
                         limits = c(min(secondary_plot_df$threshold),
                                    max(secondary_plot_df$threshold))) + 
      theme(
        axis.title.y.left = element_text(colour = "black",
                                         size = 20,
                                         face = "bold"),
        axis.text.y.left = element_text(colour = "black",
                                        size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.75, 'cm'),
        legend.position = "right",
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 2)
      ) + 
      coord_fixed(ratio=50)
    
    # Combine the main and secondary plots
    combined_plot <- secondary_plot + main_plot +
      plot_layout(ncol=1)
    
    # Store the combined plot in the list with group_id as name
    all_combined_plots[[paste0("Group", group_id)]] <- combined_plot
    
    # Export the combined plot
    pdf(file= file.path(output_path,
                        paste0("Group",
                               group_id,
                               "_MSTP.pdf")),
        width=15,
        height=10)
    print(combined_plot)
    dev.off()
  }
  return(all_combined_plots)
}

# Import from snakemake
output_path <- snakemake@params[["MSTP_dir"]]
system(paste0("mkdir -p ",output_path))
determine_strains_input <- readRDS(snakemake@input[["thresher_input"]])
peak_strains <- readRDS(snakemake@input[["peak_strains_rds"]])
plateau_strains <- readRDS(snakemake@input[["plateau_strains_rds"]])
global_strains <- readRDS(snakemake@input[["global_strains_rds"]])
discrepancy_strains <- readRDS(snakemake@input[["discrepancy_strains_rds"]])

all_combined_plots <- thresher_MSTP(output_path,
              determine_strains_input,
              plateau_strains,
              peak_strains,
              global_strains,
              discrepancy_strains)

saveRDS(all_combined_plots, 
        file = file.path(output_path,
                         "MSTP.RDS"))
