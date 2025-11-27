# Persistence plot: 
#x axis indicates the month
#y axis indicates the clusters
#circles indicates the genomes
# Libraries ----
library(ggplot2)
library(dplyr)
library(ggrepel)
# Import from snakemake ----
output_dir <- snakemake@params[["output_dir"]]
original_metadata_path <- snakemake@params[["original_metadata"]]
new_metadata_path <- snakemake@params[["new_metadata"]]
original_metadata <- read.csv(original_metadata_path,
                     header = FALSE,
                     sep = "\t") %>%
  select(V1,V4,V5)

new_metadata <- read.csv(new_metadata_path,
                         header = FALSE,
                         sep = "\t") %>%
  select(V1,V4,V5)

metadata <- rbind(original_metadata,
                  new_metadata)

colnames(metadata) <- c("genome","patientID","collection_date")
clusters_summary_path <- snakemake@input[["clusters_summary_rds"]]
clusters_summary_rds <- readRDS(clusters_summary_path)
clusters_summary <- clusters_summary_rds$clusters
output_pdf <- snakemake@output[["persistence_plot_pdf"]]
output_rds <- snakemake@output[["persistence_plot_rds"]]
# Function to generate plots  ----
persistence_plot <- function(output_dir,
                             clusters_summary,
                             metadata){
  setwd(dir = output_dir)
  # Only make the plot when length(clusters_summary) > 0
  if(length(clusters_summary) > 0){
    # Get the data frame as the input for the plot 
    all_plot_df <- do.call(rbind,
                           lapply(clusters_summary,
                                  function(cluster){
                                    
                                    data.frame(
                                      date = sapply(cluster$genomes,
                                                    function(genome) metadata$collection_date[metadata$genome == genome]),
                                      cluster = cluster$cluster,
                                      patient = sapply(cluster$genomes,
                                                       function(genome){
                                                         metadata$patientID[metadata$genome == genome]
                                                       })
                                    )
                                  }))
    # Get the data frame as the input for label
    all_label_df <- unique(all_plot_df)
    # Change the class 
    all_plot_df$date <- as.Date(all_plot_df$date)
    all_label_df$date <- as.Date(all_label_df$date)
    all_plot_df$cluster <- factor(all_plot_df$cluster,
                                  levels =  sort(unique(as.integer(all_plot_df$cluster)),decreasing = TRUE))
    
    all_label_df$cluster <- factor(all_label_df$cluster,
                                   levels =  sort(unique(as.integer(all_label_df$cluster)),decreasing = TRUE))
    
    all_persistence_plot <- ggplot(all_plot_df,
                                   aes(x = date,
                                       y = cluster)) +
      geom_line(linewidth = 1,
                aes(group = cluster),
                color = "black") +
      geom_point(size = 15,
                 color = "#011F5B",
                 alpha = 0.75) +
      scale_x_date(date_labels = "%b %Y") + 
      geom_label_repel(data = all_label_df,
                       aes(x = date,
                           y = cluster,
                           label = patient),
                       fill = "#DAEAF8",
                       color = "#25537D",
                       direction = "both",
                       label.size = 0.7,
                       nudge_y = 0.175,
                       show.legend = FALSE,
                       segment.alpha = 0.35,
                       max.overlaps = 1000) +
      labs(x = "Month",
           y = "Cluster") +
      theme(
        axis.line.x = element_line(linewidth = 2),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0.5,"cm"),
        axis.text.x = element_text(face = "bold",
                                   angle = 45,
                                   hjust = 1,
                                   size = 40),
        axis.text.y = element_text(face = "bold",
                                   size = 40),
        axis.title.x = element_text(face = "bold", size = 50),
        axis.title.y = element_text(face = "bold", size = 50),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
      )
    
    pdf(file=output_pdf,
        width=17.5,
        height=20)
    print(all_persistence_plot)
    dev.off()
    # Export the RDS of the plot
    saveRDS(all_persistence_plot,
            file = output_rds)
  }else{
    # Create an empty plot with text when no clusters are found
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = "No Transmission Cluster Found",
               size = 14, fontface = "bold") +
      theme_void() +
      theme(
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")
      )
    
    pdf(file=output_pdf,
        width=17.5,
        height=20)
    print(empty_plot)
    dev.off()
    # Export the RDS of the empty plot
    saveRDS(empty_plot,
            file = output_rds)
  }
}

# Execute the function
persistence_plot(output_dir = output_dir,
                 clusters_summary = clusters_summary,
                 metadata = metadata)
