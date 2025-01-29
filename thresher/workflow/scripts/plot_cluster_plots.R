# Libraries ----
library(ggplot2)
library(dplyr)
#geom_circle
library(ggforce)
library(shadowtext)
library(purrr)
# Functions exported from ggchicklet ----
ggname <- function(prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}
# Copy the geom_rrect() function from the ggchicklet(https://github.com/hrbrmstr/ggchicklet/blob/master/R/a-geom-rect.R)
# Because it's not in the conda env and I need to use snakemake

geom_rrect <- function(mapping = NULL, data = NULL, # nocov start
                       stat = "identity", position = "identity",
                       radius = grid::unit(6, "pt"),
                       ...,
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRrect,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      radius = radius,
      na.rm = na.rm,
      ...
    )
  )
}

GeomRrect <- ggplot2::ggproto(
  "GeomRrect", ggplot2::Geom,
  
  default_aes = ggplot2::aes(
    colour = NA, fill = "grey35", size = 0.5, linetype = 1, alpha = NA
  ),
  
  required_aes = c("xmin", "xmax", "ymin", "ymax"),
  
  draw_panel = function(self, data, panel_params, coord,
                        radius = grid::unit(6, "pt")) {
    
    coords <- coord$transform(data, panel_params)
    
    lapply(1:length(coords$xmin), function(i) {
      
      grid::roundrectGrob(
        coords$xmin[i], coords$ymax[i],
        width = (coords$xmax[i] - coords$xmin[i]),
        height = (coords$ymax[i] - coords$ymin)[i],
        r = radius,
        default.units = "native",
        just = c("left", "top"),
        gp = grid::gpar(
          col = coords$colour[i],
          fill = alpha(coords$fill[i], coords$alpha[i]),
          lwd = coords$size[i] * .pt,
          lty = coords$linetype[i],
          lineend = "butt"
        )
      )
      
    }) -> gl
    
    grobs <- do.call(grid::gList, gl)
    
    ggname("geom_rrect", grid::grobTree(children = grobs))
    
  },
  
  draw_key = ggplot2::draw_key_polygon
  
) # nocov end



# Import from snakemake ----
output_dir <- snakemake@params[["output_dir"]]
system(paste0("mkdir -p ",output_dir))
endpoint <- snakemake@params[["endpoint"]]
species <- snakemake@params[["species"]]
metadata_path <- snakemake@params[["metadata"]]
mlst_results_path <- snakemake@input[["mlst_results"]]
peak_strains_rds_path <- snakemake@input[["peak_strains_rds"]]
plateau_strains_rds_path <- snakemake@input[["plateau_strains_rds"]]
discrepancy_strains_rds_path <- snakemake@input[["discrepancy_strains_rds"]]
global_strains_rds_path <- snakemake@input[["global_strains_rds"]]

strains_rds <- switch(
  endpoint,
  "plateau" = readRDS(plateau_strains_rds_path),
  "peak" = readRDS(peak_strains_rds_path),
  "global" = readRDS(global_strains_rds_path),
  "discrepancy" = readRDS(discrepancy_strains_rds_path)
)
# read metadata 
metadata <- read.csv(metadata_path,
                     header = FALSE,
                     sep = "\t") %>%
  select(V1,V4,V5)
# read mlst 
mlst_results <- read.csv(mlst_results_path,
                         sep = "\t")

if(species == "sau"){
  mrsa_results <- read.csv(snakemake@input[["mrsa_results"]])
}
# Summarize the clusters based on the endpoint method ----
clusters_summary <- lapply(unique(strains_rds[["strains"]]$strain_id), function(strain) {
  strain_genomes <- strains_rds[["strains"]]$genome[strains_rds[["strains"]]$strain_id == strain]
  strain_patients <- unique(metadata$V4[metadata$V1 %in% strain_genomes])
  
  if(length(strain_patients) > 1) {
    list(
      cluster = NA,
      strain = strain,
      genomes = strain_genomes,
      MLST = unique(mlst_results[[3]][mlst_results$genome %in% strain_genomes]),
      AMR = if(species == "sau") mrsa_results$MRSA[mrsa_results$strain == strain] else "N/A",
      patients = strain_patients,
      first_seen = format(min(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])), "%y-%m-%d"),
      last_seen = format(max(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])), "%y-%m-%d"),
      persistence = as.integer(max(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])) - min(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])))
    )
  }
}) %>%
  compact() %>%
  {.[order(sapply(., function(x) length(x$patients)), decreasing = TRUE)]} %>%
  mapply(function(x, i) {x$cluster <- i; x}, ., seq_along(.), SIMPLIFY = FALSE)

                 

if(length(clusters_summary) > 0){
  clusters_summary_csv <- do.call(rbind, lapply(unique(strains_rds[["strains"]]$strain_id), function(strain) {
    strain_genomes <- strains_rds[["strains"]]$genome[strains_rds[["strains"]]$strain_id == strain]
    strain_patients <- unique(metadata$V4[metadata$V1 %in% strain_genomes])
    
    if(length(strain_patients) > 1) {
      data.frame(
        cluster = NA,
        strain = strain,
        MLST = paste(unique(mlst_results[[3]][mlst_results$genome %in% strain_genomes]),collapse = "|"),
        AMR = if(species == "sau") mrsa_results$MRSA[mrsa_results$strain == strain] else "N/A",
        genomes = paste(strain_genomes, collapse = "|"),
        patients = paste(strain_patients, collapse = "|"),
        first_seen = format(min(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])), "%y-%m-%d"),
        last_seen = format(max(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])), "%y-%m-%d"),
        persistence = as.integer(max(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])) - min(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])))
      )
    }
  })) %>%
    `[`(order(nchar(.$patients)), ) %>%
    within(cluster <- seq_along(cluster))
  write.csv(clusters_summary_csv,
            row.names = FALSE,
            file = snakemake@output[["clusters_summary_csv"]],
            quote = FALSE)
}else{
  clusters_summary_csv <- "No Clusters Found"
  write.csv(clusters_summary_csv,
            file = snakemake@output[["clusters_summary_csv"]],
            row.names = FALSE,
            quote = FALSE)
}


saveRDS(clusters_summary,
        snakemake@output[["clusters_summary_rds"]])



# Function to draw the cluster plot ----
clusterplot <- function(output_dir,
                        clusters_summary){
  
  setwd(dir = output_dir)
  
  #Size for export 
  pdf_export_size <- data.frame(num_circle = 2:100) %>%
    mutate(
      width = 5,
      height = case_when(
        num_circle <= 5 ~ 3.75,
        TRUE ~ ((ceiling(num_circle/5)*2.5*4 + 19 + 15)/44) * 3.75
      )
    )
  
  
  ## Iterate and make the plot for each cluster ----
  # Initialize an empty list to store results
  cluster_plots_rds <- list()
  
  for(i in seq_along(clusters_summary)){
    
    plot_members <- sort(as.integer(clusters_summary[[i]]$patients),
                         decreasing = TRUE)
    
    plot_members_num <- length(plot_members)
    
    # The data frame as ggplot input
    # There are 2 situations when making the plot:
    # 1: no more than 5 patients(1 line in the plot)
    # 2: more than 5 patients (more than 1 line in the plot)
    ### The plot_df for ggplot input ----
    if(plot_members_num <= 5){
      #situation1:
      #less than 5 patients in the cluster
      #there will always be only 1 line in the plot
      #so y of every patient is always 0.4 
      #the x of each patient will be set proportionally
      if(plot_members_num == 5){
        #If there are exactly 5 patients in the cluster
        #only 1 line
        #the x of each patient will be set proportionally(from n=1 to n=5, x = 0.1,0.3,0.5,0.7,and 0.9 respectively)
        #the y will be 0.5
        plot_df <- do.call(rbind,
                           lapply(seq_along(plot_members),
                                  function(entry){
                                    data.frame(
                                      ID = plot_members[entry],
                                      x = 1-((2*entry-1)/(2*5)),
                                      y = 0.5
                                    )
                                  }))
        
      }else{
        #less than 5 patients in the cluster
        #there will always be only 1 line in the plot
        #so y of every patient is always 0.4 
        #the x of each patient will be set proportionally
        
        plot_df <- do.call(rbind,
                           lapply(seq_along(plot_members),
                                  function(entry){
                                    data.frame(
                                      ID = plot_members[entry],
                                      x = 1-((2*(entry %% 5)-1)/(2*(plot_members_num %% 5))),
                                      y = 0.5
                                    )
                                  }))
        
      }
    }else{
      #situation2:
      #more 5 patients in the cluster
      #there will be more than 1 line in the plot
      #the x of each patient will be set proportionally
      #the y will vary depending on which line n is at
      #ceiling(nrow(plot_df)/5) is the total number of lines 
      #ceiling(n/5) indicates which line "n" is at
      #n %% 5 == 0 indicates this is the last one in this line
      plot_df <- do.call(rbind,
                         lapply(seq_along(plot_members),
                                function(entry){
                                  
                                  entry_id <- plot_members[entry]
                                  
                                  if(ceiling(entry/5) < ceiling(plot_members_num/5)){
                                    #ceiling(n/5) < ceiling(nrow(plot_df)/5
                                    #this line IS NOT the last line in the plot
                                    if(entry %% 5 == 0){
                                      data.frame(
                                        ID = entry_id,
                                        x = 0.1,
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                    }else{
                                      data.frame(
                                        ID = entry_id,
                                        x = 1-((2*(entry %% 5)-1)/10),
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                    }
                                  }else if(ceiling(entry/5) == ceiling(plot_members_num/5) & plot_members_num %% 5 == 0){
                                    
                                    #the LAST line in the plot
                                    #and also if nrow(plot_df) %% 5 == 0, which means the patient number IS 10,15,20,25, etc. 
                                    if(entry %% 5 == 0){
                                      data.frame(
                                        ID = entry_id,
                                        x = 0.1,
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                      
                                    }else{
                                      data.frame(
                                        ID = entry_id,
                                        x = 1-((2*(entry %% 5)-1)/(2*5)),
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                    }
                                  }else if(ceiling(entry/5) == ceiling(plot_members_num/5) & plot_members_num %% 5 != 0){
                                    #the LAST line in the plot
                                    #and also if nrow(plot_df) %% 5 != 0, which means the patient number IS NOT 10,15,20,25, etc. 
                                    if(entry %% 5 == 0){
                                      
                                      data.frame(
                                        ID = entry_id,
                                        x = 0.1,
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                      
                                    }else{
                                      data.frame(
                                        ID = entry_id,
                                        x = 1-((2*(entry %% 5)-1)/(2*(plot_members_num %% 5))),
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                    }
                                  }
                                }))
    }
    
    
    #now use plot_df to make the plot using ggplot2
    #In the rectangle of the plot, I want each line to have at most 5 dots
    #otherwise it will be too crowded
    #Assume the length of the rectangle is L
    #the number of the circles is N, the r of the circle is R = 2.5
    #the position of the circle center at X place is:  L(2X-1/2N)
    #I set L(the length of the rectangle) based on the number of the circles: L = 2*N*R (N < 5)
    #if N >= 5, L =  10R = 25
    #H(height of the rectangle) = ceiling(N/5)*R*4
    H <- ceiling(plot_members_num/5)*2.5*4
    
    L <- 2*10*2.5
  
    ### generate the plot ----
    #MLST of the cluster 
    cluster_mlst <- clusters_summary[[i]]$MLST[1][clusters_summary[[i]]$MLST[1] != "Unassigned"]
    cluster_amr <- clusters_summary[[i]]$AMR
    
    cluster_plot <- ggplot() +
      geom_rrect(aes(xmin = 0,
                     xmax = L,
                     ymin = 0,
                     ymax = H+15),
                 fill = "#DAEAF8",
                 color = "#25537D",
                 size = 3,
                 alpha = 1,
                 r = unit(0.1, 'npc')) +
      # upper bar showing AMR, Cluster ID, and MLST
      geom_rrect(aes(xmin = 2.5,
                     xmax = L-2.5,
                     ymin = H+2.5,
                     ymax = H+15-2.5),
                 fill = "#F2FAFF",
                 color = "transparent",
                 size = 0,
                 alpha = 1,
                 r = unit(0.175, 'npc')) + 
      geom_circle(data = plot_df,
                  aes(x0 = x*L,
                      y0 = y*H+2.5,
                      r = 4),
                  color = "transparent",
                  fill = "#011F5B",
                  linewidth = 1) + 
      guides(fill = "none") + 
      geom_text(data = plot_df,
                aes(x = x*L,
                    y = y*H+2.5,
                    label = ID),
                size = 10,
                color = "white",
                fontface = "bold") + 
      #Cluster ID
      geom_text(aes(x = L/2,
                    y = H+10,
                    label = paste0("Cluster",
                                   clusters_summary[[i]]$cluster)),
                size = 10.5,
                color = "black",
                fontface = "bold") +
      #MLST
      geom_text(aes(x = L/2,
                    y = H+2.5+3,
                    label = clusters_summary[[i]]$MLST[1]),
                size = 7.5,
                color = "black") +
      ggplot2::coord_fixed(ratio = 1,
                           xlim = c(0, L),
                           ylim = c(0, H+15)) +
      #MRSA or MSSA annotation only when species is S. aureus
      {if(species == "sau") {
        geom_shadowtext(aes(x = 0.925*L-2.5,
                            y = H+(15/2),
                            label = ifelse(clusters_summary[[i]]$AMR[1] == "MRSA",
                                           "R",
                                           "S")),
                        color = ifelse(clusters_summary[[i]]$AMR[1] == "MRSA",
                                       "#FF69B4",
                                       "#FFD700"),
                        size = 21.5,
                        fontface = "bold",
                        bg.colour = NA)
      }} +
      theme_void()
    
    
    pdf(file = file.path(output_dir,paste0("Cluster",clusters_summary[[i]]$cluster,".pdf")),
        width = as.numeric(pdf_export_size$width[pdf_export_size$num_circle == plot_members_num]),
        height = as.numeric(pdf_export_size$height[pdf_export_size$num_circle == plot_members_num]))
    print(cluster_plot)
    dev.off()
    
    # Store the plot and cluster ID in the results list
    cluster_plots_rds[[i]] <- list(
      cluster_id = clusters_summary[[i]]$cluster,
      plot = cluster_plot
    )
  }
  return(cluster_plots_rds)
}

# Execute the function ----
if(length(clusters_summary) >= 1){
  plot_rds <- clusterplot(output_dir,
                          clusters_summary)
}else{
  plot_rds <- list()
}


saveRDS(plot_rds,
        snakemake@output[["plot_rds"]])
