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
thresher_output <- snakemake@params[["thresher_output"]]
output_dir <- snakemake@params[["output_dir"]]
new_metadata_path <- snakemake@params[["new_metadata"]]
original_metadata_path <- snakemake@params[["original_metadata"]]
new_mlst_results_path <- snakemake@input[["new_mlst_results"]]
original_mlst_results_path <- snakemake@params[["original_mlst_results"]]

new_peak_strains_rds_path <- snakemake@input[["new_peak_strains_rds"]]
original_peak_strains_rds_path <- snakemake@params[["original_peak_strains_rds"]]
new_plateau_strains_rds_path <- snakemake@input[["new_plateau_strains_rds"]]
original_plateau_strains_rds_path <- snakemake@params[["original_plateau_strains_rds"]]
new_discrepancy_strains_rds_path <- snakemake@input[["new_discrepancy_strains_rds"]]
original_discrepancy_strains_rds_path <- snakemake@params[["original_discrepancy_strains_rds"]]
new_global_strains_rds_path <- snakemake@input[["new_global_strains_rds"]]
original_global_strains_rds_path <- snakemake@params[["original_global_strains_rds"]]
# Input ----
original_cluster_sum <- readRDS(file.path(thresher_output,"thresher","output","clusters_summary.RDS"))
original_cluster_csv <- read.csv(file.path(thresher_output,"thresher","output","clusters_summary.csv"))

# Get the endpoint method from the previous transmission cluster metadata
endpoint <- original_cluster_sum$method

new_strains_rds <- switch(
  endpoint,
  "plateau" = readRDS(new_plateau_strains_rds_path),
  "peak" = readRDS(new_peak_strains_rds_path),
  "global" = readRDS(new_global_strains_rds_path),
  "discrepancy" = readRDS(new_discrepancy_strains_rds_path)
)

original_strains_rds <- switch(
  endpoint,
  "plateau" = readRDS(original_plateau_strains_rds_path),
  "peak" = readRDS(original_peak_strains_rds_path),
  "global" = readRDS(original_global_strains_rds_path),
  "discrepancy" = readRDS(original_discrepancy_strains_rds_path)
)

# Combine original and new metadata
original_metadata <- read.csv(original_metadata_path,
                              header = FALSE,
                              sep = "\t")

new_metadata <- read.csv(new_metadata_path,
                         header = FALSE,
                         sep = "\t")

metadata <- rbind(original_metadata, new_metadata) %>%
  select(V1,V4,V5)
# Combine original and new MLST results
original_mlst_results <- read.csv(original_mlst_results_path,
                                   header = TRUE,
                                   sep = "\t")

new_mlst_results <- read.csv(new_mlst_results_path,
                              header = TRUE,
                              sep = "\t")

mlst_results <- rbind(original_mlst_results, new_mlst_results)
# If MRSA folder exists
original_mrsa_result_path <- file.path(thresher_output,"blastx","mrsa","output","summary","blastx_MRSA_strains.csv")
if(file.exists(original_mrsa_result_path)){
  # Read the original MRSA results
  original_mrsa_results <- read.csv(original_mrsa_result_path)
  # Read the new MRSA results
  new_mrsa_result_path <- file.path(output_dir,"blastx","mrsa","output","summary",paste0("blastx_MRSA_",endpoint,"_strains.csv"))
  new_mrsa_results <- read.csv(new_mrsa_result_path)
  # Combine both to get the temporary MRSA results
  # To determine whether a strain is MRSA positive or not
  # If a strain has no less than half of its genomes being MRSA positive, then this strain is considered MRSA positive
  mrsa_results_tmp <- rbind(original_mrsa_results, new_mrsa_results)
  
  mrsa_results <- do.call(rbind,
                          lapply(unique(mrsa_results_tmp$strain), function(strain_entry) {
                            
                            strain_mrsa_status <- mrsa_results_tmp$MRSA[mrsa_results_tmp$strain == strain_entry]
                            
                            mrsa_positive_num <- sum(strain_mrsa_status == "MRSA")
                            total_genomes_num <- length(strain_mrsa_status)
                            
                            if(mrsa_positive_num / total_genomes_num >= 0.5){
                              return(data.frame(
                                strain = strain_entry,
                                MRSA = "MRSA"
                              ))
                            }else{
                              return(data.frame(
                                strain = strain_entry,
                                MRSA = "MSSA"
                              ))
                            }
                          }))
  
}
# Update the clusters summary using original cluster summary and original endpoint method ----
new_clusters_summary <- lapply(unique(new_strains_rds[["new_strain_compositions"]]$strain_id), function(strain) {
  
  strain_genomes <- new_strains_rds[["new_strain_compositions"]]$genome[new_strains_rds[["new_strain_compositions"]]$strain_id == strain]
  strain_patients <- unique(metadata$V4[metadata$V1 %in% strain_genomes])
  
  if(length(strain_patients) > 1) {
    return(list(
      cluster = NA,
      strain = strain,
      genomes = strain_genomes,
      MLST = unique(mlst_results[[3]][mlst_results$genome %in% strain_genomes]),
      MRSA = if(file.exists(original_mrsa_result_path)) unique(mrsa_results$MRSA[mrsa_results$strain == strain]) else "N/A",
      patients = strain_patients,
      first_seen = format(min(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])), "%y-%m-%d"),
      last_seen = format(max(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])), "%y-%m-%d"),
      persistence = as.integer(max(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])) - min(as.Date(metadata$V5[metadata$V1 %in% strain_genomes])))
    ))
  }else{
    return(NULL)
  }
}) %>% 
  Filter(Negate(is.null), .)

if(original_cluster_csv[1,1] == "No Clusters Found"){
  # No clusters found in previous analysis
  # The whole analysis is created from scratch
  original_clusters_summary <- list()
  # Set max_existing_cluster_id to 0
  max_existing_cluster_id <- 0
}else{
  original_clusters_summary <- original_cluster_sum$clusters
  # Cluster ID assignment based on original clusters summary
  max_existing_cluster_id <- max(unlist(sapply(original_clusters_summary, function(cluster) { as.integer(cluster$cluster) })))
}

# Iterate through new_clusters_summary and assign cluster IDs based on original_clusters_summary and max_existing_cluster_id
for(cluster_idx in seq_along(new_clusters_summary)){
  
  cluster_entry <- new_clusters_summary[[cluster_idx]]
  # The genomes in this cluster
  cluster_genomes <-cluster_entry$genomes
  # Check if any genome in this cluster exists in the original cluster
  matched_original_cluster_id <- unlist(sapply(original_clusters_summary, function(original_cluster_entry) {
    if(all(original_cluster_entry$genomes %in% cluster_genomes)){
      return(original_cluster_entry$cluster)
    }
  }))
  
  # If there is a matched original cluster, assign its cluster ID
  if(length(matched_original_cluster_id) > 0){
    new_clusters_summary[[cluster_idx]]$cluster <- matched_original_cluster_id[1]
  }else{
    # If no match, assign a new cluster ID
    max_existing_cluster_id <- max_existing_cluster_id + 1
    new_clusters_summary[[cluster_idx]]$cluster <- as.character(max_existing_cluster_id)
  }
}

if(length(new_clusters_summary) > 0){
  # Write the new clusters summary csv
  new_clusters_summary_csv <- do.call(rbind,
                                  lapply(new_clusters_summary,
                                         function(cluster) {
                                           
                                           return(data.frame(
                                             cluster = cluster$cluster,
                                             strain = cluster$strain,
                                             MLST = paste(cluster$MLST,collapse = "|"),
                                             MRSA = cluster$MRSA,
                                             genomes = paste(cluster$genomes, collapse = "|"),
                                             patients = paste(cluster$patients, collapse = "|"),
                                             first_seen = cluster$first_seen,
                                             last_seen = cluster$last_seen,
                                             persistence = cluster$persistence
                                           ))
                                         }))
  
  write.csv(new_clusters_summary_csv,
            row.names = FALSE,
            file = snakemake@output[["clusters_summary_csv"]],
            quote = FALSE)

}else{
  new_clusters_summary_csv <- "No Clusters Found"
  write.csv(new_clusters_summary_csv,
            file = snakemake@output[["clusters_summary_csv"]],
            row.names = FALSE,
            quote = FALSE)
}


# Add the endpoint method used to cluster summary
clusters_summary_rds <- list(
  method = endpoint,
  clusters = new_clusters_summary
)

saveRDS(clusters_summary_rds,
        snakemake@output[["clusters_summary_rds"]])

# Write thie genome summary csv

genomes_summary_csv <- do.call(rbind,
                               lapply(new_strains_rds[["new_strain_compositions"]]$genome,
                                      function(genome_entry){
                                        # Is this genome a new genome
                                        is_new_genome <- genome_entry %in% new_metadata$V1
                                        # New strain ID of this genome
                                        new_strain_id <- new_strains_rds[["new_strain_compositions"]]$strain_id[new_strains_rds[["new_strain_compositions"]]$genome == genome_entry]
                                        # Original strain ID of this genome if this genome is not new
                                        original_strain_id <- if(is_new_genome){
                                          "New-Genome"
                                        }else{
                                          original_strains_rds[["strains"]]$strain_id[original_strains_rds[["strains"]]$genome == genome_entry]
                                        }
                                        # Original strain ID of this genome if this genome is not new and also belongs to a cluster
                                        
                                        
                                        if(is_new_genome){
                                          original_cluster_id <- "New-Genome"
                                        }else{
                                          original_cluster_id <- unlist(sapply(original_clusters_summary,
                                                                               function(cluster_entry){
                                                                                 if(genome_entry %in% cluster_entry$genomes){
                                                                                   return(cluster_entry$cluster)
                                                                                 }
                                                                               }))
                                          if(length(original_cluster_id) == 0){
                                            original_cluster_id <- "Non-Cluster"
                                          }
                                        }
                                        
                                        new_cluster_id <- unlist(sapply(new_clusters_summary,
                                                                        function(cluster_entry){
                                                                          if(genome_entry %in% cluster_entry$genomes){
                                                                            return(cluster_entry$cluster)
                                                                          }
                                                                        }))
                                        if(length(new_cluster_id) == 0){
                                          new_cluster_id <- "Non-Cluster"
                                        }
                                        
                                        
                                        return(
                                          data.frame(
                                            genome_name = genome_entry,
                                            genome_category = if(is_new_genome) "new" else "original",
                                            original_strain = original_strain_id,
                                            new_strain = new_strain_id,
                                            original_cluster = original_cluster_id,
                                            new_cluster = new_cluster_id
                                          )
                                        )
                                      }))

# Export the genomes summary csv
write.csv(genomes_summary_csv,
          file = snakemake@output[["genomes_summary_csv"]],
          row.names = FALSE,
          quote = FALSE)

# Function to draw the cluster plot ----
clusterplot <- function(output_dir,
                        clusters_summary){
  
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
    cluster_amr <- clusters_summary[[i]]$MRSA
    
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
      {if(file.exists(original_mrsa_result_path)) {
        geom_shadowtext(aes(x = 0.925*L-2.5,
                            y = H+(15/2),
                            label = ifelse(clusters_summary[[i]]$MRSA[1] == "MRSA",
                                           "R",
                                           "S")),
                        color = ifelse(clusters_summary[[i]]$MRSA[1] == "MRSA",
                                       "#FF69B4",
                                       "#FFD700"),
                        size = 21.5,
                        fontface = "bold",
                        bg.colour = NA)
      }} +
      theme_void()
    
    
    pdf(file = file.path(output_dir,"plots", paste0("Cluster",clusters_summary[[i]]$cluster,".pdf")),
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
# Go to plots output dir
invisible(dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE))

if(length(new_clusters_summary) >= 1){
  plot_rds <- clusterplot(output_dir = output_dir,
                          clusters_summary = new_clusters_summary)
}else{
  plot_rds <- list()
}


saveRDS(plot_rds,
        snakemake@output[["plot_rds"]])
