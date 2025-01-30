#Visualize the snp distance
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(stringr)
library(ggforce)
# Get the input from Snakemake
output_path <- snakemake@params[["output_dir"]]
# Make the directory if not existed
system(paste0("mkdir -p ",output_path))
study_matrix_path <- snakemake@input[["study_snp_matrix"]]
global_matrix_path <- snakemake@input[["global_snp_matrix"]]

# Function to visualize the gSNP distributions ----

snp_dist_visual <- function(study_matrix_path,
                            global_matrix_path,
                            output_path){
  ## Read and sort the SNP matrices ----
  
  process_matrix <- function(path,
                             category,
                             ref_col = "AlignedBases_reference",
                             query_col = "AlignedBases_query") {
    readRDS(path) %>%
      mutate(
        AlignedBases_reference = as.numeric(str_extract(.data[[ref_col]], "(?<=\\()(.*?)(?=%\\))")),
        AlignedBases_query = as.numeric(str_extract(.data[[query_col]], "(?<=\\()(.*?)(?=%\\))")),
        AlignedPercentage = (AlignedBases_reference + AlignedBases_query) / 2,
        category = category
      ) %>% 
      select(gsnp, AlignedPercentage, category)
  }
  
  matrix <- rbind(
    process_matrix(study_matrix_path, "Study"),
    process_matrix(global_matrix_path, "Global")
  )
  matrix$category <- factor(matrix$category,
                            levels = c("Study","Global"))
  ## generate the plots ----
  
  plot <- ggplot(matrix[matrix$gsnp <= 500,],
                        aes(x = category,
                            y = gsnp)) +
    geom_violin(fill = "transparent",
                color = "transparent",
                scale = "count",
                trim = TRUE) + 
    geom_sina(
      aes(color = AlignedPercentage),
      position = "dodge",
      size = 7.5
    ) +
    scale_y_continuous(
      name = "SNP Distance",
      expand = expansion(mult = c(0, 0)),
      breaks = pretty_breaks(),
    ) + 
    scale_color_viridis_c(name = "% Aligned Bases") +
    theme(
      legend.key.size = unit(5, 'cm'),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 10),
      axis.ticks = element_blank(),
      axis.title.y.left = element_text(colour = "black",
                                       size = 150,
                                       face = "bold"),
      axis.text.y.left = element_text(colour = "black",
                                      size = 100),
      axis.title.x = element_blank(),
      axis.text.x = element_text(colour = "black",
                                 size = 150,
                                 face = "bold"),
      axis.line =  element_blank(),
      legend.position = "top",
      legend.title = element_text(size=100),
      legend.text = element_text(size=75),
      plot.background = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent")
    )
  
  pdf(file=file.path(output_path,"SNP_Distance.pdf"),
      width=40,
      height=50)
  print(plot)
  dev.off() 
}
# use the function
snp_dist_visual(study_matrix_path,
                global_matrix_path,
                output_path)
