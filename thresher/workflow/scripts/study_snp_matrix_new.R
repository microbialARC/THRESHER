#use dnadiff function in mummer 
#the snp distance between 2 genomes is the average of 2 snps
#for example, when A is reference and B is query, snp = x. When B is reference and A is query, snp = y. 
#snp = (a+b)/2
library(dplyr)
library(data.table)
library(parallel)
# Helper function ----
# Replace any non-alphanumeric characters (except dash/underscore) with underscore in the first column
parse_genome_name <- function(original_name) {
  parsed_name <- gsub("[^a-zA-Z0-9._-]", "_", original_name)
  parsed_name <- gsub("_+", "_", parsed_name)
  if (parsed_name == "Reference") {
    parsed_name <- "Reference_Genome"
  }
  trimws(parsed_name, whitespace = "_")
}

# Function to summarize concatenated report and update the study snps matrix

update_study_snp_matrix <- function(output_list,
                                    original_snp_matrix_path,
                                    snp_coverage_threshold,
                                    ncores){
  
  ## Summarize the concatenated output ----
  
  sum_snp_df <- do.call(rbind,
                        mclapply(output_list,
                                 function(output_entry){
                                   
                                   ori_output_df <- read.delim(output_entry,
                                                               header = FALSE)
                                   
                                   sorted_output_df <- as.data.frame(matrix(nrow = nrow(ori_output_df)/4,
                                                                            ncol = 6))
                                   
                                   colnames(sorted_output_df) <- c("reference",
                                                                   "query",
                                                                   "snp",
                                                                   "gsnp",
                                                                   "AlignedBases_reference",
                                                                   "AlignedBases_query")
                                   
                                   
                                   for(n in 1:(nrow(ori_output_df)/4)){
                                     #reference
                                     reference_genome_path <- strsplit(ori_output_df$V1[4*n-3],
                                                                       split = " ")[[1]][1]
                                     reference_genome_name <- metadata$V1[metadata$V3 == reference_genome_path]
                                       
                                     
                                     sorted_output_df$reference[n] <- reference_genome_name
                                     #query
                                     query_genome_path <- strsplit(ori_output_df$V1[4*n-3],
                                                                   split = " ")[[1]][2]
                                     
                                     query_genome_name <- metadata$V1[metadata$V3 == query_genome_path]
                                       
                                     sorted_output_df$query[n] <- query_genome_name
                                     
                                     #snp
                                     sorted_output_df$snp[n] <- as.numeric(unique(unlist(strsplit(ori_output_df$V1[4*n-1],
                                                                                                  split = " ")[[1]][strsplit(ori_output_df$V1[4*n-1],
                                                                                                                             split = " ")[[1]] != ""])[unlist(strsplit(ori_output_df$V1[4*n-1],
                                                                                                                                                                       split = " ")[[1]][strsplit(ori_output_df$V1[4*n-1],
                                                                                                                                                                                                  split = " ")[[1]] != ""]) != "TotalSNPs"]))
                                     #gsnp
                                     sorted_output_df$gsnp[n] <- as.numeric(unique(unlist(strsplit(ori_output_df$V1[4*n],
                                                                                                   split = " ")[[1]][strsplit(ori_output_df$V1[4*n],
                                                                                                                              split = " ")[[1]] != ""])[unlist(strsplit(ori_output_df$V1[4*n],
                                                                                                                                                                        split = " ")[[1]][strsplit(ori_output_df$V1[4*n],
                                                                                                                                                                                                   split = " ")[[1]] != ""]) != "TotalGSNPs"]))
                                     
                                     #AlignedBases_reference
                                     sorted_output_df$AlignedBases_reference[n] <- strsplit(ori_output_df$V1[4*n-2],
                                                                                            split = " ")[[1]][strsplit(ori_output_df$V1[4*n-2],
                                                                                                                       split = " ")[[1]] != ""][2]
                                     #AlignedBases_query
                                     sorted_output_df$AlignedBases_query[n] <- strsplit(ori_output_df$V1[4*n-2],
                                                                                        split = " ")[[1]][strsplit(ori_output_df$V1[4*n-2],
                                                                                                                   split = " ")[[1]] != ""][3]
                                   }
                                   #return the result
                                   return(sorted_output_df)
                                   
                                 },
                                 mc.cores = ncores))
  
  sum_snp_df <- unique(sum_snp_df)
  
  ## Remove redundancy and calculate the mean ----
  
  # The new SNP matrix will include SNP comparisons of:
  # 1. New genomes vs Original genomes
  # 2. New genomes vs New genomes
  # For each comparison, the two genomes will be compared in both directions (A vs B and B vs A)
  # Meaning that each genome will be used as reference once and query once 
  total_genome <- unique(c(sum_snp_df$reference,
                           sum_snp_df$query))
  new_genome <- gsub("_concatenated.report","",basename(output_list))
  original_genome <- setdiff(total_genome,
                             new_genome)
  
  new_new_comparisons <- data.frame(
    subject = t(combn(new_genome, 2))[,1],
    query = t(combn(new_genome, 2))[,2]
  )
  
  # new genomes vs original genomes 
  new_original_comparisons <- expand.grid(
    subject = new_genome,
    query = original_genome
  )
  
  unique_comparisons <- rbind(new_new_comparisons,
                              new_original_comparisons)
  
  #now calculate the snp distance mean(SNPa,SNPb)
  
  sum_snp_df_unique <- do.call(rbind,
                               mclapply(seq_len(nrow(unique_comparisons)),
                                        function(row_entry){
                                          # If the both percentages are above the snp_coverage_threshold, the snp_quality is "good", else "poor"
                                          AlignedBases_reference_entry <- sum_snp_df$AlignedBases_reference[sum_snp_df$reference == unique_comparisons$subject[row_entry] &
                                                                              sum_snp_df$query == unique_comparisons$query[row_entry]]
                                          
                                          AlignedBases_reference_pct <- as.numeric(gsub(".*\\((.*)%\\).*", "\\1", AlignedBases_reference_entry))
                                          
                                          AlignedBases_query_entry <- sum_snp_df$AlignedBases_reference[sum_snp_df$reference == unique_comparisons$query[row_entry] &
                                                                                                         sum_snp_df$query == unique_comparisons$subject[row_entry]]
                                          
                                          AlignedBases_query_pct <- as.numeric(gsub(".*\\((.*)%\\).*", "\\1", AlignedBases_query_entry))
                                          
                                          snp_quality_entry <- if(AlignedBases_reference_pct >= snp_coverage_threshold &
                                                                  AlignedBases_query_pct >= snp_coverage_threshold){
                                            "good"
                                          } else {
                                            "poor"
                                          }
                                                                                                         
                                          return(data.table(subject = unique_comparisons$subject[row_entry],
                                                            query = unique_comparisons$query[row_entry],
                                                            snp = mean(c(sum_snp_df$snp[sum_snp_df$reference==unique_comparisons$subject[row_entry] & sum_snp_df$query == unique_comparisons$query[row_entry]],
                                                                         sum_snp_df$snp[sum_snp_df$reference==unique_comparisons$query[row_entry] & sum_snp_df$query == unique_comparisons$subject[row_entry]])),
                                                            gsnp = mean(c(sum_snp_df$gsnp[sum_snp_df$reference==unique_comparisons$subject[row_entry] & sum_snp_df$query == unique_comparisons$query[row_entry]],
                                                                          sum_snp_df$gsnp[sum_snp_df$reference==unique_comparisons$query[row_entry] & sum_snp_df$query == unique_comparisons$subject[row_entry]])),
                                                            AlignedBases_reference = AlignedBases_reference_entry,
                                                            AlignedBases_query = AlignedBases_query_entry,
                                                            snp_quality = snp_quality_entry))
                                          
                                          
                                        },
                                        mc.cores = ncores))
  
  ## Update the original SNP matrix to get the updated study SNP matrix using new_snps and new_full mode  ----
  
  original_snp_matrix <- readRDS(original_snp_matrix_path)
  
  study_snp_martix_new <- rbind(sum_snp_df_unique,
                               original_snp_matrix)
  
  # Return the matrix 
  
  return(study_snp_martix_new)
  
}



# Input from Snakemake
output_list <- list.files(path = snakemake@params[["report_dir"]],
                          pattern = "*_concatenated.report",
                          all.files = TRUE,
                          full.names = TRUE,
                          recursive = TRUE)

original_snp_matrix_path <- snakemake@params[["original_snp_matrix"]]

new_metadata_path <- snakemake@params[["new_metadata"]]
original_metadata_path <- snakemake@params[["original_metadata"]]

metadata <- rbind(read.table(new_metadata_path,
                             sep = "\t",
                             header = FALSE),
                  read.table(original_metadata_path,
                             sep = "\t",
                             header = FALSE))
metadata$V1 <- sapply(metadata$V1, parse_genome_name)
ncores <- snakemake@threads
snp_coverage_threshold <- as.numeric(snakemake@params[["snp_coverage_threshold"]])

# Execute the function to get updated study_snp_matrix
study_snp_matrix_new <- update_study_snp_matrix(output_list = output_list,
                                                original_snp_matrix_path = original_snp_matrix_path,
                                                snp_coverage_threshold = snp_coverage_threshold,
                                                ncores = ncores)

# Save the matrix in RDS format

saveRDS(study_snp_matrix_new,
        snakemake@output[["study_snp_matrix_new"]])

