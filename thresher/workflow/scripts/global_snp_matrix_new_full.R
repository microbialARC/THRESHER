#use dnadiff function in mummer 
#the snp distance between 2 genomes is the average of 2 snps
#for example, when A is reference and B is query, snp = x. When B is reference and A is query, snp = y. 
#snp = (a+b)/2
library(dplyr)
library(data.table)
library(parallel)

# Function to summarize concatenated report and update the global snps matrix
update_global_snp_matrix <- function(output_list,
                                     metadata,
                                     original_snp_matrix_path,
                                     actual_download_topgenomes_path,
                                     whatsgnu_list,
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
                                     
                                     # If the genome was downloaded to datasets_topgenomes, it is a global genome
                                     reference_genome_name <- if(grepl("GCA_",reference_genome_path) && grepl("datasets_topgenomes",reference_genome_path)){
                                       # If reference is global genome, use the base name without extension
                                       tools::file_path_sans_ext(basename(reference_genome_path))
                                       
                                     }else{
                                       # If reference is study genome, match the metadata to find the real name
                                       metadata$V1[metadata$V3 == reference_genome_path]
                                     }
                                     
                                     sorted_output_df$subject[n] <- reference_genome_name
                                     
                                     #query
                                     query_genome_path <- strsplit(ori_output_df$V1[4*n-3],
                                                                   split = " ")[[1]][2]
                                     # If the genome was downloaded to datasets_topgenomes, it is a global genome

                                     query_genome_name <- if(grepl("GCA_",query_genome_path) && grepl("datasets_topgenomes",query_genome_path)){
                                       # If reference is global genome, use the base name without extension
                                       tools::file_path_sans_ext(basename(query_genome_path))
                                       
                                     }else{
                                       # If reference is study genome, match the metadata to find the real name
                                       metadata$V1[metadata$V3 == query_genome_path]
                                     }
                                     
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
  # Unlike the study_snp_matrix, in the global_snp_matrix, the reference is always the study genome
  # While the query is always the global genome
  # Thus there is a different way to create the dataframe without redundancy
  unique_comparisons <- do.call(rbind,
                                lapply(whatsgnu_list,
                                       function(list_entry){
                                         
                                         study_genome <- gsub("_WhatsGNU_topgenomes.txt",
                                                              "",
                                                              basename(list_entry))
                                         
                                         list_df <- read.table(list_entry, skip=1)
                                         
                                         # exclude the study genomes themselves
                                         list_df <- list_df[!(list_df$V1 %in% metadata$V2),]
                                         if(nrow(list_df)>0){
                                           query_name <- sapply(list_df$V1[1:nrow(list_df)],
                                                                function(name_entry){
                                                                  
                                                                  name_split <- strsplit(name_entry,
                                                                                         split = "\\_")[[1]]
                                                                  
                                                                  gca_idx <- which(grepl("GCA",name_split))
                                                                  
                                                                  accession_entry <- paste0("GCA_",name_split[gca_idx+1])
                                                                  
                                                                  return(accession_entry)
                                                                })
                                           genome_comparison_df <- data.frame(
                                             subject = study_genome,
                                             query = query_name
                                           )
                                         }else{
                                           genome_comparison_df <- NULL
                                         }
                                         
                                         return(genome_comparison_df)
                                         
                                       }))
  
  # Only keep those global genomes that are actually downloaded
  actual_download_topgenomes <- unique(readLines(actual_download_topgenomes_path))
  unique_comparisons <- unique_comparisons[unique_comparisons$query %in% actual_download_topgenomes,]
  
  
  sum_snp_df_unique <- do.call(rbind,
                               mclapply(seq_len(nrow(unique_comparisons)),
                                        function(row_entry){
                                          
                                          return(data.table(subject = unique_comparisons$subject[row_entry],
                                                            query = unique_comparisons$query[row_entry],
                                                            snp = mean(c(sum_snp_df$snp[sum_snp_df$reference==unique_comparisons$subject[row_entry] & sum_snp_df$query == unique_comparisons$query[row_entry]],
                                                                         sum_snp_df$snp[sum_snp_df$reference==unique_comparisons$query[row_entry] & sum_snp_df$query == unique_comparisons$subject[row_entry]])),
                                                            gsnp = mean(c(sum_snp_df$gsnp[sum_snp_df$reference==unique_comparisons$subject[row_entry] & sum_snp_df$query == unique_comparisons$query[row_entry]],
                                                                          sum_snp_df$gsnp[sum_snp_df$reference==unique_comparisons$query[row_entry] & sum_snp_df$query == unique_comparisons$subject[row_entry]])),
                                                            AlignedBases_reference = sum_snp_df$AlignedBases_reference[sum_snp_df$reference == unique_comparisons$subject[row_entry] &
                                                                                                                         sum_snp_df$query == unique_comparisons$query[row_entry]],
                                                            AlignedBases_query = sum_snp_df$AlignedBases_reference[sum_snp_df$reference == unique_comparisons$query[row_entry] &
                                                                                                                     sum_snp_df$query == unique_comparisons$subject[row_entry]]))
                                          
                                          
                                        },
                                        mc.cores = ncores))
  
  ## Update the original SNP matrix to get the updated global SNP matrix using new_full mode  ----
  original_snp_matrix <- readRDS(original_snp_matrix_path)
  
  global_snp_martix_new <- rbind(sum_snp_df_unique,
                                 original_snp_matrix)
  
  # Return the matrix 
  return(global_snp_martix_new)
  
}

# Input from Snakemake
output_list <- list.files(path = snakemake@params[["report_dir"]],
                          pattern = "*_concatenated.report",
                          all.files = TRUE,
                          full.names = TRUE,
                          recursive = TRUE)
original_snp_matrix_path <- snakemake@params[["original_snp_matrix"]]
whatsgnu_list <- list.files(path = snakemake@params[["whatsgnu_dir"]],
                            pattern = "*_WhatsGNU_topgenomes.txt",
                            all.files = TRUE,
                            full.names = TRUE,
                            recursive = TRUE)

new_metadata_path <- snakemake@params[["new_metadata"]]
original_metadata_path <- snakemake@params[["original_metadata"]]

metadata <- rbind(read.table(new_metadata_path,
                             sep = "\t",
                             header = FALSE),
                  read.table(original_metadata_path,
                             sep = "\t",
                             header = FALSE))

actual_download_topgenomes_path <- snakemake@input[["actual_download_topgenomes"]]

ncores <- snakemake@threads

# Execute the function to get updated global_snp_matrix
global_snp_martix_new <- update_global_snp_matrix(output_list = output_list,
                                                  metadata = metadata,
                                                  original_snp_matrix_path = original_snp_matrix_path,
                                                  actual_download_topgenomes_path = actual_download_topgenomes_path,
                                                  whatsgnu_list = whatsgnu_list,
                                                  ncores = ncores)

# Save the matrix in RDS format
saveRDS(global_snp_martix_new,
        snakemake@output[["global_snp_matrix_new"]])

