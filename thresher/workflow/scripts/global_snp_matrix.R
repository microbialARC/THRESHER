library(dplyr)
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
# Import from snakemake
output_list <- list.files(path = snakemake@params[["report_dir"]],
                          pattern = "*_concatenated.report",
                          all.files = TRUE,
                          full.names = TRUE,
                          recursive = TRUE)

whatsgnu_list <- list.files(path = snakemake@params[["whatsgnu_dir"]],
                            pattern = "*_WhatsGNU_topgenomes.txt",
                            all.files = TRUE,
                            full.names = TRUE,
                            recursive = TRUE)

metadata_path <- snakemake@params[["metadata"]]

actual_download_topgenomes_path <- snakemake@input[["actual_download_topgenomes"]]

snp_coverage_threshold <- as.numeric(snakemake@params[["snp_coverage_threshold"]])

metadata <- read.table(metadata_path,
                       sep = "\t",
                       header = FALSE)
metadata$V1 <- sapply(metadata$V1, parse_genome_name)

actual_download_topgenomes <- unique(readLines(actual_download_topgenomes_path))

cl <- makeCluster(detectCores())
#export the dplyr library to the cluster object
clusterEvalQ(cl, {
  library(dplyr)
})
#export the data needed for the parallel process
clusterExport(cl,
              c("output_list",
                "metadata"),
              envir = environment())
#export the data needed for the parallel process
concatenate_output <- function(i){
  
  # Check if the file is empty
  # TRUE if there is at least one line with content
  # FALSE if the file is empty
  content_check <- any(nzchar(trimws(
    readLines(output_list[i])
  )))
  
  if (content_check){
    ori_output_df <- read.delim(output_list[i],
                                header = FALSE)
    
    sorted_output_df <- as.data.frame(matrix(nrow = nrow(ori_output_df)/4,
                                             ncol = 6))
    
    colnames(sorted_output_df) <- c("subject",
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
  }else{
    # return NULL for empty files
    return(NULL)
  }
  
}

sum_snp_df <- do.call(rbind,
                      parLapplyLB(cl,
                                  seq_len(length(output_list)),
                                  concatenate_output))
stopCluster(cl)
rm(cl)

#create a data frame with no redundant comparisons
unique_comparisons <- do.call(rbind,
                              lapply(whatsgnu_list,
                                     function(list_entry){
                                       
                                       study_genome <- gsub("_WhatsGNU_topgenomes.txt",
                                                            "",
                                                            basename(list_entry))
                                       # suppressing incomplete line warning
                                       file_lines <- suppressWarnings(readLines(list_entry))
                                       if(length(file_lines) <= 1){
                                         return(NULL)
                                       }else{
                                         # Parse from the already-read lines (skip header)
                                         list_df <- read.table(text = file_lines[-1])
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
                                       }
                                     }))

# Only keep those global genomes that are actually downloaded
unique_comparisons <- unique_comparisons[unique_comparisons$query %in% actual_download_topgenomes,]

#now calculate the snp distance mean(SNPa,SNPb)
cl <- makeCluster(detectCores())
#export the data needed for the parallel process
clusterExport(cl,
              c("sum_snp_df",
                "unique_comparisons",
                "snp_coverage_threshold"),
              envir = environment())

remove_redundancy <- function(i){
  # If the both percentages are above the snp_coverage_threshold, the snp_quality is "good", else "poor"
  AlignedBases_reference_entry <- sum_snp_df$AlignedBases_reference[sum_snp_df$subject == unique_comparisons$subject[i] &
                                                                      sum_snp_df$query == unique_comparisons$query[i]]
  
  AlignedBases_reference_pct <- as.numeric(gsub(".*\\((.*)%\\).*", "\\1", AlignedBases_reference_entry))
  
  
  AlignedBases_query_entry <- sum_snp_df$AlignedBases_reference[sum_snp_df$subject == unique_comparisons$query[i] &
                                                                  sum_snp_df$query == unique_comparisons$subject[i]]
  
  AlignedBases_query_pct <- as.numeric(gsub(".*\\((.*)%\\).*", "\\1", AlignedBases_query_entry))
  
  snp_quality_entry <- if(AlignedBases_reference_pct >= snp_coverage_threshold &
                          AlignedBases_query_pct >= snp_coverage_threshold){
    "good"
  } else {
    "poor"
  }
  
  return(data.frame(subject = unique_comparisons$subject[i],
                    query = unique_comparisons$query[i],
                    snp = mean(c(sum_snp_df$snp[sum_snp_df$subject==unique_comparisons$subject[i] & sum_snp_df$query == unique_comparisons$query[i]],
                                 sum_snp_df$snp[sum_snp_df$subject==unique_comparisons$query[i] & sum_snp_df$query == unique_comparisons$subject[i]])),
                    gsnp = mean(c(sum_snp_df$gsnp[sum_snp_df$subject==unique_comparisons$subject[i] & sum_snp_df$query == unique_comparisons$query[i]],
                                  sum_snp_df$gsnp[sum_snp_df$subject==unique_comparisons$query[i] & sum_snp_df$query == unique_comparisons$subject[i]])),
                    AlignedBases_reference = AlignedBases_reference_entry,
                    AlignedBases_query = AlignedBases_query_entry,
                    snp_quality = snp_quality_entry
                    ))
}

sum_snp_df_unique <- do.call(rbind,
                             parLapplyLB(cl,
                                         seq_len(nrow(unique_comparisons)),
                                         remove_redundancy))
stopCluster(cl)
rm(cl)

saveRDS(sum_snp_df_unique,
        snakemake@output[["global_snp_matrix"]])