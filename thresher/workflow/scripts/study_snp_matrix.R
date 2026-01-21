#use dnadiff function in mummer 
#the snp distance between 2 genomes is the average of 2 snps
#for example, when A is reference and B is query, snp = x. When B is reference and A is query, snp = y. 
#snp = (a+b)/2
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
# Input ----
output_list <- list.files(path = snakemake@params[["report_dir"]],
                          pattern = "*_concatenated.report",
                          all.files = TRUE,
                          full.names = TRUE,
                          recursive = TRUE)

metadata_path <- snakemake@params[["metadata"]]

snp_coverage_threshold <- as.numeric(snakemake@params[["snp_coverage_threshold"]])

metadata <- read.table(metadata_path,
                       sep = "\t",
                       header = FALSE)

metadata$V1 <- sapply(metadata$V1, parse_genome_name)

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
  
  ori_output_df <- read.delim(output_list[i],
                              header = FALSE)
  
  reference_genome <- gsub("_concatenated.report",
                           "",
                           basename(output_list[i]))
  
  sorted_output_df <- as.data.frame(matrix(nrow = nrow(ori_output_df)/4,
                                           ncol = 6))
  
  colnames(sorted_output_df) <- c("reference",
                                  "query",
                                  "snp",
                                  "gsnp",
                                  "AlignedBases_reference",
                                  "AlignedBases_query")
  
  sorted_output_df$reference <- reference_genome
  
  for(n in 1:(nrow(ori_output_df)/4)){
    #query
    query_genome_path <- strsplit(ori_output_df$V1[4*n-3],
                                   split = "\\ ")[[1]][2]
    
    # Use the query genome to match in the metadata
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
    # For the query, use the file path to match the real name in the metadata 
    
    sorted_output_df$AlignedBases_query[n] <- strsplit(ori_output_df$V1[4*n-2],
                                                       split = " ")[[1]][strsplit(ori_output_df$V1[4*n-2],
                                                                                  split = " ")[[1]] != ""][3]
  }
  #return the result
  return(sorted_output_df)
}

sum_snp_df <- do.call(rbind,
                      parLapplyLB(cl,
                                  seq_len(length(output_list)),
                                  concatenate_output))
stopCluster(cl)
rm(cl)
#create a data frame with no redundant comparisons
unique_comparisons <- sum_snp_df %>%
  filter(reference != query) %>%
  mutate(subject = pmin(reference,query),
         query = pmax(reference,query)) %>%
  distinct(subject,
           query)

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
  AlignedBases_reference_entry <- sum_snp_df$AlignedBases_reference[sum_snp_df$reference == unique_comparisons$subject[i] &
                                      sum_snp_df$query == unique_comparisons$query[i]]
  
  AlignedBases_reference_pct <- as.numeric(gsub(".*\\((.*)%\\).*", "\\1", AlignedBases_reference_entry))
  
  
  AlignedBases_query_entry <- sum_snp_df$AlignedBases_reference[sum_snp_df$reference == unique_comparisons$query[i] &
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
                    snp = mean(c(sum_snp_df$snp[sum_snp_df$reference==unique_comparisons$subject[i] & sum_snp_df$query == unique_comparisons$query[i]],
                                 sum_snp_df$snp[sum_snp_df$reference==unique_comparisons$query[i] & sum_snp_df$query == unique_comparisons$subject[i]])),
                    gsnp = mean(c(sum_snp_df$gsnp[sum_snp_df$reference==unique_comparisons$subject[i] & sum_snp_df$query == unique_comparisons$query[i]],
                                  sum_snp_df$gsnp[sum_snp_df$reference==unique_comparisons$query[i] & sum_snp_df$query == unique_comparisons$subject[i]])),
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
gc()

saveRDS(sum_snp_df_unique,
        snakemake@output[[1]])
