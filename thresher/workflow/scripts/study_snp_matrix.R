#use dnadiff function in mummer 
#the snp distance between 2 genomes is the average of 2 snps
#for example, when A is reference and B is query, snp = x. When B is reference and A is query, snp = y. 
#snp = (a+b)/2
library(dplyr)
library(parallel)
output_list <- list.files(path = snakemake@params[["report_dir"]],
                          pattern = "*_concatenated.report",
                          all.files = TRUE,
                          full.names = TRUE,
                          recursive = TRUE)

cl <- makeCluster(detectCores())
#export the dplyr library to the cluster object
clusterEvalQ(cl, {
  library(dplyr)
})
#export the data needed for the parallel process
clusterExport(cl,
              c("output_list"),
              envir = environment())
#export the data needed for the parallel process
concatenate_output <- function(i){
  
  ori_output_df <- read.delim(output_list[i],
                              header = FALSE)
  
  ori_output_df$V1 <- gsub(".fasta|.fna|.fa",
                           "",
                           basename(ori_output_df$V1))
  
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
    sorted_output_df$query[n] <- ori_output_df$V1[4*n-3]
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
                "unique_comparisons"),
              envir = environment())

remove_redundancy <- function(i){
  return(data.frame(subject = unique_comparisons$subject[i],
                    query = unique_comparisons$query[i],
                    snp = mean(c(sum_snp_df$snp[sum_snp_df$reference==unique_comparisons$subject[i] & sum_snp_df$query == unique_comparisons$query[i]],
                                 sum_snp_df$snp[sum_snp_df$reference==unique_comparisons$query[i] & sum_snp_df$query == unique_comparisons$subject[i]])),
                    gsnp = mean(c(sum_snp_df$gsnp[sum_snp_df$reference==unique_comparisons$subject[i] & sum_snp_df$query == unique_comparisons$query[i]],
                                  sum_snp_df$gsnp[sum_snp_df$reference==unique_comparisons$query[i] & sum_snp_df$query == unique_comparisons$subject[i]])),
                    AlignedBases_reference = sum_snp_df$AlignedBases_reference[sum_snp_df$reference == unique_comparisons$subject[i] &
                                                                                 sum_snp_df$query == unique_comparisons$query[i]],
                    AlignedBases_query = sum_snp_df$AlignedBases_reference[sum_snp_df$reference == unique_comparisons$query[i] &
                                                                             sum_snp_df$query == unique_comparisons$subject[i]]))
  
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
