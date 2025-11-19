library(dplyr)
library(parallel)
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

metadata <- read.table(metadata_path,
                       sep = "\t",
                       header = FALSE)

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
    
    
    reference_genome_name <- if(grepl("GCA_",reference_genome_path)){
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
    
    query_genome_name <- if(grepl("GCA_",query_genome_path)){
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
                                       list_df <- read.table(list_entry, skip=1)
                                       return(data.frame(
                                         subject = study_genome,
                                         query = list_df$V1[grepl("GCA_", list_df$V1)][1:10]
                                       ))
                                     }))

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
                    snp = mean(c(sum_snp_df$snp[sum_snp_df$subject==unique_comparisons$subject[i] & sum_snp_df$query == unique_comparisons$query[i]],
                                 sum_snp_df$snp[sum_snp_df$subject==unique_comparisons$query[i] & sum_snp_df$query == unique_comparisons$subject[i]])),
                    gsnp = mean(c(sum_snp_df$gsnp[sum_snp_df$subject==unique_comparisons$subject[i] & sum_snp_df$query == unique_comparisons$query[i]],
                                  sum_snp_df$gsnp[sum_snp_df$subject==unique_comparisons$query[i] & sum_snp_df$query == unique_comparisons$subject[i]])),
                    AlignedBases_reference = sum_snp_df$AlignedBases_reference[sum_snp_df$subject == unique_comparisons$subject[i] &
                                                                                 sum_snp_df$query == unique_comparisons$query[i]],
                    AlignedBases_query = sum_snp_df$AlignedBases_reference[sum_snp_df$subject == unique_comparisons$query[i] &
                                                                             sum_snp_df$query == unique_comparisons$subject[i]]
                    ))
}

sum_snp_df_unique <- do.call(rbind,
                             parLapplyLB(cl,
                                         seq_len(nrow(unique_comparisons)),
                                         remove_redundancy))
stopCluster(cl)
rm(cl)

saveRDS(sum_snp_df_unique,
        snakemake@output[[1]])