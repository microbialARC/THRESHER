# The script to generate the profile of the provided reference genome
# 1. Get the entropy as proxy to imply mutation possibility
# 2. Infer the MGEs(Mobile genetic elements) using coverage
# Eight helper functions to complete the full analysis. The manual and detail instructions are coming soon
# This script is used as a command line tool 
# Libraries ----
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
suppressMessages(library(parallel))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(entropy))
suppressMessages(library(cowplot))
suppressMessages(library(ggnewscale))
suppressMessages(library(ggExtra))

# Input from arguments ----
# Command-line Rscript tool. Accepts these flags when run standalone:
#   --input_genome <path>   : reference assembly (FASTA; .fasta/.fna). Required.
#   --gff <path>            : annotation file (GFF/GFF3). Required.
#   --snp <dir>             : directory with SNP/coord outputs (expects *.snps and *.coords from MUMmer). Required.
#   --cpus <int>            : number of CPU cores to use (optional; defaults to detectCores()).
#   --output <dir>          : directory where all outputs (RDS/CSV/PDF) will be written. Required.
#
# Snakemake integration:
#   When run inside a Snakemake rule, Snakemake should invoke Rscript with the same flags
#   (see the pipeline's genome_profiler.smk rule). The script reads commandArgs(trailingOnly=TRUE).
#
# Outputs (written to --output):
#   *_concat.fasta, *_new_pos.RDS, *_position_coverage.RDS, *_snps_sum.RDS, *_chr_bins.RDS,
#   *_entropy.RDS / .csv, *_mges.RDS / .csv / _mges_seq.fasta, *_bin_summary.csv, *_profiler_plot.pdf

# When run inside the Thresher pipeline, Snakemake calls this script and provides the necessary command-line arguments.
# Therefore the argument-parsing code (args and related blocks) is commented out.
# To run the script standalone, edit the file and uncomment those sections or supply the arguments manually.
#args <- commandArgs(trailingOnly = TRUE)
#for (i in seq_along(args)) {
#  if (args[i] == "--input_genome") {
#    input_genome_path <- args[i + 1]
#  } else if (args[i] == "--gff") {
#    gff_path <- args[i + 1]
#  }else if (args[i] == "--snp") {
#    snp_dir <- args[i + 1]
#  }else if (args[i] == "--cpus") {
#    ncores <- as.integer(args[i + 1])
#  }else if (args[i] == "--output") {
#    output_dir <- args[i + 1]
#  }
#}

# Directly get the input from snakemake when run inside the pipeline
# Path to the input genome assembly file
input_genome_path <- snakemake@input[["fna_path"]]
# Path to the gff file of the input genome assembly 
gff_path <- snakemake@input[["gff3_path"]]
# Path to directory where the output files will be saved
output_dir <- snakemake@params[["output_dir"]]
# Make directory if not exists
if(!dir.exists(output_dir)){
  dir.create(output_dir,recursive = TRUE)
}
# Path to the directory of corresponding snp results 
snp_dir <- snakemake@params[["snp_dir"]]
# how many cores will be used
ncores <- snakemake@threads

# Check if required arguments were provided
if (is.null(input_genome_path)) {
  stop("Missing required argument: --input_genome")
}

if (is.null(gff_path)) {
  stop("Missing required argument: --gff")
}

if (is.null(snp_dir)) {
  stop("Missing required argument: --snp")
}

if (is.null(output_dir)) {
  stop("Missing required argument: --output")
}

cat("\n")
cat("Input genome:", input_genome_path, "\n")
cat("GFF file:", gff_path, "\n")
cat("SNP directory:", snp_dir, "\n")
cat("Output Directory:", output_dir, "\n")
cat("Threads: ", ncores, "\n")
cat("\n")

# Helper Functions ----
## Concatenate the reference genome into one contig if not complete genome ----

concat_genome <- function(input_genome_name,
                          input_genome_path,
                          output_dir){
  
  input_genome_fasta <- readDNAStringSet(input_genome_path)
  
  if(length(input_genome_fasta) > 1){
    
    names(input_genome_fasta) <- sapply(seq_along(names(input_genome_fasta)),
                                        function(ctg_idx){
                                          strsplit(names(input_genome_fasta)[ctg_idx],split = " ")[[1]][1]
                                        })
    # get data frame recording the new positions matching to the original contig positions
    
    ctg_pos_count <- 0
    ctg_pos_df <- data.frame()
    
    for(ctg_idx in seq_along(input_genome_fasta)){
      
      ctg_name <- names(input_genome_fasta[ctg_idx])
      ctg_length <- width(input_genome_fasta[ctg_idx])
      
      
      ctg_pos_df <- rbind(ctg_pos_df,
                          data.frame(
                            ctg = ctg_name,
                            length = ctg_length,
                            new_start = ctg_pos_count + 1,
                            new_end = ctg_pos_count + ctg_length 
                          ))
      
      ctg_pos_count <- ctg_pos_count + ctg_length
    }
    # export the reference fasta
    input_genome_fasta_concat <- DNAStringSet(unlist(input_genome_fasta))
    
    names(input_genome_fasta_concat) <- input_genome_name
    
  }else{
    # get the position matching data frame
    # no change to position
    ctg_pos_df <- data.frame(
      ctg = names(input_genome_fasta),
      length = width(input_genome_fasta),
      new_start = 1,
      new_end = width(input_genome_fasta)
    )
    
    input_genome_fasta_concat <- input_genome_fasta
  }
  
  
  writeXStringSet(input_genome_fasta_concat,
                  file.path(output_dir,paste0(input_genome_name,"_concat.fasta")))
  
  
  new_pos <- list(
    concat_size = width(input_genome_fasta_concat),
    ctg_pos_df = ctg_pos_df
  )
  
  saveRDS(new_pos,
          file.path(output_dir,paste0(input_genome_name,"_new_pos.RDS")))
  
  return(new_pos)
  
}

## Calculate the coverage for every site in the concatenated contig ----
get_position_coverage <- function(input_genome_name,
                                  new_pos,
                                  snp_dir,
                                  output_dir,
                                  ncores){
  
  
  ctg_pos_df <- new_pos$ctg_pos_df
  
  chr_size <- new_pos$concat_size
  
  # get the alignment coverage count for every site
  # first concatenate all aln_coord to one data frame
  # and then read all files at once
  aln_coord_list <- list.files(path = snp_dir,
                               pattern = "\\.coords$",
                               all.files = TRUE,
                               full.names = TRUE,
                               recursive = FALSE)
  
  aln_coord <- rbindlist(
    mclapply(seq_along(aln_coord_list), function(list_idx) {
      dt <- fread(aln_coord_list[list_idx],
                  header = FALSE)
      dt[, query_id := list_idx]
      return(dt)
    }, mc.cores = ncores)
  )
  
  calc_cover <- function(ctg_alns, ctg_pos) {
    vapply(ctg_pos, function(pos_entry) {
      ctg_alns[V1 <= pos_entry & V2 >= pos_entry, uniqueN(query_id)]
    }, integer(1))
  }
  
  position_coverage <- rbindlist(
    mclapply(ctg_pos_df$ctg, function(ctg_idx) {
      ctg_info <- ctg_pos_df[ctg_pos_df$ctg == ctg_idx, ]
      ctg_alns <- aln_coord[aln_coord$V3 == ctg_idx, ]
      ctg_pos <- seq_len(ctg_info$length)
      
      if (nrow(ctg_alns) > 0) {
        coverage <- calc_cover(
          ctg_alns,
          ctg_pos
        )
      } else {
        coverage <- rep(0L, length(ctg_pos))
      }
      
      data.table(
        position = (ctg_info$new_start -1)+ ctg_pos ,
        coverage = coverage
      )
    },
    mc.cores = ncores))
  
  saveRDS(position_coverage,
          file.path(output_dir,paste0(input_genome_name,"_position_coverage.RDS")))
  
  return(position_coverage)
}

## Summarize the SNPs into the data frame ----
get_snps_sum <- function(input_genome_name,
                         new_pos,
                         snp_dir,
                         output_dir,
                         ncores){
  
  ctg_pos_df <- new_pos$ctg_pos_df
  # get the list of files containing snps called by mummer
  snps_list <- list.files(path = snp_dir,
                          pattern = "\\.snps$",
                          all.files = TRUE,
                          full.names = TRUE,
                          recursive = FALSE)
  
  snps_sum <- rbindlist(mclapply(snps_list,
                                 function(list_idx){
                                   if(file.exists(list_idx) && file.size(list_idx) > 0){
                                     
                                     snps_entry <- fread(list_idx, header = FALSE)
                                     
                                     if(nrow(snps_entry) == 0){
                                       return(NULL)
                                     }
                                     
                                     genome_snps_sum <- rbindlist(lapply(seq_len(nrow(snps_entry)),
                                                                         function(row_idx){
                                                                           
                                                                           row_ctg <- snps_entry$V5[row_idx]
                                                                           
                                                                           # V1: SNP position in the reference.
                                                                           # V2: Character in the reference.
                                                                           # V3: Character in the query.
                                                                           # V4: SNP Buffer
                                                                           # V5: query FastA ID
                                                                           
                                                                           old_snp_pos <- snps_entry$V1[row_idx]
                                                                           new_snp_pos <- (ctg_pos_df$new_start[ctg_pos_df$ctg == row_ctg] - 1) + old_snp_pos
                                                                           ref_site <- snps_entry$V2[row_idx]
                                                                           new_site <- snps_entry$V3[row_idx]
                                                                           
                                                                           data.frame(
                                                                             snp_pos = new_snp_pos,
                                                                             ref_site = ref_site,
                                                                             new_site = new_site
                                                                           )
                                                                         }))
                                     
                                     return(genome_snps_sum)
                                     
                                   }else{
                                     return(NULL)
                                   }
                                 },
                                 mc.cores = ncores))
  # In this analysis
  # Insertion and deletion will be ignored
  snps_sum <- dplyr::filter(snps_sum,
                            ref_site != "." & 
                              new_site != ".")
  
  saveRDS(snps_sum,
          file.path(output_dir,paste0(input_genome_name,"_snps_sum.RDS")))
  
  return(snps_sum)
}

## Binning the concatenated contig ----
get_chr_bins <- function(input_genome_name,
                         new_pos,
                         gff_path,
                         output_dir){
  
  ctg_pos_df <- new_pos$ctg_pos_df
  
  chr_size <- new_pos$concat_size
  
  # Create use the coding sequences from annotation as bins
  # Annotate the bins using gff files
  # Load and filter gene annotations
  ref_gff <- read.delim(gff_path,
                        comment.char = "#",
                        header = FALSE) %>%
    filter(V3 == "CDS")
  
  # Add new contig positions
  ref_gff <- ref_gff %>%
    mutate(gene_id = paste0("gene", 1:n())) %>%
    rowwise() %>%
    mutate(
      ctg_offset = ctg_pos_df$new_start[ctg_pos_df$ctg == V1],
      new_start = ctg_offset + V4 - 1,
      new_end = ctg_offset + V5 - 1
    ) %>%
    select(-ctg_offset)
  
  # The overlapping regions generate the new bins that don't belong any of the genes
  
  chr_bins <- data.frame()
  bin_idx <- 1
  
  cds_count <- c()
  
  for(ctg_row in seq_len(nrow(ctg_pos_df))){
    
    
    ctg_start <- ctg_pos_df$new_start[ctg_row]
    ctg_end <- ctg_pos_df$new_end[ctg_row]
    ctg_pos_range <- ctg_start : ctg_end
    
    ref_gff_ctg <- ref_gff[ref_gff$V1 == ctg_pos_df$ctg[ctg_row],]
    
    for(nt_idx in ctg_pos_range){
      
      # check if nt is within cds
      cds_idx <- which(nt_idx >= ref_gff_ctg$new_start & nt_idx < ref_gff_ctg$new_end)
      cds_id <- ref_gff_ctg$gene_id[cds_idx]
      
      if(length(cds_id) == 0){
        # check if this nt is already in any bin
        if(any(nt_idx >= chr_bins$start & nt_idx <= chr_bins$end)){
          next
        }else{
          # check if this chunk is already a bin
          # non-cds by annotation (Scenario1)
          if(nrow(ref_gff_ctg) == 0){
            # if there is no gene annotated in this contig
            # the whole contig is a bin
            
            bin_gene <- "non-cds"
            bin_locus_tag <- "non-cds"
            bin_start <- ctg_start
            bin_end <- ctg_end
            bin_gene <- "non-cds"
            bin_locus_tag <- "non-cds"
            
          }else{
            if(nt_idx < min(ref_gff_ctg$new_start)){
              
              # before all cds in this contig
              next_cds_idx <- min(which(ref_gff_ctg$new_start > nt_idx), na.rm = TRUE)
              next_cds_start <- ref_gff_ctg$new_start[next_cds_idx]
              bin_start <- ctg_start
              bin_end <- next_cds_start - 1
              bin_gene <- "non-cds"
              bin_locus_tag <- "non-cds"
              
            }else if(nt_idx > max(ref_gff_ctg$new_start)){
              
              # after all cds in this contig
              
              before_cds_idx <- max(which(ref_gff_ctg$new_start < nt_idx), na.rm = TRUE)
              before_cds_end <- ref_gff_ctg$new_end[before_cds_idx]
              bin_start <- before_cds_end + 1
              bin_end <- ctg_end
              bin_gene <- "non-cds"
              bin_locus_tag <- "non-cds"
              
            }else{
              
              # in gap between cds
              # next cds position
              next_cds_idx <- min(which(ref_gff_ctg$new_start > nt_idx), na.rm = TRUE)
              next_cds_start <- ref_gff_ctg$new_start[next_cds_idx]
              
              # before cds position
              before_cds_idx <- max(which(ref_gff_ctg$new_start < nt_idx), na.rm = TRUE)
              before_cds_end <- ref_gff_ctg$new_end[before_cds_idx]
              
              bin_start <- before_cds_end + 1
              bin_end <- next_cds_start - 1
              bin_gene <- "non-cds"
              bin_locus_tag <- "non-cds"
            }
          }
          
          chr_bins <- rbind(chr_bins,
                            data.frame(
                              bin_index = bin_idx,
                              start = bin_start,
                              end = bin_end,
                              length = bin_end - bin_start + 1,
                              gene = bin_gene,
                              locus_tag = bin_locus_tag,
                              original_contig = ctg_pos_df$ctg[ctg_row]
                            ))
          
          bin_idx <- bin_idx + 1
        }
      }else{
        
        # if this nt is within cds
        # check if this gene is already marked as bin
        # or if there is overlapping
        # check which gene(s) with the overlapping are not mark
        
        no_bin_cds_id <- setdiff(cds_id,cds_count)
        
        if(length(no_bin_cds_id) == 0){
          # length(no_bin_cds_id) == 0 means all genes are binned
          next
        }else{
          for(gene_id_idx in no_bin_cds_id){
            # if no, use this gene as a bin
            # cds by annotation
            # cds info 
            # the row of the gene in ref_gff_ctg
            idx <- which(ref_gff_ctg$gene_id == gene_id_idx)
            gene_start <- ref_gff_ctg$new_start[idx]
            gene_end <- ref_gff_ctg$new_end[idx]
            
            bin_start <- ref_gff_ctg$new_start[idx]
            bin_end <- ref_gff_ctg$new_end[idx]
            
            gene_info <- strsplit(ref_gff_ctg$V9[idx],
                                  split = "\\;")[[1]]
            
            bin_locus_tag <- gsub("locus_tag=","",grep("locus_tag=",gene_info,value = TRUE))
            bin_gene <- gsub("Name=","",grep("Name=",gene_info,value = TRUE))
            
            chr_bins <- rbind(chr_bins,
                              data.frame(
                                bin_index = bin_idx,
                                start = bin_start,
                                end = bin_end,
                                length = bin_end - bin_start + 1,
                                gene = bin_gene,
                                locus_tag = bin_locus_tag,
                                original_contig = ctg_pos_df$ctg[ctg_row]
                              ))
            
            bin_idx <- bin_idx + 1
            cds_count <- unique(c(cds_count,cds_id))
          }
        }
      }
    }
  }
  
  saveRDS(chr_bins,
          file.path(output_dir,paste0(input_genome_name,"_chr_bins.RDS")))
  
  return(chr_bins)
}

## Calculate the entropy for each site in the in the concatenated contig ----
get_entropy <- function(input_genome_name,
                        new_pos,
                        snps_sum,
                        position_coverage,
                        output_dir,
                        ncores){
  
  # Input used in the function
  snps_sum <- as.data.table(snps_sum)
  position_coverage <- as.data.table(position_coverage)
  ctg_pos_df <- new_pos$ctg_pos_df
  chr_size <- new_pos$concat_size
  
  # Calculate the Shannon Entropy using R package entropy with James-Stein-type shrinkage estimator
  # https://www.jmlr.org/papers/volume10/hausser09a/hausser09a.pdf
  
  # Convert the position_coverage to vector for quicker process
  position_coverage_vector <- setNames(as.integer(position_coverage$coverage), as.integer(position_coverage$position))
  
  entropy_df <- rbindlist(lapply(seq_len(chr_size),
                                 function(pos_entry){
                                   tryCatch({
                                   # If no SNPs in position, entropy is 0
                                   if (sum(snps_sum$snp_pos == pos_entry) == 0) {
                                     return(data.frame(
                                       position = pos_entry,
                                       entropy = 0
                                     ))
                                     
                                   }else{
                                     
                                     pos_cov <- position_coverage_vector[pos_entry]
                                     
                                     pos_snp <- snps_sum[snps_sum$snp_pos == pos_entry,]
                                     # Aggregate SNP counts by position and composition
                                     pos_snp_counts <- aggregate(
                                       rep(1, nrow(pos_snp)),
                                       by = list(pos = pos_snp$snp_pos,
                                                 ref_nt = pos_snp$ref_site,
                                                 new_nt = pos_snp$new_site),
                                       FUN = sum
                                     )
                                     names(pos_snp_counts)[4] <- "count"
                                     
                                     pos_ref_nt <- unique(pos_snp_counts$ref_nt)
                                     
                                     pos_all_nt <- c(pos_ref_nt,pos_snp_counts$new_nt)
                                     
                                     pos_observed_nt_count <- sapply(pos_all_nt, function(nt) {
                                       
                                       if(nt == pos_ref_nt) {
                                         # Reference nt count
                                         1 + (pos_cov - sum(pos_snp_counts$count))
                                       } else {
                                         pos_snp_counts$count[pos_snp_counts$new_nt == nt]
                                       }
                                     })
                                     
                                     pos_entropy <- entropy(pos_observed_nt_count, 
                                                            method = "shrink",
                                                            verbose = FALSE)[1]
                                     return(
                                       data.frame(
                                         position = pos_entry,
                                         entropy = pos_entropy
                                       ))
                                   }
                                   }, error = function(e) {
                                     # Debug output for error
                                     cat("ERROR at position:", pos_entry, "\n")
                                     cat("Error message:", e$message, "\n")
                                     cat("pos_cov:", ifelse(exists("pos_cov"), pos_cov, "NOT_DEFINED"), "\n")
                                     if(exists("pos_snp")) {
                                       cat("pos_snp rows:", nrow(pos_snp), "\n")
                                       cat("pos_snp head:\n")
                                       print(head(pos_snp))
                                     }
                                     if(exists("pos_observed_nt_count")) {
                                       cat("pos_observed_nt_count:\n")
                                       print(pos_observed_nt_count)
                                     }
                                     # Return a placeholder result to continue processing
                                     return(data.frame(
                                       position = pos_entry,
                                       entropy = NA
                                     ))
                                   })
                                 }))
  
  # double check the nrow(entropy_df) == chr_size
  while(nrow(entropy_df) != chr_size){
    # check which position is missing
    missing_entropy_position <- setdiff(1:chr_size,entropy_df$position)
    
    complement_entropy_df <- rbindlist(lapply(missing_entropy_position,
                                              function(pos_entry){
                                                tryCatch({
                                                # If no SNPs in position, entropy is 0
                                                if (sum(snps_sum$snp_pos == pos_entry) == 0) {
                                                  return(data.table(
                                                    position = pos_entry,
                                                    entropy = 0
                                                  ))
                                                  
                                                }else{
                                                  
                                                  pos_cov <- position_coverage_vector[pos_entry]
                                                  
                                                  pos_snp <- snps_sum[snps_sum$snp_pos == pos_entry,]
                                                  # Aggregate SNP counts by position and composition
                                                  pos_snp_counts <- aggregate(
                                                    rep(1, nrow(pos_snp)),
                                                    by = list(pos = pos_snp$snp_pos,
                                                              ref_nt = pos_snp$ref_site,
                                                              new_nt = pos_snp$new_site),
                                                    FUN = sum
                                                  )
                                                  names(pos_snp_counts)[4] <- "count"
                                                  
                                                  pos_ref_nt <- unique(pos_snp_counts$ref_nt)
                                                  
                                                  pos_all_nt <- c(pos_ref_nt,pos_snp_counts$new_nt)
                                                  
                                                  pos_observed_nt_count <- sapply(pos_all_nt, function(nt) {
                                                    
                                                    if(nt == pos_ref_nt) {
                                                      # Reference nt count
                                                      1 + (pos_cov - sum(pos_snp_counts$count))
                                                    } else {
                                                      pos_snp_counts$count[pos_snp_counts$new_nt == nt]
                                                    }
                                                  })
                                                  
                                                  pos_entropy <- entropy(pos_observed_nt_count, 
                                                                         method = "shrink",
                                                                         verbose = FALSE)[1]
                                                  return(
                                                    data.frame(
                                                      position = pos_entry,
                                                      entropy = pos_entropy
                                                    ))
                                                }
                                                }, error = function(e) {
                                                  # Debug output for error in complement section
                                                  cat("ERROR in complement section at position:", pos_entry, "\n")
                                                  cat("Error message:", e$message, "\n")
                                                  cat("pos_cov:", ifelse(exists("pos_cov"), pos_cov, "NOT_DEFINED"), "\n")
                                                  if(exists("pos_snp")) {
                                                    cat("pos_snp rows:", nrow(pos_snp), "\n")
                                                    cat("pos_snp head:\n")
                                                    print(head(pos_snp))
                                                  }
                                                  if(exists("pos_observed_nt_count")) {
                                                    cat("pos_observed_nt_count:\n")
                                                    print(pos_observed_nt_count)
                                                  }
                                                  # Return a placeholder result to continue processing
                                                  return(data.frame(
                                                    position = pos_entry,
                                                    entropy = NA
                                                  ))
                                                })
                                              }))
    
    entropy_df <- unique(rbind(entropy_df,
                        complement_entropy_df))
  }
  
  saveRDS(entropy_df,
          file.path(output_dir,paste0(input_genome_name,"_entropy.RDS")))
  
  # export the csv for the input of evolution simulation
  write.csv(entropy_df,
            file.path(output_dir,paste0(input_genome_name,"_entropy.csv")),
            quote = FALSE,
            row.names = FALSE)
  
  return(entropy_df)
}

## Infer and annotate MGEs (Mobile genetic elements) ----

get_mges <- function(input_genome_name,
                     entropy_df,
                     new_pos,
                     position_coverage,
                     chr_bins,
                     min_mge_size = 100,
                     cov_cutoff = 0.7,
                     output_dir,
                     ncores = ncores){
  
  setDTthreads(ncores)
  
  # To avoid artifact
  # DO NOT go across different contigs and only limits MGE identification to the same contig
  # Consecutive min_mge_size * TRUE is considered a mge
  # Get the coverage value for each position
  # < cov_cutoff TRUE
  
  ctg_pos_df <- new_pos$ctg_pos_df
  
  # Create a idx data table for all MGEs 
  # This might seem redundant but it serves debug purpose
  
  mges_idx_df <- rbindlist(lapply(seq_len(nrow(ctg_pos_df)),
                                  function(row_idx) {
                                    ctg_start <- ctg_pos_df$new_start[row_idx]
                                    ctg_end <- ctg_pos_df$new_end[row_idx]
                                    
                                    # Identify low coverage positions and find runs
                                    position_status <- position_coverage$coverage[ctg_start:ctg_end] / max(position_coverage$coverage) < cov_cutoff
                                    status_rle <- rle(position_status)
                                    
                                    # Calculate run positions
                                    cumulative_end <- cumsum(status_rle$lengths)
                                    cumulative_start <- c(1, cumulative_end[-length(cumulative_end)] + 1)
                                    
                                    # Find MGE indices that meet size criteria
                                    mges_idx <- which(status_rle$values & status_rle$lengths >= min_mge_size)
                                    
                                    if (length(mges_idx) == 0) return(NULL)
                                    
                                    # Create MGE data with adjusted coordinates
                                    data.table(
                                      mge_ctgs = ctg_pos_df$ctg[row_idx],
                                      mge_idx = seq_along(mges_idx),
                                      mge_start = cumulative_start[mges_idx] + ctg_start - 1L,
                                      mge_end = cumulative_end[mges_idx] + ctg_start - 1L
                                    )
                                  }))
  
  # Annotate mges
  # And only mges_df is return to the environment
  
  mges_df <- rbindlist(lapply(seq_len(nrow((mges_idx_df))),function(row_idx){
    
    mge_start <- as.integer(mges_idx_df$mge_start[row_idx])
    mge_end <- as.integer(mges_idx_df$mge_end[row_idx])
    
    mge_bin_idx <- which(chr_bins$start <= mge_end & chr_bins$end >= mge_start)
    mge_bin_idx_report <- if(length(mge_bin_idx) == 1) mge_bin_idx else paste0(min(mge_bin_idx),"-",max(mge_bin_idx))
    mge_bin_gene <- chr_bins$gene[mge_bin_idx]
    # ignore "non-cds" in the report
    mge_bin_gene <- setdiff(mge_bin_gene,"non-cds")
    
    return(data.frame(
      mge_index = paste0(input_genome_name,"_MGE_",row_idx),
      start = mge_start,
      end = mge_end,
      length = mge_end - mge_start + 1,
      bin = mge_bin_idx_report,
      bin_gene = paste(mge_bin_gene,collapse = "; ")
    ))
  }))
  
  # Extract the sequences of mges as the pool for simulation
  # And export the entropy for each site for every MGE
  # The entropy values for each MGE are stored in the individual file while the fasta file contains all MGEs found in the genome 
  
  dir.create(file.path(output_dir,"mge_entropy"), recursive = TRUE, showWarnings = FALSE)
  input_genome_fasta_concat <- readDNAStringSet(file.path(output_dir,paste0(input_genome_name,"_concat.fasta")))
  
  mges_seq <- as.character(sapply(seq_len(nrow(mges_df)),function(row_idx){
    
    
    mge_start <- mges_df$start[row_idx]
    mge_end <- mges_df$end[row_idx]
    
    # Extract entropy values for this MGE
    entropy_values <- entropy_df$entropy[mge_start:mge_end]
    mge_positions <- seq_along(entropy_values)
    
    # Export the data 
    write.csv(
      data.frame(mge_position = mge_positions, entropy = entropy_values),
      file.path(output_dir,"mge_entropy",paste0(input_genome_name,"_MGE_",row_idx,"_entropy.csv")),
      quote = FALSE,
      row.names = FALSE
    )
    
    saveRDS(mges_df,
            file.path(output_dir,paste0(input_genome_name,"_mges.RDS")))
    
    write.csv(mges_df,
              file.path(output_dir,paste0(input_genome_name,"_mges.csv")),
              quote = FALSE,
              row.names = FALSE)
    
    
    mge_fasta <- c(paste0(">",mges_df$mge_index[row_idx]),
                   substr(as.character(input_genome_fasta_concat),mge_start,mge_end))
    
    
    return(mge_fasta)
  }))
  
  writeLines(mges_seq,
             file.path(output_dir,paste0(input_genome_name,"_mges_seq.fasta")))
  
  return(mges_df)
}
## Summarize entropy, coverage, and MGEs, by bins ----
get_bin_sum <- function(input_genome_name,
                        chr_bins,
                        position_coverage,
                        snps_sum,
                        entropy_df,
                        mges_df,
                        output_dir){
  
  # Summarize entropy and coverage by bins
  bin_stats <- rbindlist(lapply(seq_len(nrow(chr_bins)),function(bin_row){
                         
                         bin_idx <- chr_bins$bin_index[bin_row]
                         bin_start <- chr_bins$start[bin_row]
                         bin_end <- chr_bins$end[bin_row]
                         # bin includes the start position but excludes the end position
                         
                         # coverage of every site in the bin 
                         bin_coverage <- as.numeric(position_coverage$coverage[position_coverage$position >= bin_start & position_coverage$position <= bin_end])
                         bin_coverage_mean = as.numeric(mean(bin_coverage))
                         bin_coverage_sd = as.numeric(sd(bin_coverage))
                         bin_coverage_median = as.numeric(median(bin_coverage))
                         bin_coverage_q1 <- as.numeric(quantile(bin_coverage,0.25))
                         bin_coverage_q3 <- as.numeric(quantile(bin_coverage,0.75))
                         
                         # entropy of every site in the bin
                         bin_entropy <- as.numeric(entropy_df$entropy[entropy_df$position >= bin_start &  entropy_df$position <= bin_end])
                         bin_entropy_mean = as.numeric(mean(bin_entropy))
                         bin_entropy_sd = as.numeric(sd(bin_entropy))
                         bin_entropy_median = as.numeric(median(bin_entropy))
                         bin_entropy_q1 <- as.numeric(quantile(bin_entropy,0.25))
                         bin_entropy_q3 <- as.numeric(quantile(bin_entropy,0.75))
                         
                         # count of total snps found in all query genome at every site in the bin 
                         bin_snp_count <- sum(snps_sum$snp_pos >= bin_start & snps_sum$snp_pos <= bin_end)
                         
                         
                         # check if the bin belongs to any MGE
                         mge_idx <- which(mges_df$start <= bin_end & mges_df$end >= bin_start)
                         
                         return(
                           data.frame(
                             bin_index = bin_idx,
                             coverage_max = max(bin_coverage),
                             coverage_min = min(bin_coverage),
                             coverage_mean = bin_coverage_mean,
                             coverage_sd = bin_coverage_sd,
                             coverage_median = bin_coverage_median,
                             coverage_q1 = bin_coverage_q1,
                             coverage_q3 = bin_coverage_q3,
                             entropy_max = if(length(bin_entropy) == 0) NA else max(bin_entropy),
                             entropy_min = if(length(bin_entropy) == 0) NA else min(bin_entropy),
                             entropy_mean = if(length(bin_entropy) == 0) NA else bin_entropy_mean,
                             entropy_sd = if(length(bin_entropy) == 0) NA else bin_entropy_sd,
                             entropy_median = if(length(bin_entropy) == 0) NA else bin_entropy_median,
                             entropy_q1 = if(length(bin_entropy) == 0) NA else bin_entropy_q1,
                             entropy_q3 = if(length(bin_entropy) == 0) NA else bin_entropy_q3,
                             total_snp_count = bin_snp_count,
                             mge_id = if (length(mge_idx) == 0) "non-mge" else mge_idx
                           )
                         )
                       }))
  
  # Export the summary as csv file
  bin_summary <- merge(chr_bins,
                       bin_stats,
                       by = "bin_index",
                       all = TRUE)
  
  write.csv(bin_summary,
            file.path(output_dir,paste0(input_genome_name,"_bin_summary.csv")),
            quote = FALSE,
            row.names = FALSE)
}

## Visualize the coverage and entropy ----
get_visual <- function(input_genome_name,
                       new_pos,
                       chr_bins,
                       position_coverage,
                       entropy_df,
                       mges_df,
                       output_dir){
  
  # Visualize the coverage and entropy
  ctg_pos_df <- new_pos$ctg_pos_df
  # Coverage plot
  
  coverage_plot_df <- position_coverage %>%
    mutate(cov_pct = 100*coverage/ max(position_coverage$coverage))
  
  coverage_plot <- ggplot() +
    # the rect showing the MGEs
    geom_rect(data = mges_df,
              aes(xmin = start,
                  xmax = end,
                  ymin = -Inf,
                  ymax = Inf),
              fill = "#41b6e6",
              color = "transparent",
              linewidth = 1,
              alpha = 0.65) +
    geom_point(data = coverage_plot_df,
               aes(x = position,
                   y = cov_pct),
               alpha = 0) + 
    geom_line(data = coverage_plot_df,
              aes(x = position,
                  y = cov_pct),
              linewidth = 0.1,
              color = "black",
              alpha = 0.85) + 
    scale_y_continuous(name = "% Coverage") +
    scale_x_continuous(expand = c(0,0)) + 
    theme(
      axis.title.y.left = element_text(colour = "black",
                                       size = 20,
                                       face = "bold"),
      axis.text.y.left = element_text(colour = "black",
                                      size = 15),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      legend.position = "none",
      panel.border = element_rect(colour = "black",
                                  fill = NA,
                                  linewidth = 1)
    )
  
  coverage_plot_margin <- ggMarginal(p = coverage_plot,
                                     type = "histogram",
                                     margins = "y",
                                     groupFill = FALSE,
                                     groupColour = FALSE,
                                     size = 5,
                                     yparams = list(linewidth = 0.15,
                                                    binwidth = 1,
                                                    position = "identity"))
  
  # Entropy Plot
  entropy_plot_df <- entropy_df %>%
    mutate(
      category = sapply(entropy_df$position,
                        function(pos_entry){
                          
                          row_idx <- which(pos_entry >= chr_bins$start &
                                             pos_entry <= chr_bins$end)
                          
                          row_gene <- chr_bins$gene[row_idx]
                          
                          if(all(row_gene %in% "non-cds")){
                            "Intergenic"
                          }else{
                            "Coding"
                          }
                        })
    )
  
  
  # change the negative entropy values to 0
  # not sure why there would be negative entropy values ...
  entropy_plot_df$entropy[entropy_plot_df$entropy < 0] <- 0
  
  entropy_plot <- ggplot() +
    # the rect showing the contig segments
    geom_rect(data = ctg_pos_df %>%
                mutate(
                  x_start = lag(new_end, default = 0),
                  x_end = new_end,
                  fill_color = ifelse(row_number() %% 2 == 1, "#808285", "white")
                ),
              aes(xmin = x_start,
                  xmax = x_end,
                  ymin = -Inf,
                  ymax = 0,
                  fill = fill_color),
              color = "transparent",
              alpha = 0.75) +
    scale_fill_identity() +
    ggnewscale::new_scale_fill() + 
    geom_col(data = entropy_plot_df,
             aes(x = position,
                 y = entropy,
                 fill = category),
             position = "identity",
             width = 10,
             alpha = 0.75) +
    geom_point(data = entropy_plot_df,
               aes(x = position,
                   y = entropy,
                   color = category),
               size = 0.1,
               alpha = 0,
               show.legend = FALSE) +
    # the vertical solid lines indicate the edges of contig concatenation
    geom_segment(data = ctg_pos_df,
                 aes(x = new_end,
                     xend = new_end,
                     y = 0,
                     yend = -Inf), 
                 linetype = "solid",
                 color = "transparent",
                 linewidth = 0.5,
                 alpha = 0.75) +
    scale_color_manual(values = c("Intergenic" = "#E2A4C6", "Coding" = "#91a01e")) +
    scale_fill_manual(values = c("Intergenic" = "#E2A4C6", "Coding" = "#91a01e")) +
    labs(title = paste0(),
         x = "Concatenated Contig Position (bp)",
         y = "Entropy",
         fill = "Position Type") +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous() + 
    theme(axis.text.x = element_text(size = 15),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.text.y = element_text(angle = 0,
                                     size = 15),
          axis.title.y = element_text(size = 20,
                                      face = "bold"),
          axis.title.x = element_text(size = 20,
                                      face = "bold"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.size = unit(0.5,"cm"),
          legend.text = element_text(size = 7.5),
          legend.title = element_text(size = 10),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      linewidth = 1),
          legend.position = "inside",
          legend.position.inside = c(0.95,0.9),
          legend.background = element_rect(colour = NA,
                                           fill = NA)
    )
  
  entropy_plot_margin <- ggMarginal(p = entropy_plot,
                                    type = "histogram",
                                    margins = "y",
                                    groupFill = TRUE,
                                    groupColour = TRUE,
                                    size = 5,
                                    yparams = list(linewidth = 0.15,
                                                   binwidth = max(entropy_plot_df$entropy)/200,
                                                   position = "identity"))
  
  combine_plot <- plot_grid(
    coverage_plot_margin,
    entropy_plot_margin,
    ncol = 1,
    align = "v",
    axis = "tblr",
    rel_heights = c(1, 3)
  )
  
  pdf(file=file.path(output_dir,paste0(input_genome_name,"_profiler_plot.pdf")),
      width=15,
      height=7.5)
  print(combine_plot)
  dev.off()
  
  
  # Also put the 2 data frames in a list and save as RDS just in case later modifications
  plot_df <- list(
    coverage_plot_df = coverage_plot_df,
    entropy_plot_df = entropy_plot_df
  )
  
  saveRDS(plot_df,
          file.path(output_dir,paste0(input_genome_name,"_plot_df.RDS")))
  
  # clean up the Rplots.pdf
  # the way to get rid of this file is too tedious so I just do this
  file.remove("Rplots.pdf")
}

## Main Function ----
profiler <- function(input_genome_path,
                     output_dir, 
                     gff_path,
                     snp_dir,
                     ncores){
  
  setDTthreads(ncores)
  setwd(dir = output_dir)
  
  ## Concatenate the reference genome into 1 if not complete genome ----
  
  if(grepl("GCA_",input_genome_path)){
    # If GenBank genome
    input_genome_name <- paste0("GCA_",strsplit(gsub("\\.fna","",basename(input_genome_path)),
                                                split = "\\_")[[1]][2])
  }else{
    input_genome_name <- gsub("\\.fna|\\.fasta","",basename(input_genome_path))
  }
  
  # Use the helper function to generate new_pos if new_pos does not exist
  
  if(file.exists(file.path(output_dir,paste0(input_genome_name,"_new_pos.RDS")))){
    cat("new_pos Found. Load RDS","\n")
    new_pos <- readRDS(file.path(output_dir,paste0(input_genome_name,"_new_pos.RDS")))
    
  }else{
    cat("new_pos Not Found. Use helper function to generate new_pos","\n")
    new_pos <- concat_genome(input_genome_name = input_genome_name,
                             input_genome_path = input_genome_path,
                             output_dir = output_dir)
    
    cat("Finished generating new_pos","\n")
  }
  
  ## Use the helper function to generate position_coverage if position_coverage does not exist ----
  
  if(file.exists(file.path(output_dir,paste0(input_genome_name,"_position_coverage.RDS")))){
    cat("position_coverage Found. Load RDS","\n")
    position_coverage <- readRDS(file.path(output_dir,paste0(input_genome_name,"_position_coverage.RDS")))
    
  }else{
    cat("position_coverage Not Found. Use helper function to generate position_coverage","\n")
    position_coverage <- get_position_coverage(input_genome_name = input_genome_name,
                                               new_pos = new_pos,
                                               snp_dir = snp_dir,
                                               output_dir = output_dir,
                                               ncores = ncores)
    cat("Finished generating position_coverage","\n")
  }
  
  ## Use the helper function to generate snps_sum if snps_sum does not exist ----
  
  if(file.exists(file.path(output_dir,paste0(input_genome_name,"_snps_sum.RDS")))){
    cat("snps_sum Found. Load RDS","\n")
    snps_sum <- readRDS(file.path(output_dir,paste0(input_genome_name,"_snps_sum.RDS")))
    
  }else{
    cat("snps_sum Not Found. Use helper function to generate snps_sum","\n")
    snps_sum <- get_snps_sum(input_genome_name = input_genome_name,
                             new_pos = new_pos,
                             snp_dir = snp_dir,
                             output_dir = output_dir,
                             ncores = ncores)
    cat("Finished generating snps_sum","\n")
  }
  
  ## Use the helper function to generate chr_bins if chr_bins does not exist ----
  
  if(file.exists(file.path(output_dir,paste0(input_genome_name,"_chr_bins.RDS")))){
    cat("chr_bins Found. Load RDS","\n")
    chr_bins <- readRDS(file.path(output_dir,paste0(input_genome_name,"_chr_bins.RDS")))
    
  }else{
    cat("chr_bins Not Found. Use helper function to generate chr_bins","\n")
    chr_bins <- get_chr_bins(input_genome_name = input_genome_name,
                             new_pos = new_pos,
                             gff_path = gff_path,
                             output_dir = output_dir)
    cat("Finished generating chr_bins","\n")
  }
  
  ## Use the helper function to generate entropy_df if entropy_df does not exist ----
  
  if(file.exists(file.path(output_dir,paste0(input_genome_name,"_entropy.RDS")))){
    cat("entropy_df Found. Load RDS","\n")
    entropy_df <- readRDS(file.path(output_dir,paste0(input_genome_name,"_entropy.RDS")))
    
  }else{
    cat("entropy_df Not Found. Use helper function to generate entropy_df","\n")
    entropy_df <- get_entropy(input_genome_name = input_genome_name,
                               new_pos = new_pos,
                               snps_sum = snps_sum,
                               position_coverage = position_coverage,
                               output_dir = output_dir,
                               ncores = ncores)
    cat("Finished generating entropy_df","\n")
  }
  ## Use the helper function to infer MGEs -----
  
  if(file.exists(file.path(output_dir,paste0(input_genome_name,"_mges.RDS")))){
    cat("mges_df Found. Load RDS","\n")
    
    mges_df <- readRDS(file.path(output_dir,paste0(input_genome_name,"_mges.RDS")))
    
  }else{
    cat("mges_df Not Found. Use helper function to generate mges_df.RDS","\n")
    mges_df <- get_mges(input_genome_name = input_genome_name,
                        new_pos = new_pos,
                        entropy_df = entropy_df, 
                        position_coverage = position_coverage,
                        chr_bins = chr_bins,
                        min_mge_size = 100,
                        cov_cutoff = 0.7,
                        output_dir = output_dir,
                        ncores = ncores)
    cat("Finished generating mges_df","\n")
  }
  
  ## Use the helper function to summarize entropy and coverage ---- 
  
  if(!file.exists(file.path(output_dir,paste0(input_genome_name,"_bin_summary.csv")))){
    cat("bin_summary.csv Not Found. Use helper function to generate bin_summary.csv","\n")
    get_bin_sum(input_genome_name = input_genome_name,
                chr_bins = chr_bins,
                position_coverage = position_coverage,
                snps_sum = snps_sum,
                mges_df = mges_df,
                entropy_df = entropy_df,
                output_dir = output_dir)
    cat("Finished generating bin_summary.csv","\n")
  }
  
  ## Use the helper function to visualize the coverage and entropy if not exist ---- 
  
  if(!file.exists(file.path(output_dir,paste0(input_genome_name,"_profiler_plot.pdf")))){
    cat("Visualization Not Found. Use helper function to create visualization","\n")
    get_visual(input_genome_name = input_genome_name,
               new_pos = new_pos,
               chr_bins = chr_bins,
               position_coverage = position_coverage,
               entropy_df = entropy_df,
               mges_df = mges_df,
               output_dir = output_dir)
    cat("Finished Visualization","\n")
  }
  
}

# Execute the function----
profiler(input_genome_path = input_genome_path,
         output_dir = output_dir, 
         gff_path = gff_path,
         snp_dir = snp_dir,
         ncores = ncores)

