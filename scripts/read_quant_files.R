
library(data.table)
library(reshape2)


# =======================================
# salmon
# =======================================

### get log information
read_qc_info_salmon <- function(dir) {
  
  pattern <- "meta_info.json"
  salmon_qc_files <- list.files(dir, pattern=pattern, recursive = T, full.names = T)
  
  qc_list <- lapply(salmon_qc_files, function(file) {
    
    dat <- readLines(file)
    dat <- gsub(",$", "", dat)
    
    log_df <- colsplit(dat, ": ", c("info", "value"))
    
    num_processed  <- as.numeric(log_df$value[grep("num_processed", log_df$info)])
    num_mapped     <- as.numeric(log_df$value[grep("num_mapped", log_df$info)])
    percent_mapped <- as.numeric(log_df$value[grep("percent_mapped", log_df$info)])
    
    sample_name <- gsub("/aux_info/meta_info.json", "", file)
    
    log_df <- data.frame(num_processed = num_processed,
                         num_mapped = num_mapped,
                         percent_mapped = percent_mapped,
                         row.names = sample_name)
    
    return(log_df)
  })
  
  sample_info <- do.call(rbind, qc_list)
  sample_info <- as.data.frame(sample_info)
  
  return(sample_info)
  
}



### get quantification values
read_quants_salmon <- function(dir, feat="genes", counts="TPM") {
  
  if (feat == "genes") {
    pattern <- "quant.genes.sf"
  } else if (feat == "transcripts") {
    pattern <- "quant.sf"
  } else {
    stop("feat needs to be genes or transcripts")
  }
  
  salmon_quant_files <- list.files(dir, pattern=pattern, recursive = T, full.names = T)
  
  count_list <- lapply(salmon_quant_files, function(file) {
    
    dat <- fread(file, data.table = F)
    if (counts == "TPM") {
      dat <- subset(dat, select = c(Name, TPM))
    } else {
      dat <- subset(dat, select = c(Name, NumReads))
    }
    
    rownames(dat) <- dat$Name
    dat$Name <- NULL
    
    if (!nrow(dat)) { dat <- NA }
    
    return(dat)
  })
  
  expr <- do.call(cbind, count_list)
  # expr <- expr[,! colSums(is.na(expr)) == nrow(expr)]
  colnames(expr) <- gsub(paste0("/", pattern), "", salmon_quant_files)
  expr <- as.data.frame(t(expr))
  expr <- expr[,order(colnames(expr))]
  
  return(expr)
  
}


# =======================================
# featureCounts
# =======================================


### get log information
read_qc_info_featureCounts <- function(dir) {
  
  log_ext <- ".counts.summary$"
  qc_files <- list.files(dir, pattern = log_ext, full.names = T)
  
  qc_list <- lapply(qc_files, function(file) {
    
    sample <- gsub(log_ext, "", basename(file))
    
    log_df <- fread(file, data.table = F)
    colnames(log_df)[ncol(log_df)] <- sample
    
    num_processed  <- sum(log_df[,sample])
    num_mapped     <- log_df[log_df$Status == "Assigned", sample]
    percent_mapped <- num_mapped / num_processed * 100
    
    log_df <- data.frame(num_processed = num_processed,
                         num_mapped = num_mapped,
                         percent_mapped = percent_mapped,
                         row.names = sample)
    
    return(log_df)
  })
  
  sample_info <- do.call(rbind, qc_list)
  sample_info <- as.data.frame(sample_info)
  
  return(sample_info)
  
}



### get quantification values
read_quants_featureCounts <- function(dir, mode="fpkm") {
  
  count_ext <- ".counts$"
  count_files <- list.files(dir, pattern = count_ext, full.names = T)
  
  count_list <- lapply(count_files, function(file) {
    
    sample <- gsub(count_ext, "", basename(file))
    
    dat <- fread(file, data.table = F)
    colnames(dat)[ncol(dat)] <- sample
    rownames(dat)    <- dat$Geneid
    
    if (mode != "raw") {
      if (mode == "fpkm"){
        # scale to FPKM
        dat[,sample] <- dat[,sample] / dat[,"Length"]
      } else if (mode == "tpm") {
        # scale to TPM
        dat[,sample] <- dat[,sample] / sum(dat[,sample]) * 1e6
      } 
    }
    dat <- subset(dat, select=sample)
    
    return(dat)
  })
  
  all_count <- do.call(cbind, count_list)
  all_count <- as.data.frame(t(all_count))
  all_count <- all_count[,order(colnames(all_count))]
  
  return(all_count)
  
}


# =======================================
# star
# =======================================

### get log information
read_qc_info_star <- function(dir) {
  
  pattern <- "Log.final.out"
  qc_files <- list.files(dir, pattern=pattern, recursive = T, full.names = T)
  
  qc_list <- lapply(qc_files, function(file) {
    
    dat <- readLines(file)
    dat <- trimws(dat)
    
    log_df <- colsplit(dat, " \\|\t", c("info", "value"))
    
    num_processed  <- as.numeric(log_df$value[grep("Number of input reads", log_df$info)])
    num_mapped     <- as.numeric(log_df$value[grep("Uniquely mapped reads number", log_df$info)])
    percent_mapped <- as.numeric(sub("%", "", log_df$value[grep("Uniquely mapped reads %", log_df$info)]))
    num_chimeric_reads <- as.numeric(log_df$value[grep("Number of chimeric reads", log_df$info)])
    
    sample_name <- gsub("_star_out_Log.final.out", "", basename(file))
    
    log_df <- data.frame(num_processed = num_processed,
                         num_mapped = num_mapped,
                         percent_mapped = percent_mapped,
                         num_chimeric_reads = num_chimeric_reads,
                         row.names = sample_name)
    
    return(log_df)
  })
  
  sample_info <- do.call(rbind, qc_list)
  sample_info <- as.data.frame(sample_info)
  
  return(sample_info)
  
}


