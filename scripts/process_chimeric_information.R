# ====================================
# calculate distance of chimeric reads and append to data frame
# ====================================
get_chim_read_distance <- function(chim_dat) {
  
  chim_dat$distance <- ifelse(chim_dat$chr_don == chim_dat$chr_acc &
                                chim_dat$chr_don_str == chim_dat$chr_acc_str, 
                              ifelse(chim_dat$chr_don_str == "+",
                                     chim_dat$chr_acc_start - chim_dat$chr_don_start + 1,
                                     chim_dat$chr_don_start - chim_dat$chr_acc_start + 1),
                              NA)
  return(chim_dat)
}

# ====================================
# parse chimeric.out.junction files
# ====================================
get_chim_info <- function(dir, sample_names = F) {
  
  pat <- "*_R1_R2_H1N1_Chimeric.out.junction"
  chim_files <- list.files(dir, pattern = pat, full.names = T)
  
  dat_numReads <- lapply(chim_files, function(file) {
    
    finfo <- file.info(file)
    empty <- finfo$size == 0
    
    if (!empty) {
      dat <- fread(file)
      dat <- subset(dat, grepl("H1N1", dat$V1) & grepl("H1N1", dat$V4))
      
      if(nrow(dat)) {
        colnames(dat)[1:6] <- c("chr_don",
                                "chr_don_start",
                                "chr_don_str",
                                "chr_acc",
                                "chr_acc_start",
                                "chr_acc_str")
        colnames(dat)[10] <- "rname"
        dat_numReads <- aggregate(rname ~ chr_don +
                                    chr_don_start +
                                    chr_don_str + 
                                    chr_acc +
                                    chr_acc_start + 
                                    chr_acc_str, data = dat, FUN = length)
        colnames(dat_numReads)[7] <- "num_reads"
        
        dat_numReads$sample <- gsub("_R1", "", dat$V15[1])
        
        return(dat_numReads)
      }
    }
    
    return(NULL)
    
    
  })
  
  dat_numReads <- do.call(rbind, dat_numReads)
  
  if(!sample_names) {
    chim_num_samples <- aggregate(sample ~ chr_don +
                                    chr_don_start +
                                    chr_don_str + 
                                    chr_acc +
                                    chr_acc_start + 
                                    chr_acc_str, dat_numReads, FUN = length)
    
    chim_sum_reads <- aggregate(num_reads ~ chr_don +
                                  chr_don_start +
                                  chr_don_str + 
                                  chr_acc +
                                  chr_acc_start + 
                                  chr_acc_str, dat_numReads, FUN = sum)
    chim_stat <- merge(chim_num_samples, chim_sum_reads)
  } else {
    chim_stat <- dat_numReads
  }
  
  chim_stat <- get_chim_read_distance(chim_stat)
  
  # filter out chimeric reads spanning over different H1N1 segments
  # and negative distal junction sites
  chim_stat <- subset(chim_stat, subset= !is.na(distance) & distance > 0)
  
  return(chim_stat)
}



# ====================================
# get read counts per junction and sample
# ====================================
get_junction_counts <- function(dir, junctions) {
  
  # junctions: data.frame with chr, start, end
  
  if(grep("junction_merged", colnames(junctions))) {
    merge_j <- T
  } else { merge_j <- F }
  
  pat <- "*_R1_R2_H1N1_Chimeric.out.junction"
  chim_files <- list.files(dir, pattern = pat, full.names = T)
  
  dat_count_reads <- lapply(chim_files, function(file) {
    
    finfo <- file.info(file)
    empty <- finfo$size == 0
    
    if (!empty) {
      dat <- fread(file)
      dat <- subset(dat, grepl("H1N1", dat$V1) & grepl("H1N1", dat$V4))
      
      if(nrow(dat)) {
        colnames(dat)[1:6] <- c("chr_don",
                                "chr_don_start",
                                "chr_don_str",
                                "chr_acc",
                                "chr_acc_start",
                                "chr_acc_str")
        colnames(dat)[10] <- "rname"
        
        res <- merge(dat, junctions,
                     by.x = c("chr_don", "chr_don_start", "chr_acc_start"),
                     by.y = c("chr", "start", "end"))
        
        # make sure chimeric read is on same segment
        res <- subset(res, chr_don == chr_acc)
        
        if(nrow(res)) {
          # read counts per junction
          counts <- aggregate(rname ~ junction_merged + chr_don,
                              data = res, FUN = length)
          
          counts$sample <- gsub("_R.$", "", res$V15)[1]
          
          setnames(counts, "rname", "num_reads")
          
          return(counts)
        }
      } 
    }
    return(NULL)
  })
  
  dat_count_reads <- do.call(rbind, dat_count_reads)
  
  dat_count_reads$junctions <- paste0(dat_count_reads$chr_don,
                                      ":",
                                      dat_count_reads$junction_merged)
  
  
  
  dat_count_reads <- dcast(dat_count_reads, sample ~ junctions,
                           value.var = "num_reads", fill = 0)
  
  return(dat_count_reads)
}


# ====================================
# merge chimeric junctions of forward and reverse strands
# ====================================
merge_loci <- function(df, window=3) {
  
  df_by_segments <- split(df, list(df$chr_don, df$chr_acc), drop=TRUE)
  
  merged_loci_list <- lapply(df_by_segments, function(ddf) {
    
    ddf <- ddf[order(ddf$num_reads, decreasing = T),]
    
    df_edit <- data.frame()
    
    while (nrow(ddf) >= 1) {
      
      start <- as.numeric(ddf$chr_don_start[1])
      end   <- as.numeric(ddf$chr_acc_start[1]) 
      str   <- ddf$chr_don_str[1]
      
      curr_row <- ddf[1,]
      
      ddf <- ddf[-1,]
      
      
      if (!nrow(ddf)) { # last row of ddf being deleted
        df_edit <- rbind(df_edit, curr_row)
        break
      }
      
      if(str == "+") {
        
        ### get subset which has only complementary strand and matching positions
        
        comp_ddf <- subset(ddf, subset =
                             chr_acc_str   == "-" &
                             chr_acc_start <= start + window &
                             chr_don_start >= end   - window)
        
      } else {
        
        comp_ddf <- subset(ddf, subset =
                             chr_acc_str   == "+" &
                             chr_acc_start >= start - window &
                             chr_don_start <= end   + window) 
      }
      
      if (nrow(comp_ddf)) {
        
        if(nrow(comp_ddf) > 1) {
          # if there are multiple entries take entry
          # with the minimum distance to current row entry
          min_dist <- min(abs(comp_ddf$chr_acc_start - start))
          comp_ddf <- subset(comp_ddf, subset =
                               chr_acc_start == start + min_dist | 
                               chr_acc_start == start - min_dist )
        }
        
        # found matching entry -> merge both rows
        
        final_row <- curr_row
        
        
        
        final_row["chr_don_start"] <- paste(sort(unique(c(curr_row$chr_don_start,
                                                          comp_ddf$chr_don_start))),
                                            collapse = "/", sep="/")
        final_row["chr_acc_start"] <- paste(sort(unique(c(curr_row$chr_acc_start,
                                                          comp_ddf$chr_acc_start)), decreasing = T),
                                            collapse = "/", sep="/")
        final_row["chr_don_str"]   <- "+/-"
        final_row["chr_acc_str"]   <- "+/-"
        final_row["num_reads"]     <- sum(final_row$num_reads, comp_ddf$num_reads)
        
        df_edit <- rbind(df_edit, final_row)
        
        # delete entry in ddf
        ddf <- anti_join(ddf, comp_ddf,
                         by = c("chr_don",
                                "chr_don_start",
                                "chr_don_str",
                                "chr_acc",
                                
                                "chr_acc_start",
                                "chr_acc_str",
                                "sample",
                                "num_reads",
                                "distance"))
        
      } else {
        
        df_edit <- rbind(df_edit, curr_row)
      }
    }
    return(df_edit)
  }) 
  
  merged_loci_df <- do.call(rbind, merged_loci_list)
  rownames(merged_loci_df) <- NULL
  
  return(merged_loci_df)
}



merge_junctions <- function(chim_dat) {
  
  merge_junctions_df <- lapply(split(chim_dat, chim_dat$distance), merge_loci)
  merge_junctions_df <- do.call(rbind, merge_junctions_df)  
  
  return(merge_junctions_df)
}

# ====================================
# splits annotated junctions 
# ====================================
split_junctions <- function(my_junction) {
  
  # split chr:start_don/start_acc-end_don/end_acc annotation
  my_junction <- colsplit(my_junction, ":", names = c("chr", "junction_merged"))
  
  tmp   <- colsplit(my_junction$junction_merged, "-", names = c("start", "end"))
  start <- colsplit(tmp$start, "/", names=c("start1", "start2"))
  end   <- colsplit(tmp$end, "/", names=c("end1", "end2"))
  
  my_junction <- rbind(
    cbind(my_junction, start = start$start1, end = end$end1),
    cbind(my_junction, start = start$start2, end = end$end2))
  
  # remove NA's
  my_junction <- my_junction[rowSums(is.na(my_junction)) == 0,]
  
  return(my_junction)
}


# ====================================
# compare junctions
# ex: H1N1_seg2:1/100 - 99/2 == H1N1_seg2:1-99 | H1N1_seg2:100-2
# ====================================
compare_junctions <- function(junction1, junction2, keep.x = F) {
  
  junction1_split <- split_junctions(junction1)
  junction2_split <- split_junctions(junction2)
  
  common_junctions <- merge(junction1_split, junction2_split,
                            by = c("chr", "start", "end"),
                            all.x = keep.x)
  
  return(common_junctions)
  
}
