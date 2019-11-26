# ========== EXECUTION =======================
# R CMD BATCH --no-save --no-restore \
# '--args dir="/project/influenza/lamha/star/aligned_to_genome/H1N1_only/" splice_type="canonical"' \
# extract_chimeric_reads.R
# ============================================


# ========== DESCRIPTION =====================

# =============================================

library(data.table)
library(reshape2)
library(ggplot2)


args <- (commandArgs(TRUE))
if (length(args) == 0) {
  print("No arguments supplied.")
} else {
  for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
options(echo=TRUE) # if you want see commands in output file


read_chimJunctions_star <- function(dir, splice_type) {
  
  if (! splice_type %in% c("canonical", "all")) {
    stop("unknown splice_type, splice_type must be either canonical or all")
  } 
  
  pattern <- "Chimeric.out.junction"
  star_chim_files <- list.files(dir, pattern=pattern, recursive = T, full.names = T)
  
  colNames  <- c("chr_donor",
                 "start_donor",
                 "strand_donor",
                 "chr_acceptor",
                 "start_acceptor",
                 "strand_acceptor",
                 "junction_type",
                 "repeat_length_left",
                 "repeat_length_right",
                 "read_name",
                 "first_bp_first_segment",
                 "CIGAR_first_segment",
                 "first_bp_second_segment",
                 "CIGAR_second_segment")
  
  chim_dat <- lapply(star_chim_files, function(file) {
    
    info <-  file.info(file)
    
    if (info$size) {
      dat <- fread(file)
      colnames(dat) <- colNames
      dat$sample <- gsub("_h1n1_star_out_Chimeric.out.junction", "", basename(file))
      
      return(dat)
    } else {
      return(NULL)
    }
  })
  
  chim_dat <- do.call(rbind, chim_dat)
  
  if (splice_type == "canonical") {
    chim_dat <- subset(chim_dat, junction_type > 0)
  }
  
  return(chim_dat)
  
}

#### execute #### 

dir <- "/project/influenza/lamha/star/aligned_to_genome/H1N1_only/"
dat2 <- read_chimJunctions_star(dir, "canonical")


