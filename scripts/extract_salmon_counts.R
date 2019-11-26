# ========== EXECUTION =======================
# R CMD BATCH --no-save --no-restore \
# '--args dir="/project/influenza/lamha/salmon/" feat="genes"' \
# extract_salmon_counts.R
# ============================================


# ========== DESCRIPTION =====================
# parses salmon files (requires "quant.genes.sf" files)
# get qc information        -> output: sample_info.csv
# get tpm matrix w/o ERCCs  -> output: tpm.csv
# get tpm matrix with ERCCs -> output: all_tpm.csv
# =============================================

args <- (commandArgs(TRUE))
if (length(args) == 0) {
  print("No arguments supplied.")
} else {
  for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
options(echo=TRUE) # if you want see commands in output file


source("/project/influenza/Scripts/read_quant_files.R")


sample_info <- read_qc_info_salmon(dir)
tpm   <- read_quants_salmon(dir, feat)
count <- read_quants_salmon(dir, feat, counts="NumReads")

### read annotation
gene_annotation <- read.csv("/project/influenza/Annotation/Canis_familiaris_CanFam3.1_E86_cDNA_influenza_H1N1_gene_annotation.csv",
                            check.names = F)

### remove ERCC and rescale to TPM

gene_tpm <- tpm[,grep("ERCC-", colnames(tpm), invert=T)]
gene_tpm <- gene_tpm / rowSums(gene_tpm) * 1e6

### get some statistics

sample_info$total_ERCC <- rowSums(tpm[,grep("ERCC-", colnames(tpm))])
sample_info$total_cell <- rowSums(gene_tpm[,grep("H1N1", colnames(gene_tpm), invert = T)])
sample_info$detected_cell_genes <- rowSums(gene_tpm[,grep("H1N1", colnames(gene_tpm), invert=T)] > 1)
sample_info$num_detected_h1n1_transcripts <- rowSums(gene_tpm[,grep("H1N1", colnames(gene_tpm))] > 1)
sample_info$total_H1N1 <- rowSums(gene_tpm[,grep("H1N1", colnames(gene_tpm))])

# get mitochondrial info
idx    <- match(colnames(gene_tpm), gene_annotation$`Ensembl Gene ID`)
mt_idx <- which(gene_annotation[idx,"Chromosome Name"] == "MT")
sample_info$MT_content <- rowSums(gene_tpm[,mt_idx])


gene_tpm <- as.data.frame(t(gene_tpm))
tpm      <- as.data.frame(t(tpm))
count    <- as.data.frame(t(count))


### write output

write.table(sample_info,
            file = "sample_info.csv",
            quote = FALSE,
            sep = ",",
            col.names = T,
            row.names = T)

write.table(gene_tpm,
            file = paste0("tpm_", feat, ".csv"),
            quote = FALSE,
            sep = ",",
            col.names = T,
            row.names = T)

write.table(tpm,
            file = paste0("all_tpm_", feat, ".csv"),
            quote = FALSE,
            sep = ",",
            col.names = T,
            row.names = T)

write.table(count,
            file = "count.csv",
            quote = FALSE,
            sep = ",",
            col.names = T,
            row.names = T)
