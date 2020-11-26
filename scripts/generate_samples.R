#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[[1]]
files <- list.files(path, full.names = TRUE, recursive = TRUE, pattern = "fastq")
fq1 <- sort(grep(files, pattern = "_R1", value = TRUE))
fq2 <- sort(grep(files, pattern = "_R2", value = TRUE))
sample <- sapply(strsplit(basename(fq1), "_"), function(x) x[1])
sample <- gsub("-[^V]+", "", sample)
run <- sapply(strsplit(basename(fq1), "_R[1,2]"), function(x) x[1])
samples <- data.frame(sample, run, fq1, fq2, platform = rep("ILLUMINA", length(sample)), stringsAsFactors = FALSE)
# rewrite old samples file
old <- list.files(".", pattern = "samples.tsv")
stopifnot(length(old) <= 1)
file.copy(old, paste0("samples_", gsub("\\s+", "_", Sys.time()), ".tsv"))
write.table(samples, "samples.tsv", row.names=FALSE, quote = FALSE)

