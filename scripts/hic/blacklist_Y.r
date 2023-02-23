#!/bin/usr/env Rscript

# Grab arguments from commandline
args <- commandArgs(trailingOnly = TRUE)

# Read in files, convert chromosome names to appropriate type, blacklist chrY
# Per dcHiC, blacklisted regions = 1
for (file in args) {
    bedfile <- read.delim(file, header = FALSE)
    bedfile$V1 <- paste0("chr", bedfile$V1) # convert to UCSC genome name
    bedfile$blacklist <- ifelse(bedfile$V1 == "chrY", yes = 1, no = 0)

    # Write to new .bed file
    filename <- gsub(pattern = ".bed", x = file, replacement = "_blacklisted.bed")
    write.table(bedfile, filename, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}