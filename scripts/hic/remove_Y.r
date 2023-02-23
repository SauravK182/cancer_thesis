#!/bin/usr/env Rscript
library(optparse)
library(readr)

# See https://stackoverflow.com/questions/2151212/ for more

# Make options
option_list <- list(
    make_option(c("-m", "--matrix"), type = "character",
        help = "Matrix file from HiC-Pro for Y chromosome interactions to be removed."),
    make_option(c("-b", "--bed"), type = "character",
        help = "BED file associated with the HiC-Pro matrix that indicates chr and location of interactions.")
)
options(scipen = 999) # to turn off scientific notation

# Get options
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opts <- args$options

# Read matrix and bedfile
hic.mat <- as.data.frame(read_table(opts$matrix, col_names = FALSE))
hic.bed <- as.data.frame(read_table(opts$bed, col_names = FALSE))

# Find all genome indices corresponding to chrY
chr.y.indices <- ifelse(test = hic.bed$X1 == "Y" | hic.bed$X1 == "chrY", yes = TRUE, no = FALSE)
chr.y <- hic.bed[chr.y.indices, 4]

# Filter Hi-C matrix and bed file
hic.mat.y <- hic.mat$X1 %in% chr.y.indices | hic.mat$X2 %in% chr.y.indices
hic.mat <- hic.mat[-hic.mat.y, ]
hic.bed <- hic.bed[-chr.y.indices, ]

# Write matrix and bed files
mat.filename <- gsub(pattern = ".matrix$", replacement = "_filtered.matrix", x = opts$matrix)
bed.filename <- gsub(pattern = ".bed", replacement = "_filtered.bed", x = opts$bed)

write.table(hic.mat, mat.filename, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(hic.bed, bed.filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)