# Load necessary packages
library(tidyverse)
library(DESeq2)
library(EnsDb.Hsapiens.v79)

# Set up data directories
main.data.dir <- "D:/SK/data"
ren.data.dir <- "ren-panc"
rod.data.dir <- "rod-rcc"
cai.data.dir <- "cai-bc"
rna.data.dir <- "rna-seq/aligned"

# Scripts directory
scripts.dir <- "C:/Users/jvons/Documents/NCF/Thesis/scripts"
rna.scripts.dir <- "rna-seq"

# Set directory to main data directory
setwd(main.data.dir)

# Source necessary files
source(file.path(scripts.dir, rna.scripts.dir, "dge_analysis.r"))
