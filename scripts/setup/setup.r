library(automateR)  # custom package

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
load(file.path(scripts.dir, "rna-seq/dge.RData"))
load(file.path(scripts.dir, "chip-seq/chip-dba.RData"))
load(file.path(scripts.dir, "atac-seq/atac-dba.RData"))
