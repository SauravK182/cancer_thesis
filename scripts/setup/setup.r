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
tryCatch(
    {
        setwd(main.data.dir)
    },
    error = function(e) {
        message("Data directory does not exist/cannot be found. Moving on...")
    }
)

# Source necessary files
load(file.path(scripts.dir, "rna-seq/dge.RData"))
load(file.path(scripts.dir, "chip-seq/chip-dba.RData"))
load(file.path(scripts.dir, "atac-seq/atac-dba.RData"))
