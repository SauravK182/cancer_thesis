library(automateR)  # custom package
library(cowplot)
library(ggsignif)
library(viridis)
library(memes)

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

# Colors for plotting
n <- length(dge.list.full)
col.vec <- viridis_pal(option = "C")(n)
names(col.vec) <- names(dge.list.full)
# pie(rep(1, n), col = col.vec)

# Names
names.comp <- c("Pancreatic System",
                "786 ccRCC System",
                "OS ccRCC System",
                "BrM2 Brain vs. Primary",
                "LM2 Lung vs. Primary")
