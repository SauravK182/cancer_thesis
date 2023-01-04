library(strawr)
library(tidyverse)
library(HiCcompare)

# Load dataset directory
myData <- file.path("D:/SK/data/ren-panc/hic")
capan.1.test <- file.path(myData, "Capan-1/Capan-1.chr1.hic")

# Use straw() to read in the test file with KR normalization
capan.df <- straw(fname = capan.1.test,
                  norm = "KR",
                  chr1loc = "1",
                  chr2loc = "1",
                  unit = "BP",
                  binsize = 5000)

# Read Capan-1 Hi-C files into R
## Pre-create master list, which will be a list of 3 lists
### Each list will pertain to one cell line, and contain all 23 .hic matrices
### Each of the 3 cell-line lists will have key names chr<#>

# Note that cellLines is ordered normal | primary tumor | metastatic tumor
cellLines <- c("HPNE", "PANC-1", "Capan-1")
chromNumber <- 23
hicListCulture.keys <- c(paste0("chr", 1:22), "chrX")

## Make lists for each of the three conditions, create the keys for each
hpne <- vector(mode = "list", length = chromNumber)
panc <- vector(mode = "list", length = chromNumber)
capan <- vector(mode = "list", length = chromNumber)

## Make the master list and re-name the keys of each sub-list
hicList <- list(hpne, panc, capan)
names(hicList) <- cellLines
for (i in 1:length(cellLines)) {
    names(hicList[[i]]) <- hicListCulture.keys
}

## Iterate through Hi-C matrices in each of the Supp Files, read each .hic
## Note that files are stored as, e.g., /HPNE/HPNE.chr1.hic
### I.e., files are stored as /<cell_line>/<cell_line>.chr<#>.hic
for (chrm in c(1:22, "X")) {
    for (culture in cellLines) {
        hic.df <- straw(fname = paste("./", culture, "/", culture, ".chr", chrm, ".hic", sep = ""),
                    norm = "NONE",
                    chr1loc = chrm,
                    chr2loc = chrm,
                    unit = "BP",
                    binsize = 5000)
    chrom <- paste0("chr", i)
    hicList[[culture]][[chrom]] <- hic.df
    rm(hic.df)
    }
}

# Convert sparse (i.e. N x 3) upper triangular matrices to a hic.table from HiCcompare