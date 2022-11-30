library(tidyverse)
library(org.Hs.eg.db)
library(annotate)
library(DESeq2)
library(gprofiler2)
library(ComplexHeatmap)

# Set directory to data in D: drive
project.dir <- "D:/SK/data"
ren.dir <- "ren-panc"
rna.dir <- "rna-seq-aligned"

setwd(file.path(project.dir, ren.dir, rna.dir))

# Read in featureCounts file
## Skip first line, which contains the command used
featurecounts <- read.delim("ren_rna_rawcounts.txt", header = TRUE, skip = 1)

# First column to the right of length is where the counts start
## Counts go to end of matrix - extract these cols only!
countcol.start <- which(colnames(featurecounts) == "Length") + 1
countcol.end <- ncol(featurecounts)

# Isolate actual counts and rename rows to the Ensembl IDs
counts.df <- featurecounts[, countcol.start:countcol.end]
rownames(counts.df) <- featurecounts[, "Geneid"]

# Remove the ".bam" extension after each column in the counts.df matrix
colnames(counts.df) <- gsub(colnames(counts.df), pattern = ".bam", replacement = "")

# Set up coldata df, which tells DESeq info about each column of the count matrix
## Will also serve as the "design" - how to model the samples, estimate dispersions and log2-fold changes b/w samples
coldata <- data.frame(Culture = rep(c("CAPAN1", "HPNE", "PANC1"), each = 3))

# Set factor level so DESeq2 knows which is the reference (lowest level)
coldata$Culture <- factor(coldata$Culture, levels = c("CAPAN1", "PANC1", "HPNE"))

# Run DESeq to get the DESeqDataSet, keep only features with > 10 counts
ren.dds <- DESeqDataSetFromMatrix(countData = counts.df,
                                  colData = coldata,
                                  design = ~ Culture)
keepCount <- rowSums(counts(ren.dds)) >= 10
ren.dds <- ren.dds[keepCount, ]
ren.dge <- DESeq(ren.dds)

# Save size factors
ren.sf <- sizeFactors(ren.dge)
ren.sf
summary(ren.sf)

# The size factors range from 0.9592 - 1.08 - very good
## The fact that size factors are clustered indicates the sequencing of the libraries
## Was similiar in depth - our comparisons will be stronger
## Normalization technique:
### Step 1: for each gene, calculate the row-wise geometric means (i.e., geometric mean of its expression lvl across all samples)
### Step 2: For each gene in all samples, divide the raw count by its corresponding geometric mean
### Step 3: For each sample, take the median of normalized gene counts - these are the dfs
### Since the DESeq2 assumption is that most genes are not DE, the median normalized expression
### Should be around 1 by contruction

# Pairwise comparison of all 3 cases
## Note that DESeq2's sample_1 vs. sample_2 indicates the LFC values are log2(sample_1/sample_2)
contrasts <- list(c("PANC1", "HPNE"), c("CAPAN1", "HPNE"), c("CAPAN1", "PANC1"))
alphaTest <- 0.05
ren.hpne.panc <- results(ren.dge, contrast = c("Culture", contrasts[[1]]), alpha = alphaTest)
ren.hpne.capan <- results(ren.dge, contrast = c("Culture", contrasts[[2]]), alpha = alphaTest)
ren.panc.capan <- results(ren.dge, contrast = c("Culture", contrasts[[3]]), alpha = alphaTest)

# Plot PCA
rld <- rlog(ren.dge, blind = TRUE)
plotPCA(rld, intgroup = "Culture")

# Plot dispersion
plotDispEsts(ren.dge)
