txt2counts <- function(featurecounts.txt) {
    # Read file, skip first line which contains terminal command
    featurecounts <- read.delim(featurecounts.txt, header = TRUE, skip = 1)

    # First column to the right of length is where counts start
    countcol.start <- which(colnames(featurecounts) == "Length") + 1
    countcol.end <- ncol(featurecounts)

    # Counts go to end of matrix, extract these cols only
    counts.df <- featurecounts[, countcol.start:countcol.end]
    rownames(counts.df) <- featurecounts[, "Geneid"]    # Sets rownames to Ensembl IDs

    # Remove the ".bam" or "_filtered.bam" extension after each column in the counts.df matrix
    colnames(counts.df) <- gsub(colnames(counts.df), pattern = "_trimmed[.]bam|[.]bam", replacement = "")
    return(counts.df)
}

# Performs DGE analysis on a given featurecounts text file in current WD
# Requires:
    # featurecounts.txt: Name/path to a featurecounts count file
    # num.reps: Number of replicates per sample
    # culture.facts: A character vector of the cultures - NOTE THAT THE ORDER OF THE VECTOR WILL BE THE ORDER OF THE FACTOR LEVEL!
        ## Assumes that the character vector is already leveled such that, e.g., metastasis > primary cancer > control
    # alphaTest: A FDR-adjusted p-value to determine significance for DE. Default is 0.05.
    # lfc: A threshold to pass to results() from DESeq2 in determining the cutoff for what consitutes a differentially expressed gene.
        ## Default is 2.
# Return value: A list of DESeqResults objects
dge_analysis <- function(featurecounts, culture.factors, num.reps, alphaTest = 0.05, lfc = 2) {
    require(DESeq2)
    require(org.Hs.eg.db)

    if (class(featurecounts) == "character") {
        counts.df <- txt2counts(featurecounts)
    } else if (class(featurecounts) != "data.frame") {
        stop("Argument 'featurecounts' must be either a valid filename or data frame.")
    }

    # Set up coldata df, which tells DESeq info about each column of the count matrix
    ## Will also serve as the "design" - how to model the samples, estimate dispersions and log2-fold changes b/w samples
    coldata <- data.frame(Culture = rep(culture.factors, each = num.reps))

    # Set factor level so DESeq2 knows which is the reference (lowest level)
    coldata$Culture <- factor(coldata$Culture, levels = culture.factors)

    # Run DESeq to get the DESeqDataSet, keep only features with > 10 counts
    canc.dds <- DESeqDataSetFromMatrix(countData = counts.df,
                                    colData = coldata,
                                    design = ~ Culture)
    keepCount <- rowSums(counts(canc.dds)) >= 10
    canc.dds <- canc.dds[keepCount, ]
    canc.dge <- DESeq(canc.dds)

    # Save size factors
    canc.sf <- sizeFactors(canc.dge)
    canc.sf
    print(summary(canc.sf))
    
    # Use combn() to get all possible combinations of factors for DGE testing - must do it pairwise
    factorpairs.mat <- combn(culture.factors, 2)

    # Use lapply() to transform the matrix into a list of column vectors
    ## Seq_len(ncol()) will produce a vector of 1:ncol()
    ## function(x) x[, i] will be iteratively applied to each element of the vector above
    ## and the return value will be placed in a list in the corresponding location
    contrasts <- lapply(seq_len(ncol(factorpairs.mat)), function(x) factorpairs.mat[, x])
    # dge.results <- lapply(1:length(contrasts), function(x) results(canc.dge, contrast = c("Culture", contrasts[[x]]), alpha = alphaTest))

    dge.results <- lapply(1:length(contrasts), function(x) results(canc.dge, contrast = c("Culture", contrasts[[x]]), alpha = alphaTest, lfcThreshold = lfc))

    # Now add HGNC Gene Symbols from Ensembl IDs
    for (i in 1:length(dge.results)) {
        dge.results[[i]]$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(dge.results[[i]]), keytype = "ENSEMBL", column = "SYMBOL")
    }

    return(dge.results)
}

# Keeps only significant DGE genes (p-adj < 0.05)
# Arguments:
    # results: a DESeqResults object
    # threshold (optional): the p-adjusted threshold for significance. Default is 0.05
# Return value: A DESeqResults object with only significant DE features
signifDE <- function(results, threshold = 0.05) {
    i <- is.na(results$padj)
    res.nona <- results[!i, ]
    return(res.nona[res.nona$padj < threshold, ])
}

#
volcano.plot <- function(de.results, lfc = 2, threshold = 0.05) {
    require(ggplot2)

    deseq.results <- as.data.frame(de.results)
    deseq.results$"Upregulated" <- deseq.results$log2FoldChange >= lfc & deseq.results$padj < threshold
    deseq.results$"Downregulated" <- deseq.results$log2FoldChange <= -lfc & deseq.results$padj < threshold
    deseq.results$"Threshold" <- as.factor(abs(deseq.results$log2FoldChange >= lfc) & deseq.results$padj < threshold)

    plt <- ggplot(data = deseq.results, aes(x = log2FoldChange, y = -log10(padj))) +
           geom_point(data = deseq.results, size = 1, color = "gray") +
           geom_point(data = deseq.results[deseq.results$"Upregulated" == TRUE, ], size = 2, color = "#CC0000") +
           geom_point(data = deseq.results[deseq.results$"Downregulated"  == TRUE, ], size = 2, colour = "#000099") +
           xlab("log2 fold change") +
           ylab("-log10 p-value adjusted") +
           ggtitle("My Title") +
           scale_x_continuous() +
           scale_y_continuous(limits = c(0, 15)) +
           theme_bw() +
           theme(axis.title.y = element_text(face = "bold", size = 16),
           axis.title.x = element_text(face = "bold", size = 16, colour = "black"),
           axis.text = element_text(size = 12),
           legend.title = element_blank() ,
           legend.text = element_text(size = 12)) +
           theme(plot.title = element_text(hjust = 0.5))

    return(plt)
}