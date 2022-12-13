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

ensembl_to_gene <- function(ensembl.genes) {
    geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79,
                                 keys = ensembl.genes,
                                 keytype = "GENEID",
                                 columns = c("SYMBOL", "GENEID"))
    hgnc.sym <- c()
    for (i in seq_len(nrow(geneIDs))) {
        geneID.ind <- match(geneIDs[i, 2], ensembl.genes)
        hgnc.sym[geneID.ind] <- geneIDs[i, 1]
    }
    return(hgnc.sym)
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
dge_analysis <- function(featurecounts, coldata, formula.vec = colnames(coldata), contrasts, alphaTest = 0.05, lfc = 2) {
    require(DESeq2)
    require(EnsDb.Hsapiens.v79)

    if (class(featurecounts) == "character") {
        counts.df <- txt2counts(featurecounts)      # assume file input
    } else if ("data.frame" %in% class(featurecounts)) {
        counts.df <- as.data.frame(featurecounts)   # standardize input
    } else {
        stop("Argument 'featurecounts' must be a data frame-like object or valid input file.")
    }

    if ("data.frame" %in% class(coldata)) {
        design.mat <- as.data.frame(coldata)
    } else {
        stop("Argument 'coldata' must be a data frame-like object for passing to DESeq2.")
    }

    # Generate formula for design from formula.vec
    collapsed.vec <- paste(formula.vec, collapse = "+") # will collapse vector "a b c" into "a+b+c"
    des.formula <- as.formula(paste0("~", collapsed.vec))   # paste as formula object for use in DDS

    # Run DESeq to get the DESeqDataSet, keep only features with > 10 counts
    canc.dds <- DESeqDataSetFromMatrix(countData = counts.df,
                                    colData = design.mat,
                                    design = des.formula)
    keepCount <- rowSums(counts(canc.dds)) >= 10
    canc.dds <- canc.dds[keepCount, ]
    canc.dge <- DESeq(canc.dds)

    # Save size factors
    canc.sf <- sizeFactors(canc.dge)
    print(summary(canc.sf))

    # Use combn() to get all possible combinations of factors for DGE testing - must do it pairwise
    factorpairs.mat <- combn(contrasts, 2)

    # Use lapply() to transform the matrix into a list of column vectors
    ## Seq_len(ncol()) will produce a vector of 1:ncol()
    ## function(x) x[, i] will be iteratively applied to each element of the vector above
    ## and the return value will be placed in a list in the corresponding location
    contrasts <- lapply(seq_len(ncol(factorpairs.mat)), function(x) factorpairs.mat[, x])
    # dge.results <- lapply(1:length(contrasts), function(x) results(canc.dge, contrast = c("Culture", contrasts[[x]]), alpha = alphaTest))

    dge.results <- lapply(1:length(contrasts),
                          function(x) results(canc.dge,
                                              contrast = c("Culture", contrasts[[x]]),
                                              alpha = alphaTest,
                                              lfcThreshold = lfc))

    # Now add HGNC Gene Symbols from Ensembl IDs
    for (i in 1:length(dge.results)) {
        ensembl.genes <- rownames(dge.results[[i]])
        geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79,
                              keys = ensembl.genes,
                              keytype = "GENEID",
                              columns = c("SYMBOL", "GENEID"))
        hgnc.sym <- c()
        for (i in seq_len(nrow(geneIDs))) {
            geneID.ind <- match(geneIDs[i, 2], ensembl.genes)
            hgnc.sym[geneID.ind] <- geneIDs[i, 1]
        }
        dge.results[[i]]$Symbol <- hgnc.sym
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