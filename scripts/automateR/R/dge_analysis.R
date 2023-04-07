#' @title Reading in raw RNA-seq data
#'
#' @description
#' Converts a featureCounts count file into a data frame for downstream analysis.
#'
#' @param featurecounts.txt Path to valid `featureCounts` filename
#'
#' @return A data frame containing columns corresponding to samples and rows corresponding to raw gene counts.
#' @export
#'
#' @examples
#' # Reading in a featureCounts file from directory "~/counts"
#' # featurecounts.txt <- file.path("~/counts", "counts.txt")
#' # txt2counts(featurecounts.txt)
txt2counts <- function(featurecounts.txt) {
    # Read file, skip first line which contains terminal command
    featurecounts <- utils::read.delim(featurecounts.txt, header = TRUE, skip = 1)

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

#' @title Ensembl ID conversion
#'
#' @description
#' Takes a vector of Ensembl IDs and converts them to the corresponding HGNC gene symbols using the `EnsDb.Hsapiens.v79` package.
#'
#' @param ensembl.genes A character vector of Ensembl IDs.
#'
#' @return A character vector containing the matching HGNC symbol for each Ensembl ID.
#' NA indicates the symbol could not be matched to a gene feature.
#' @export
#'
#' @examples
#' ens.vec <- c("ENSG00000166035", "ENSG00000102554", "ENSG00000141510")
#' ensembl_to_gene(ens.vec)
ensembl_to_gene <- function(ensembl.genes) {
    # See https://stackoverflow.com/questions/28543517/ for converting ENSEMBL to HGNC
    geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79,
                                 keys = ensembl.genes,
                                 keytype = "GENEID",
                                 columns = c("SYMBOL", "GENEID"))
    # Ensembl IDs mapping to NA are excluded
    # Need to manually match each original Ensembl ID to corresponding HGNC symbol
    # If no Ensembl ID is included in geneIDs, assign it NA for HGNC symbol
    hgnc.sym <- c()
    for (i in seq_len(length(ensembl.genes))) {
        id.index <- match(ensembl.genes[i], geneIDs[, 2])
        if (is.na(id.index)) {
            hgnc.sym <- c(hgnc.sym, NA)
        } else {
            hgnc.sym <- c(hgnc.sym, geneIDs[id.index, 1])
        }
    }
    return(hgnc.sym)
}


#' @title Automate DGE from `DESeq2`
#'
#' @description
#' Takes in raw count data and a design matrix and automates repeated calls to `DESeq2`` for user convenience.
#' @param featurecounts Can either be a character string representing path to valid `featureCounts` matrix file,
#' or an isolated data frame of raw gene counts.
#' @param coldata A df-like design matrix object for `DESeq2` - this assumes the user has already factored and leveled the columns.
#' @param contrast.var Character string indicating which variable to perform DGE for.
#' Default: Assumes last variable in `coldata`.
#' @param contrasts Character string for levels to be pairwise tested. Levels listed first will be in the numerator.
#' Default: will perform pairwise DGE tests on all levels for the last column in the design matrix.
#' @param formula.vec Character vector indicating the columns or interactions to use for the design for `DESeq2`.
#' Default: Will use all columns available in the design matrix, with the last column being the main effect of interest.
#' @param alphaTest A FDR-adjusted p-value to determine significance for DE. Default is 0.05.
#' @param lfc A threshold to pass to `lfcShrink()` from `DESeq2` to determine the fold-change threshold for significance.
#' Default: 0.
#'
#' @return A list of 1) list of `DESeqResults` objects, 2) the original `DGE` object from `DESeq()`
#' @export
#'
#' @examples
#' # For analyzing data from GSE149103
#' featurecounts.ren <- ("D:/SK/data/ren-panc/rna-seq/aligned/ren_rna_rawcounts.txt")
#' coldata.ren <- data.frame(Culture = rep(c("CAPAN1", "HPNE", "PANC1"), each = 3))
#' coldata.ren$Culture <- factor(coldata.ren$Culture, levels = c("CAPAN1", "PANC1", "HPNE"))
#'
#' # Call DESeq2, isolate capan vs. panc comparison
#' ren.dge <- dge_analysis(featurecounts.ren, coldata.ren, contrasts = c("CAPAN1", "PANC1", "HPNE"))
#' capan.panc <- ren.dge[[1]][[1]]
dge_analysis <- function(featurecounts,
                         coldata,
                         contrast.var = colnames(coldata)[length(colnames(coldata))],
                         contrasts = coldata[, colnames(coldata) == contrast.var] %>% levels(),
                         formula.vec = colnames(coldata),
                         alphaTest = 0.05,
                         lfc = 0) {
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
    keepCount <- rowSums(counts(canc.dds) >= 10) >= ncol(counts.df) / 3
    canc.dds <- canc.dds[keepCount, ]
    canc.dge <- DESeq(canc.dds)

    # Save size factors
    canc.sf <- sizeFactors(canc.dge)
    print(canc.sf)
    print(summary(canc.sf))

    # Use combn() to get all possible combinations of factors for DGE testing - must do it pairwise
    factorpairs.mat <- combn(contrasts, 2)

    # Use lapply() to transform the matrix into a list of column vectors
    ## Seq_len(ncol()) will produce a vector of 1:ncol()
    ## function(x) x[, i] will be iteratively applied to each element of the vector above
    ## and the return value will be placed in a list in the corresponding location
    contrasts <- lapply(seq_len(ncol(factorpairs.mat)), function(x) factorpairs.mat[, x])
    dge.results <- lapply(1:length(contrasts),
                          function(x) lfcShrink(canc.dge,       # need to include call to results() to set alphaTest manually
                                                res = results(  # for lfcShrink, alpha = 0.10
                                                              canc.dge,
                                                              contrast = c(contrast.var, contrasts[[x]]),
                                                              alpha = alphaTest,
                                                              lfcThreshold = lfc
                                                ),
                                                type = "normal",
                                                contrast = c(contrast.var, contrasts[[x]])))

    for (i in 1:length(dge.results)) {
        ensembl.genes <- rownames(dge.results[[i]])
        dge.results[[i]]$Symbol <- ensembl_to_gene(ensembl.genes)
    }

    return(list(dge.results, canc.dge))
}

#' @title Isolate DGE genes
#'
#' @description
#' Takes in a `results` object produced by `DESeq2` and isolates only differentially expressed genes.
#' @param results Object produced by a call to `results` or `lfcShrink` in `DESeq2`
#' @param threshold FDR-adjusted p-value threshold. Only genes below this threshold will be kept. Default: 0.05.
#' @param lfc Log-fold change minimum to be imposed to extract DE genes. Default = 0 (i.e., no threshold)
#'
#' @return Subset of the `results` object, containing only genes below the user-defined threshold.
#' @export
signifDE <- function(results, threshold = 0.05, lfc = 0) {
    i <- is.na(results$padj)
    res.nona <- results[!i, ]
    signif <- res.nona[res.nona$padj < threshold, ]
    return(signif[abs(signif$log2FoldChange) > lfc, ])
}

#' @title Isolate up- and down-regulated genes
#'
#' @description
#' Splits a `DESeq2 results` object into genes that are upregulated in the contrast (LFC > 0) and
#' genes that are downregulated in the contrast (LFC < 0). Makes an internal call to `signifDE()`
#' such that only genes whose change in regulation is statistically significant is returned.
#' @param ... Arguments to be passed to `signifDE` (see `?signifDE` for more)
#'
#' @return A list of 2 objects, containing 1) `results` subset that has significantly upregulated genes,
#' and 2) a `results` subset that has significantly downregulated genes.
#' @export
splitDE <- function(...) {
    signif.res <- signifDE(...)
    upreg.res <- signif.res[signif.res$log2FoldChange > 0, ]
    downreg.res <- signif.res[signif.res$log2FoldChange < 0, ]
    return(list(upreg.res, downreg.res))
}

#' @title Volcano plot
#'
#' @description
#' Produces a volcano plot displaying -log10(p-adj) vs. LFC for genes given a `DESeq2 results` object.
#' @param de.results Object produced by a call to DESeq2 `results`.
#' @param lfc Threshold for determining significant change in gene expression.
#' E.g., lfc = 1 indicates only genes with a p-adjusted < `threshold` and `abs(LFC)` > `lfc` should be deemed
#' significantly different in gene expression across the comparison samples. Default = 0 (i.e., no LFC threshold)
#' @param threshold FDR-adjusted p-value cutoff for statistical significance. Default = 0.05.
#'
#' @return Plot object produced by a call to `ggplot` from the `ggplot2` package.
#' @export
volcano.plot <- function(de.results, lfc = 0, threshold = 0.05) {
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
