require(BSDA)
require(cowplot)
require(DepecheR)

setClass("ProximalGeneExp", representation = representation(
    upregProx = "data.frame",
    downregProx = "data.frame",
    upregIntPercent = "numeric",
    downregIntPercent = "numeric",
    totalIntersect = "numeric",
    totalGenes = "numeric",
    binomTestUp = "numeric",
    binomTestDown = "numeric"
))

get_proximal_genes <- function(ranges, deseqObject, minlfc = 0, maxlfc = Inf) {
    # Get up/down changing features
    anno.up <- ranges[ranges$Fold > minlfc & ranges$Fold < maxlfc, ]
    anno.down <- ranges[ranges$Fold < -minlfc & ranges$Fold > -maxlfc, ]
    anno.up.features <- anno.up$feature
    anno.down.features <- anno.down$feature

    # Subset DESeq2 results object to isolate these features in particular
    anno.up.subset <- intersect(anno.up.features, rownames(deseqObject))
    anno.down.subset <- intersect(anno.down.features, rownames(deseqObject))
    gene.anno.up <- deseqObject[anno.up.subset, ] %>%
                        as.data.frame() %>%
                        select(log2FoldChange)
    gene.anno.down <- deseqObject[anno.down.subset, ] %>%
                        as.data.frame() %>%
                        select(log2FoldChange)

    # Calculate intersection stats
    upreg.int <- length(anno.up.subset) / length(anno.up.features)
    downreg.int <- length(anno.down.subset) / length(anno.down.features)
    total.int <- length(anno.up.subset) + length(anno.down.subset)

    # Calculate sign test for median > 0 or median < 0
    binomUp <- BSDA::SIGN.test(x = gene.anno.up$log2FoldChange, md = 0, alternative = "greater")$p.value
    binomDown <- BSDA::SIGN.test(x = gene.anno.down$log2FoldChange, md = 0, alternative = "less")$p.value

    # Instantiate ProximalGeneExp class and return
    pge.obj <- new("ProximalGeneExp",
                    upregProx = gene.anno.up,
                    downregProx = gene.anno.down,
                    upregIntPercent = upreg.int,
                    downregIntPercent = downreg.int,
                    totalIntersect = total.int,
                    totalGenes = length(rownames(deseqObject)),
                    binomTestUp = binomUp,
                    binomTestDown = binomDown)
    return(pge.obj)
}

make_bplot_pge <- function(pgeObject, plot = "up", fillCol = "white", 
                           notch = TRUE, name = "0", xpos = 1) {
    if (plot == "up") {
        log.df <- pgeObject@upregProx
        pval <- as.numeric(pgeObject@binomTestUp)
    } else if (plot == "down") {
        log.df <- pgeObject@downregProx
        pval <- as.numeric(pgeObject@binomTestDown)
    } else {
        stop("Parameter plot must be 'up' or 'down'.")
    }

    bplot <- geom_boxplot(data = log.df, 
                          aes(x = factor(name), y = log2FoldChange, fill = fillCol),
                          color = "black",
                          notch = TRUE)
    vplot <- geom_violin(data = log.df,
                         aes(x = factor(name), y = log2FoldChange, fill = fillCol),
                         alpha = 0.4)
    return(list(bplot, vplot))
}

#-----BASE PLOT-----
proximal.chip.list <- list()
for (i in 1:length(dge.list.full)) {
    prox.genes <- get_proximal_genes(anno.chip.list.full[[i]], dge.list.full[[i]])
    proximal.chip.list[[names(dge.list.full)[i]]] <- prox.genes
}

# Format p-values from exact binomial test for upregulated genes
pval.vec.up <- c()
for (i in 1:length(proximal.chip.list)) {
    pval.vec.up <- c(pval.vec.up, formatC(proximal.chip.list[[i]]@binomTestUp, digits = 3))
}

# Create boxplot
col.vec <- dColorVector(1:length(proximal.chip.list))
bplot.up.list <- lapply(1:length(proximal.chip.list), function(i) {
    make_bplot_pge(proximal.chip.list[[i]], plot = "up", fillCol = col.vec[i], name = names(proximal.chip.list[i]), xpos = i)
})

# Create and order annotations
names.sorted <- sort(names(dge.list.full), decreasing = FALSE)
label.list <- lapply(seq_len(length(names.sorted)), function(i) {
    j <- match(names.sorted[i], names(dge.list.full))
    annotate("text",
             x = i,
             y = max(proximal.chip.list[[j]]@upregProx$log2FoldChange) + 1,
             size = 10,
             label = "***")
})
bplot.up <- ggplot()
for (i in 1:length(bplot.up.list)) {
    bplot.up <- bplot.up + bplot.up.list[[i]][[1]] + bplot.up.list[[i]][[2]] + label.list[[i]]
}
bplot.up <- bplot.up +
             theme_bw(base_size = 11) +
             xlab("Comparison") +
             ylab("Log2FC") +
             theme(axis.text.x = element_text(size = 14),
                   axis.text.y = element_text(size = 14),
                   legend.position = "none") +
             geom_hline(yintercept = 0)

#------PLOT STRATIFIED BY H3K27AC LFC------
chip.zerotwo.list <- list()
for (i in 1:length(dge.list.full)) {
    prox.genes <- get_proximal_genes(anno.chip.list.full[[i]], dge.list.full[[i]], minlfc = 0, maxlfc = 2)
    chip.zerotwo.list[[names(dge.list.full)[i]]] <- prox.genes
}

chip.twofour.list <- list()
for (i in 1:length(dge.list.full)) {
    prox.genes <- get_proximal_genes(anno.chip.list.full[[i]], dge.list.full[[i]], minlfc = 2, maxlfc = 4)
    chip.twofour.list[[names(dge.list.full)[i]]] <- prox.genes
}

# Seems like there is a huge sample size difference between the 0 and 2 and 2 and 4 groups, particularly
# for the Rodrigues et al. data. We will see how this maps out

# Test for Ren et al
wilcox.test(x = chip.zerotwo.list[[1]]@upregProx$log2FoldChange,
            y = chip.twofour.list[[1]]@upregProx$log2FoldChange,
            paired = FALSE,
            alternative = "two.sided")
ggplot() +
    geom_boxplot(data = chip.zerotwo.list[[1]]@upregProx, aes(x = factor(0), y = log2FoldChange)) +
    geom_boxplot(data = chip.twofour.list[[1]]@upregProx, aes(x = factor(1), y = log2FoldChange))
