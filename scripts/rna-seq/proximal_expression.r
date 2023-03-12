require(BSDA)
require(cowplot)
require(DepecheR)
require(reshape2)
require(ggsignif)
require(gridExtra)

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
                          aes(x = factor(name), y = log2FoldChange),
                          fill = fillCol,
                          color = "black",
                          notch = TRUE)
    vplot <- geom_violin(data = log.df,
                         aes(x = factor(name), y = log2FoldChange),
                         fill = fillCol,
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

#----SAVE THE BASE PLOT----
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/k27_all_binom.pdf")
bplot.up
dev.off()

#------PLOT STRATIFIED BY H3K27AC LFC------
chip.zerotwo.list <- list()
chip.twofour.list <- list()
merged.list <- list()
for (i in 1:length(dge.list.full)) {
    prox.genes.zero <- get_proximal_genes(anno.chip.list.full[[i]], dge.list.full[[i]], minlfc = 0, maxlfc = 2)
    prox.genes.two <- get_proximal_genes(anno.chip.list.full[[i]], dge.list.full[[i]], minlfc = 2, maxlfc = 4)
    chip.zerotwo.list[[names(dge.list.full)[i]]] <- prox.genes.zero
    chip.twofour.list[[names(dge.list.full)[i]]] <- prox.genes.two
    merged.df <- merge(prox.genes.zero@upregProx, prox.genes.two@upregProx, by = 0, all = TRUE) %>%
                    remove_rownames() %>%
                    column_to_rownames(var = "Row.names")
    colnames(merged.df) <- c("H3K27ac LFC 0-2", "H3K27ac LFC 2-4")
    merged.list[[names(dge.list.full)[i]]] <- merged.df
}

# Reshape merged data frames to long format to be used for ggplot
merged.list <- lapply(merged.list, reshape2::melt)
names.comp <- c("Pancreatic System",
                "786 ccRCC System",
                "OS ccRCC System",
                "BrM2 Brain vs. Primary",
                "LM2 Lung vs. Primary")
comp.plot.list <- lapply(1:length(merged.list), function(i) {
    ggplot(data = merged.list[[i]], aes(x = variable, y = value, fill = variable)) +
        geom_boxplot(color = "black") +
        geom_violin(alpha = 0.4) +
        geom_signif(comparisons = list(c("H3K27ac LFC 0-2", "H3K27ac LFC 2-4")),
                    map_signif_level = TRUE,
                    test = "wilcox.test",
                    test.args = list(alternative = "two.sided", paired = FALSE),
                    size = 0.4,
                    textsize = 4,
                    vjust = 0.5) +
        theme_bw(base_size = 9) +
        ylab("Log2 FC") +
        theme(axis.title.x = element_blank(),
              axis.text.y = element_text(size = 9),
              plot.title = element_text(size = 9),
              legend.position = "none") +
        ggtitle(names.comp[i])
})

# Save plot
bp.k27.comp <- plot_grid(plotlist = comp.plot.list, ncol = 2, labels = "AUTO")
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/k27_plot_compare.pdf")
bp.k27.comp
dev.off()

#------PLOT STRATIFIED BY H3K27AC LFC W/ MORE REFINED THRESHOLDS-------
merged.list.refined <- list()
for (i in 1:length(dge.list.full)) {
    prox.genes.zero <- get_proximal_genes(anno.chip.list.full[[i]], dge.list.full[[i]], minlfc = 0, maxlfc = 1)
    prox.genes.one <- get_proximal_genes(anno.chip.list.full[[i]], dge.list.full[[i]], minlfc = 1, maxlfc = 2)
    prox.genes.two <- get_proximal_genes(anno.chip.list.full[[i]], dge.list.full[[i]], minlfc = 2, maxlfc = Inf)
    prox.genes <- list(prox.genes.zero@upregProx, prox.genes.one@upregProx, prox.genes.two@upregProx)
    merged.df <- purrr::reduce(prox.genes, .f = function(df1, df2) {
        merge(df1, df2, by = 0, all = TRUE) %>%
        remove_rownames() %>%
        column_to_rownames("Row.names")
    })
    colnames(merged.df) <- c("H3K27ac LFC 0-1", "H3K27ac LFC 1-2", "H3K27ac LFC > 2")
    merged.list.refined[[i]] <- merged.df
}

melted.refined.list <- lapply(merged.list.refined, reshape2::melt)
refined.plot.list <- lapply(1:length(melted.refined.list), function(i) {
    comp.list <- split(t(combn(levels(melted.refined.list[[i]]$variable), 2)), 
                       seq(nrow(t(combn(levels(melted.refined.list[[i]]$variable), 2)))))
    ggplot(data = melted.refined.list[[i]], aes(x = variable, y = value, fill = variable)) +
        geom_boxplot(color = "black", aes(group = variable)) +
        geom_violin(alpha = 0.4) +
        geom_signif(comparisons = comp.list,
                    map_signif_level = c("***" = 0.001 / length(comp.list),
                                         "**" = 0.01 / length(comp.list),
                                         "*" = 0.05 / length(comp.list)),
                    test = wilcox.test,
                    test.args = list(alternative = "less", paired = FALSE),
                    step_increase = 0.1,
                    size = 0.4,
                    textsize = 2.5,
                    vjust = 0.2) +
        theme_bw(base_size = 9) +
        ylab("Log2 FC Gene Expression") +
        theme(axis.title.x = element_blank(),
              axis.text.y = element_text(size = 9),
              plot.title = element_text(size = 9),
              legend.position = "none") +
        ggtitle(names.comp[i])
})

#---SAVE PLOTS---
k27.comp.trio.plot <- plot_grid(plotlist = refined.plot.list, ncol = 2, labels = "AUTO")
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/k27_triple_plot_compare.pdf")
k27.comp.trio.plot
dev.off()

# Seems like there is a huge sample size difference between the 0 and 2 and 2 and 4 groups, particularly
# for the Rodrigues et al. data. We will see how this maps out

# Test for Ren et al
wilcox.test(x = chip.zerotwo.list[[1]]@upregProx$log2FoldChange,
            y = chip.twofour.list[[1]]@upregProx$log2FoldChange,
            paired = FALSE,
            alternative = "less")
