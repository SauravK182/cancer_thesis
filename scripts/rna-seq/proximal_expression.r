require(reshape2)

#-----BASE PLOT-----
# Get LFCs of genes annotated to be near a significant increase in H3K27ac occupancy
proximal.up.list <- lapply(seq_len(length(dge.list.full)), function(i) {
    prox.genes <- get_proximal_genes(anno.chip.list.full[[i]], dge.list.full[[i]])@upregProx
})
names(proximal.up.list) <- names(dge.list.full)

# Get LFCs of genes annotated to be near a significant decrease in H3K27ac occupancy
proximal.down.list <- lapply(seq_len(length(dge.list.full)), function(i) {
    prox.genes <- get_proximal_genes(anno.chip.list.full[[i]], dge.list.full[[i]])@downregProx
})
names(proximal.down.list) <- names(dge.list.full)

# Within each list, merge all the dataframes in the list
merged.df.list <- lapply(list(up = proximal.up.list, down = proximal.down.list), function(ls) {
    merged.df <- purrr::reduce(.x = ls, .f = function(df1, df2) {
        merge(df1, df2, by = 0, all = TRUE) %>%
        remove_rownames() %>%
        column_to_rownames("Row.names")
    })
    colnames(merged.df) <- names.comp
    return(merged.df)
})
melted.df.list <- lapply(merged.df.list, reshape2::melt)

# Creat plots
bplot.list <- lapply(seq_len(length(melted.df.list)), function(i) {
    bplot <- ggplot(data = melted.df.list[[i]], aes(x = variable, y = value, fill = variable)) +
                geom_boxplot(color = "black") +
                stat_boxplot(geom = "errorbar", width = 0.3) +
                geom_violin(alpha = 0.4) +
                theme_cowplot() +
                ylab("Log2 FC Gene Expression") +
                geom_hline(yintercept = 0) +
                theme(axis.title.x = element_blank(),
                      legend.position = "none",
                      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
                annotate("text",
                         x = seq_len(ncol(merged.df.list[[i]])),
                         y = sapply(merged.df.list[[i]], function(x) max(x[!is.na(x)]) + 1),
                         size = 7,
                         label = "***") +
                scale_fill_manual(values = unname(col.vec))
})

#----SAVE THE BASE PLOTS----
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/k27_all_binom_up.pdf")
bplot.list[[1]]
dev.off()

cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/k27_all_binom_down.pdf")
bplot.list[[2]]
dev.off()

cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/k27_both_directions.pdf", height = 10, width = 6)
# Remove x-labels of first plot since they're redundant
bplot.up.noname <- bplot.list[[1]] + theme(axis.text.x = element_blank())
cowplot::plot_grid(plotlist = list(bplot.up.noname, bplot.list[[2]]), ncol = 1, nrow = 2, labels = "AUTO",
                   rel_heights = c(0.70, 1))
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
col.vec2 <- c("blue", "red")
comp.plot.list <- lapply(1:length(merged.list), function(i) {
    ggplot(data = merged.list[[i]], aes(x = variable, y = value, fill = variable)) +
        geom_boxplot(color = "black", width = 0.5) +
        stat_boxplot(geom = "errorbar", width = 0.25) +
        geom_violin(alpha = 0.4, width = 0.65) +
        geom_signif(comparisons = list(c("H3K27ac LFC 0-2", "H3K27ac LFC 2-4")),
                    map_signif_level = TRUE,
                    test = "wilcox.test",
                    test.args = list(alternative = "two.sided", paired = FALSE),
                    size = 0.4,
                    textsize = 4,
                    vjust = 0.5) +
        theme_cowplot() +
        ylab("Log2 FC Gene Expression") +
        theme(axis.title.x = element_blank(),
              axis.text.y = element_text(size = 9),
              axis.title.y = element_text(size = 9),
              axis.text.x = element_text(size = 9),
              plot.title = element_text(size = 9),
              legend.position = "none") +
        ggtitle(names.comp[i]) +
        scale_fill_manual(values = col.vec2)
})

# Save plot
# Use nested plot_grids to center bottom sole plot - modify the widths of the 2nd plot to
# center the lone graph by using NULL space
bp.k27.comp <- cowplot::plot_grid(
                        cowplot::plot_grid(comp.plot.list[[1]], 
                                           comp.plot.list[[2]], 
                                           comp.plot.list[[3]], 
                                           comp.plot.list[[4]], ncol = 2, nrow = 2, labels = "AUTO"),
                        cowplot::plot_grid(NULL, comp.plot.list[[5]], NULL, rel_widths = c(0.5, 1, 0.5), nrow = 1, labels = "E", hjust = -14),
                        nrow = 2,
                        rel_heights = c(2, 1))
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/k27_plot_compare.pdf")
bp.k27.comp
dev.off()

#------PLOT STRATIFIED BY H3K27AC LFC W/ MORE REFINED THRESHOLDS-------
h3.comps <- c("H3K27ac LFC 0-1", "H3K27ac LFC 1-2", "H3K27ac LFC > 2")
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
    colnames(merged.df) <- h3.comps
    iscol.na <- c()
    for (j in seq_len(ncol(merged.df))) {
        if (sum(is.na(merged.df[, j])) == length(merged.df[, j])) {
            iscol.na <- c(iscol.na, TRUE)
        } else {
            iscol.na <- c(iscol.na, FALSE)
        }
    }
    merged.list.refined[[i]] <- merged.df[, !iscol.na]
}

# Create color vector and melt
melted.refined.list <- lapply(merged.list.refined, reshape2::melt)
col.vec3 <- c("blue", "green", "red")
# Make plots
refined.plot.list <- lapply(1:length(melted.refined.list), function(i) {
    comp.list <- split(t(combn(levels(melted.refined.list[[i]]$variable), 2)), 
                       seq(nrow(t(combn(levels(melted.refined.list[[i]]$variable), 2)))))
    colors <- col.vec3[match(levels(melted.refined.list[[i]]$variable), h3.comps)]
    ggplot(data = melted.refined.list[[i]], aes(x = variable, y = value, fill = variable)) +
        geom_boxplot(color = "black", aes(group = variable)) +
        stat_boxplot(geom = "errorbar", width = 0.25) +
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
        theme_cowplot() +
        ylab("Log2 FC Gene Expression") +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_text(size = 9),
              axis.text.x = element_text(size = 7),
              plot.title = element_text(size = 9),
              legend.position = "none") +
        ggtitle(names.comp[i]) +
        scale_fill_manual(values = colors)
})

#---SAVE PLOTS---
# Use multiple plot_grids to center lone 5th graph
k27.comp.trio.plot <- cowplot::plot_grid(
                        cowplot::plot_grid(refined.plot.list[[1]], 
                                           refined.plot.list[[2]], 
                                           refined.plot.list[[3]], 
                                           refined.plot.list[[4]], ncol = 2, nrow = 2, labels = "AUTO"),
                        cowplot::plot_grid(NULL, refined.plot.list[[5]], NULL, rel_widths = c(0.5, 1, 0.5), nrow = 1, labels = "E", hjust = -14),
                        nrow = 2,
                        rel_heights = c(2, 1)
)
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
