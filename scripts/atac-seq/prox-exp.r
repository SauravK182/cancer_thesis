require(cowplot)

color_df <- function(melted.df, color.vec) {
    if (! ("variable" %in% colnames(melted.df))) {
        stop("Must pass a melted data frame with column 'variable' indicating original colnames.")
    }

    groups <- unique(melted.df$variable)
    print(groups)
    if (length(groups) != length(color.vec)) {
        stop("Must pass a color vector of equal size as the number of groups in the 'variable' column.")
    }

    melted.df$color <- color.vec[match(melted.df$variable, groups)]
    return(melted.df)
}

comp.names <- names(anno.atac.list.full)

#----BASE PLOT----
atac.prox.list <- lapply(comp.names, function(i) get_proximal_genes(anno.atac.list.full[[i]], dge.list.full[[i]])@upregProx)
atac.merged.up <- purrr::reduce(atac.prox.list, .f = function(df1, df2) {
    merge(df1, df2, by = 0, all = TRUE) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Row.names")
})

# Code in names for ggplot and get color vector
ggnames <- c("Pancreatic System", "BrM2 Brain vs. Primary", "LM2 Lung vs. Primary")
colnames(atac.merged.up) <- ggnames
atac.color <- col.vec[match(comp.names, names(col.vec))]


# Melt the merged df
atac.melted <- reshape2::melt(atac.merged.up)

# Make ggplot object
names(color.vec) <- ggnames
atac.plot.up <- ggplot(data = atac.melted, aes(x = variable, y = value, fill = factor(variable))) +
    geom_boxplot(color = "black") +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_violin(alpha = 0.4) +
    theme_cowplot() +
    scale_fill_manual(values = unname(atac.color)) +
    ylab("Log2 FC Gene Expression") +
    theme(axis.title.x = element_blank(),
          legend.position = "none") +
    annotate("text",
             x = seq_len(ncol(atac.merged.up)),
             y = sapply(atac.merged.up, function(x) max(na.omit(x)) + 0.5),
             size = 8,
             label = "***") +
    geom_hline(yintercept = 0)

#----SAVE PLOT----
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/atac_all_binom.pdf")
atac.plot.up
dev.off()




#---PLOT FOR H3K27ac vs. ATAC-seq----
atac.prox.objs <- lapply(comp.names, function(x) get_proximal_genes(anno.atac.list.full[[x]], dge.list.full[[x]]))
chip.prox.objs <- lapply(comp.names, function(x) get_proximal_genes(anno.chip.list.full[[x]], dge.list.full[[x]]))

# Get list of melted dfs for plots
atac.chip.df <- lapply(seq_len(length(atac.prox.objs)), function(i) {
    df.prox <- intersect_prox(atac.prox.objs[[i]], chip.prox.objs[[i]])
    colnames(df.prox) <- c("Differential Chromatin Accessibility", "Differential H3K27ac", "Both")
    return(df.prox)
})
atac.chip.df.melted <- lapply(atac.chip.df, reshape2::melt)

# Make plot with unpaired Wilcoxon rank-sum test for population distributions
xlabels <- expression(
    Delta~Accesibility,
    Delta~H3K27ac,
    Both
)
col.vec3 <- c("blue", "green", "red")
ggplot.list <- lapply(seq_len(length(atac.chip.df.melted)), function(i) {
    comp.list <- split(t(combn(levels(atac.chip.df.melted[[i]]$variable), 2)),
                       seq(nrow(t(combn(levels(atac.chip.df.melted[[i]]$variable), 2)))))
    bplot <- ggplot(data = atac.chip.df.melted[[i]], aes(x = variable, y = value, fill = variable)) +
                geom_boxplot(color = "black", aes(group = variable)) +
                stat_boxplot(geom = "errorbar", width = 0.2) +
                geom_violin(alpha = 0.4) +
                theme_cowplot() +
                geom_signif(
                    comparisons = comp.list,
                    map_signif_level = c("*" = 0.05 / length(comp.list),
                                         "**" = 0.01 / length(comp.list),
                                         "***" = 0.001 / length(comp.list)),
                    test = wilcox.test,
                    test.args = list(paired = FALSE, alternative = "less"),
                    step_increase = 0.1,
                    size = 0.4,
                    textsize = 2.5,
                    vjust = 0.2
                ) +
                annotate("text",
                         x = seq_len(length(levels(atac.chip.df.melted[[i]]$variable))) + 0.33,
                         y = sapply(atac.chip.df[[i]], function(x) min(na.omit(x)) - 0.5),
                         label = paste("n =", sapply(atac.chip.df[[i]], function(x) sum(!is.na(x)))),
                         size = 3) +
                ggtitle(ggnames[i]) +
                ylab("Log2 FC Gene Expression") +
                theme(axis.title.x = element_blank(),
                      axis.title.y = element_text(size = 10),
                      axis.text.x = element_text(size = 9),
                      plot.title = element_text(size = 9),
                      legend.position = "none") +
                scale_x_discrete(labels = xlabels) +
                scale_fill_manual(values = col.vec3)
    return(bplot)
})

# Use plot_grid and save the plots
arranged.plots <- cowplot::plot_grid(
    cowplot::plot_grid(ggplot.list[[1]], ggplot.list[[2]], nrow = 1, ncol = 2, labels = "AUTO"),
    cowplot::plot_grid(NULL, ggplot.list[[3]], NULL, rel_widths = c(0.5, 1, 0.5), nrow = 1, labels = "C", hjust = -13),
    nrow = 2
)
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/atac_chip_comp.pdf")
arranged.plots
dev.off()

# Note that p-value for Brain vs. Primary of both vs h3k27ac is 0.04359
wilcox.test(x = atac.chip.df[[2]][, 2], y = atac.chip.df[[2]][, 3], paired = FALSE)

# p-value for Lung vs. Primary of both vs h3k27ac is 0.0192
wilcox.test(x = atac.chip.df[[3]][, 2], y = atac.chip.df[[3]][, 3], paired = FALSE, alternative = "less")
