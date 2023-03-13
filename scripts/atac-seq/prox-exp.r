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

# Code in names for ggplot
ggnames <- c("Pancreatic System", "MDA-MB Brain/Primary", "MDA-MB Lung/Primary")
colnames(atac.merged.up) <- ggnames

# Melt the merged df
atac.melted <- reshape2::melt(atac.merged.up)

# Make ggplot object
color.vec <- col.vec[c("panc", "brain_mb", "lung_mb")]
names(color.vec) <- ggnames
atac.plot.up <- ggplot(data = atac.melted, aes(x = variable, y = value, fill = factor(variable))) +
    geom_boxplot(color = "black") +
    geom_violin(alpha = 0.4) +
    theme_bw(base_size = 10) +
    scale_fill_manual(values = color.vec) +
    ylab("Log2 FC Gene Expression") +
    theme(axis.title.x = element_blank(),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          axis.text.x = element_text(size = 9),
          legend.position = "none") +
    annotate("text",
             x = seq_len(ncol(atac.merged.up)),
             y = sapply(atac.merged.up, function(x) max(na.omit(x)) + 0.5),
             size = 10,
             label = "***")

#----SAVE PLOT----
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/atac_all_binom.pdf")
atac.plot.up
dev.off()




#---PLOT FOR H3K27ac vs. ATAC-seq----
atac.prox.objs <- lapply(comp.names, function(x) get_proximal_genes(anno.atac.list.full[[x]], dge.list.full[[x]]))
chip.prox.objs <- lapply(comp.names, function(x) get_proximal_genes(anno.chip.list.full[[x]], dge.list.full[[x]]))

# Test on Ren data
ren.mutual <- intersect_prox(atac.prox.objs[[1]], chip.prox.objs[[1]], type = "up")
colnames(ren.mutual) <- c("Differential Chromatin Accessibility", "Differential H3K27ac", "Both")
ren.melted <- reshape2::melt(ren.mutual)

# ggplot
ggplot(data = ren.melted, aes(x = variable, y = value, fill = variable)) +
    geom_boxplot(color = "black") +
    geom_violin(alpha = 0.4) +
    theme_bw()
