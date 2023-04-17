keep_levels <- function(x) {
    return(factor(x, levels = unique(x)))
}

names.comp <- c("Pancreatic", "786 ccRCC", "OS ccRCC", "BrM2", "LM2")

#-----------PRELIMINARY-----------
# Read in pluripotency .bed file
pluri.table <- read.table("pluri_wnt_genes.bed",
                          col.names = c("Chromosome", "Start", "End", "Gene", "Strand"))

# Get Ensembl IDs with ensembldb and EnsDb package
pluri.ensembl <- ensembldb::select(EnsDb.Hsapiens.v79,
                                   keys = pluri.table$Gene,
                                   keytype = "SYMBOL",
                                   column = "GENEID") %>%
                                   dplyr::select(GENEID) %>% deframe()

# Get only IDs recognized by featureCounts
pluri.ensfinal <- intersect(pluri.ensembl, rownames(full.df))
num.genes.mapped <- ensembldb::select(EnsDb.Hsapiens.v79,
                                      keys = pluri.ensfinal,
                                      keytype = "GENEID",
                                      column = "SYMBOL") %>%
                                      select(SYMBOL) %>%
                                      deframe() %in% pluri.table$Gene %>% sum()


#----------DGE ANALYSIS----------
# Get genes that are changing DE vs. not
de.pluri <- lapply(dge.list.full, function(df) {
    tested <- intersect(pluri.ensfinal, rownames(df))
    upreg <- intersect(pluri.ensfinal, rownames(splitDE(df)[[1]]))
    downreg <- intersect(pluri.ensfinal, rownames(splitDE(df)[[2]]))
    not.de <- length(tested) - length(upreg) - length(downreg)
    return(c("Not Tested" = num.genes.mapped - length(tested), "No Change" = not.de,
             "Downregulated" = length(downreg), "Upregulated" = length(upreg)))
})

# Transform into proportions out of num.genes.mapped
de.pluri.prop <- lapply(de.pluri, function(x) x / num.genes.mapped)

# Create df for ggplot
names(de.pluri.prop) <- names.comp
de.pluri.prop %<>% as.data.frame(check.names = FALSE) %>%
                    rownames_to_column(var = "Expr") %>% 
                    reshape2::melt()


# Make barplot
pluri.dge.bar <- ggplot(data = de.pluri.prop, aes(x = variable, y = value, fill = keep_levels(Expr))) +
                    geom_bar(stat = "identity", color = "black", width = 0.75) +
                    theme_cowplot() +
                    ylab("Fraction of geneset") +
                    theme(axis.title.x = element_blank(),
                          axis.text.x = element_text(size = 9),
                          axis.title.y = element_text(size = 9)) +
                    scale_fill_manual(values = c("black", "gray", "darkblue", "red"),
                                      name = expression(Delta~Expression)) +
                    scale_y_continuous(expand = c(0,0))




#----------CHIP/ATAC ANALYSIS---------
epi.list <- list(chipseq = anno.chip.list.full, atacseq = anno.atac.list.full)
epi.list %<>% lapply(function(list) lapply(list, function(gr) {
    df <- gr %>% as.data.frame(row.names = NULL) %>% dplyr::select(feature, Fold)
    upreg <- intersect(pluri.ensfinal, filter(df, Fold > 0)$feature)
    downreg <- intersect(pluri.ensfinal, filter(df, Fold < 0)$feature)
    return(c("No Change" = num.genes.mapped - length(upreg) - length(downreg),
             "Depletion" = length(downreg), "Enrichment" = length(upreg)))
}))
names(epi.list[[1]]) <- names.comp
names(epi.list[[2]]) <- names.comp[c(1, 4, 5)]
# Create dfs & plots from epigenomic data
epi.dfs <- lapply(epi.list, function(list) {
    prop.list <- lapply(list, function(x) x / num.genes.mapped)
    return(
        as.data.frame(prop.list, check.names = FALSE) %>%
            rownames_to_column(var = "Change") %>%
            reshape2::melt()
    )
})
names <- c(expression(Delta~H3K27ac, Delta~Accessibility))
epi.plots <- lapply(seq_along(epi.dfs), function(x) {
    ggplot(data = epi.dfs[[x]], aes(x = variable, y = value, fill = keep_levels(Change))) +
        geom_bar(stat = "identity", color = "black", width = 0.75) +
        theme_cowplot() +
        ylab("Fraction of geneset") +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 9),
              axis.title.y = element_text(size = 9)) +
        scale_fill_manual(values = c("gray", "darkblue", "red"),
                          name = names[x]) +
        scale_y_continuous(expand = c(0, 0))
})




#--------HI-C ANALYSIS----------
hic.diffcomp <- lapply(hic.anno.list, function(gr) {
    df <- gr %>%
            as.data.frame(row.names = NULL) %>%
            dplyr::select(feature, transition)
    total.signif <- length(intersect(pluri.ensfinal, df$feature))
    b.to.a <- intersect(pluri.ensfinal, unique(filter(df, transition == "B to A Transition")$feature))
    a.to.b <- intersect(pluri.ensfinal, unique(filter(df, transition == "A to B Transition")$feature))
    within.comp <- total.signif - length(b.to.a) - length(a.to.b)
    no.change <- num.genes.mapped - total.signif
    return(c("No Change" = no.change, "Within Compartment Change" = within.comp,
             "A to B Transition" = length(a.to.b), "B to A Transition" = length(b.to.a)) / num.genes.mapped)
})
names(hic.diffcomp) <- c("Pancreatic", "786 ccRCC")
hic.df <- as.data.frame(hic.diffcomp, check.names = FALSE) %>%
            rownames_to_column(var = "transition") %>%
            reshape2::melt()
hic.barplot <- ggplot(data = hic.df, aes(x = variable, y = value, fill = keep_levels(transition))) +
                geom_bar(stat = "identity", color = "black", width = 0.75) +
                theme_cowplot() +
                ylab("Fraction of geneset") +
                theme(axis.title.x = element_blank(),
                      axis.title.y = element_text(size = 9)) +
                scale_fill_manual(values = c("gray", "orange", "blue", "red"),
                                  name = expression(Delta~Compartment)) +
                scale_y_continuous(expand = c(0, 0))

#----------SAVE PLOT----------
plotlist <- list(pluri.dge.bar, epi.plots[[1]], epi.plots[[2]], hic.barplot)
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/wnt_pluri_combined.pdf", height = 9, width = 6)
cowplot::plot_grid(plotlist = plotlist, ncol = 1, labels = "AUTO")
dev.off()