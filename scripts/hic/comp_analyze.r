# Read in diff comp files
get_comp <- function(name, file = "filtered", getsubcomp = TRUE) {
    compdir <- file.path("E:/SK/data/", name, "dchic_res")
    comp <- file.path(compdir, "40kb/DifferentialResult/sample/fdr_result/")
    if (file == "full") {
        compfile <- file.path(comp, "differential.intra_sample_combined.pcQnm.bedGraph")
    } else if (file == "filtered") {
        compfile <- file.path(comp, "differential.intra_sample_group.Filtered.pcQnm.bedGraph")
    } else {
        stop("Variable 'file' must be either 'filtered' or 'raw' to read in associated dcHiC bedGraph file.")
    }
    subcompfile <- file.path(comp, "intra_sample_group.subcompartments.bedGraph")

    tryCatch(
        {
            message(paste0("Reading in from ", compfile, "..."))
            diffcomp <- read.delim(compfile)
            message("Done!")

            if (getsubcomp) {
                message(paste0("Reading in from ", subcompfile, "..."))
                subcomp <- read.delim(subcompfile)
                message("Done!")
                return(list(diffcomp = diffcomp, subcomp = subcomp))
            }
        },
        error = function(e) {
            stop("Differential compartment score file not found. Please try again")
        }
    )

    return(diffcomp)

}

# To assign either within or between compartment transitions
compswitch <- function(compdf, refcol = 5, expcol = 4) {
    aref <- compdf[, refcol] > 0
    aexp <- compdf[, expcol] > 0

    # Use mutate with case_when to add a column for compartment transition in vectorized format
    compdf.mut <- compdf %>%
                    mutate(transition = case_when(
                        aref & aexp ~ "Within A Transition",
                        aref & !aexp ~ "A to B Transition",
                        !aref & aexp ~ "B to A Transition",
                        !aref & !aexp ~ "Within B Transition"
                    ))
    
    return(compdf.mut)
}

anno_hic <- function(hic.list) {
    hic.gr.list <- lapply(hic.list, function(df) {
    df %>%
        compswitch() %>%
        (function(x) {
            for (i in seq_len(nrow(x))) {
                x$chr[i] <- gsub(pattern = "chr", replacement = "", x = x$chr[i]) # convert to GRCh37 chr names
            }
            return(x)
        }) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    })

    # Annotate diff comp regions, reporting all genes overlapping each region
    hic.anno.list <- lapply(hic.gr.list, function(gr) {
        gr.anno <- anno_peak_gr37(gr, output = "overlapping", region = NULL, select = "all")
        gr.omit <- gr.anno[!is.na(gr.anno$feature), ]
    })

    return(hic.anno.list)
}

panc.comp <- get_comp("ren-panc")
rcc.comp <- get_comp("rod-rcc", getsubcomp = FALSE)

# Assign transition states to compartment scores
# Note diff compartmentalization is called at FDR 0.10 threshold
hic.list <- list(panc = panc.comp[["diffcomp"]], m1a_o = rcc.comp)

# Get LFC for genes in each of the 4 subgroups
hic.anno.list <- anno_hic(hic.list)
genexp.comp <- lapply(names(hic.anno.list), function(name) {
    splitcomp <- split(as.data.frame(hic.anno.list[[name]]), hic.anno.list[[name]]$transition)
    comp.exp.df <- lapply(splitcomp, function(df) {
        features <- unique(df$feature)
        lfc <- as.data.frame(dge.list.full[[name]])[features, "log2FoldChange"]
        return(data.frame(row.names = features,
                          lfc = lfc,
                          transition = as.factor(df$transition[seq_len(length(features))])))
    }) %>%
        purrr::reduce(rbind)
})

# Plot gene expression by compartment change
names.comp <- c("Pancreatic Cell Line", "786 Renal Cell Line")
col.vec4 <- c("blue", "red", "gray", "gray")
compexp.plots <- lapply(seq_len(length(genexp.comp)), function(i) {
    comp.list <- list(c("A to B Transition", "B to A Transition"),
                      c("Within A Transition", "B to A Transition"),
                      c("Within B Transition", "B to A Transition"))
    ggplot(data = genexp.comp[[i]], aes(x = transition, y = lfc, fill = transition)) +
        geom_boxplot(color = "black", width = 0.4) +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_violin(alpha = 0.4) +
        geom_signif(map_signif_level = c("***" = 0.001 / length(comp.list),
                                         "**" = 0.01 / length(comp.list),
                                         "*" = 0.05 / length(comp.list)),
                    comparisons = comp.list,
                    test = "wilcox.test",
                    test.args = list(alternative = "less", paired = FALSE),
                    step_increase = 0.1,
                    size = 0.4,
                    textsize = 5,
                    vjust = 0.5
                    ) +
        ylab("Log2 FC Gene Expression") +
        theme_cowplot() +
        theme(axis.title.x = element_blank(),
              axis.text.y = element_text(size = 9),
              plot.title = element_text(size = 9),
              axis.text.x = element_text(size = 8),
              legend.position = "none") +
        ggtitle(names.comp[i]) +
        scale_fill_manual(values = col.vec4)
})

# Save plot
bp.genexp.comp <- plot_grid(plotlist = compexp.plots, labels = "AUTO", ncol = 1)
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/compartment_expr.pdf", height = 11, width = 8)
bp.genexp.comp
dev.off()

# Get diff in comp scores and LFC gene exp
panc.full <- get_comp("ren-panc", file = "full", getsubcomp = FALSE)
rcc.full <- get_comp("rod-rcc", file = "full", getsubcomp = FALSE)
hic.full <- list(panc = panc.full, m1a_o = rcc.full)
hic.anno.full <- anno_hic(hic.full)

# NOTE - many genes overlap multiple windows
# In some cases, one window is changing compartment, the other isnt
# Important to distinguish this in the results - do not label compartment transitions, since we get duplicate genes
# Instead, just rely on the average change in compartment score for a gene (metastasis - primary)
dcomp.lfc <- lapply(names(hic.anno.full), function(name) {
    df <- as.data.frame(hic.anno.full[[name]])
    colnames(df)[c(10, 11)] <- c("metastasis", "primary")
    features.dcomp <- df %>%
                            select(metastasis, primary, feature) %>%
                            dplyr::group_by(feature) %>%
                            dplyr::summarize(dComp = mean(metastasis - primary))
    lfc <- as.data.frame(dge.list.full[[name]])[features.dcomp$feature, "log2FoldChange"]
    dcomp <- data.frame(row.names = features.dcomp$feature,
                         genelfc = lfc,
                         dComp = features.dcomp$dComp)
})

# Create bivariate density plot of LFC vs dComp
# Note: for drawing panel around ggplot object, see https://stackoverflow.com/questions/26191833/
# Use panel.grid options to blank out background, draw a rectangular border
names.vec <- c("Pancreatic System", "786 ccRCC System")

# Need to hardcore p-values since dynamic expressions within lapply() are not working
# %*% renders as times symbol in plotmath
pval.panc <- expression(
    "p = 1.7" %*%~10^{-185}
)
pval.rod <- expression(
    "p = 4.7" %*%~10^{-94}
)
pval.list <- list(pval.panc, pval.rod)

xlabel <- expression(Delta~Eigenvector~Signal)
pdensity.list <- lapply(seq_len(length(dcomp.lfc)), function(i) {
    pval <- cor.test(dcomp.lfc[[i]]$genelfc, dcomp.lfc[[i]]$dComp)
    pearson <- paste("r =", formatC(pval$estimate, digits = 2, format = "f"))
    # Note: trying to format p-value using plotmath and expression-type objects doesn't work
    # Either doesn't render or errors out when trying to parse
    # pval.format <- formatC(pval$p.value, digits = 1, format = "e")
    # pval.nums <- strsplit(pval.format, "e")
    # pval.final <- paste("'p ='~", pval.nums[[1]][1], "~\u00D7~'10'^", pval.nums[[1]][2])
    density.plot <- ggplot(data = dcomp.lfc[[i]], aes(x = dComp, y = genelfc)) +
                        stat_density_2d(aes(fill = ..level..), geom = "polygon", color = "white") +
                        theme_bw() +
                        scale_fill_viridis_c() +
                        xlab(xlabel) +
                        ylab("Log2 FC Gene Expression") +
                        ggtitle(names.vec[i]) +
                        annotate("text",
                                 x = rep(quantile(dcomp.lfc[[i]]$dComp, 0.08, na.rm = TRUE), 2),
                                 y = c(quantile(dcomp.lfc[[i]]$genelfc, 0.9, na.rm = TRUE) + 0.1,
                                       quantile(dcomp.lfc[[i]]$genelfc, 0.89, na.rm = TRUE)),
                                 label = c(pearson, pval.list[[i]])) +
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_rect(color = "black", linewidth = 1),   # draw border
                              plot.title = element_text(face = "bold"),     # boldface title
                              axis.title.x = element_text(size = 13))
})

# Save plot
density.grid <- cowplot::plot_grid(plotlist = pdensity.list, nrow = 2, ncol = 1, labels = "AUTO")
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/comp_genelfc_density.pdf", height = 10, width = 8)
density.grid
dev.off()


# Pie chart to visualize % of windows changing compartment
percent.change <- lapply(names(hic.list), function(name) {
    total.windows <- nrow(hic.full[[name]])
    df.comp <- compswitch(hic.list[[name]])
    a.to.b <- sum(df.comp$transition == "A to B Transition")
    b.to.a <- sum(df.comp$transition == "B to A Transition")
    stable <- total.windows - (a.to.b + b.to.a)
    return(c(a.to.b, b.to.a, stable) / total.windows)
})
