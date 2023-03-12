require(automateR)

calculate_overlap <- function(deseq.one, deseq.two) {
    upreg.one <- splitDE(deseq.one)[[1]] %>% rownames()
    upreg.two <- splitDE(deseq.two)[[1]] %>% rownames()

    downreg.one <- splitDE(deseq.one)[[2]] %>% rownames()
    downreg.two <- splitDE(deseq.two)[[2]] %>% rownames()

    return(list(intersect(upreg.one, upreg.two), intersect(downreg.one, downreg.two)))
}

overlap_signif <- function(deseq.interest, deseq.pop, nPerm = 10000) {
    # Initialize experiment - define population, observed success states, total success states, gene list of interest
    de.overlap <- calculate_overlap(deseq.interest, deseq.pop)
    num.observed.suc <- length(de.overlap[[1]])
    num.total.suc <- length(rownames(splitDE(deseq.pop)[[1]]))
    upreg.interest <- rownames(splitDE(deseq.interest)[[1]])
    gene.pop <- rownames(deseq.pop)

    # Permute the population success states and compute overlap with our gene list of interest
    empty.vec <- rep(NA, nPerm)
    permute.df <- data.frame(overlap = empty.vec, greaterObs = empty.vec)
    for (i in 1:nPerm) {
        gene.sample <- sample(x = gene.pop, size = num.total.suc, FALSE)
        perm.overlap <- length(intersect(gene.sample, upreg.interest))
        permute.df$overlap[i] <- perm.overlap
        permute.df$greaterObs[i] <- (perm.overlap > num.observed.suc)
    }

    # Compute the expected overlap and p-value
    exp.overlap <- mean(permute.df$overlap)
    pval <- sum(permute.df$greaterObs) / nPerm

    # Instantiate overlap class with computed resutls and return
    olap.obj <- new("Overlap",
                    permutedf = permute.df,
                    observed = num.observed.suc,
                    expected = exp.overlap,
                    pval = pval)
    return(olap.obj)

}

mda.m1a.overlap <- overlap_signif(lung.vs.primary, m1a.o.comp, nPerm = 1e6)

mda.m1a.list <- list(MDA = lung.vs.primary, M1A = m1a.o.comp) %>%
                    lapply(function(de) {
                        dplyr::select(as.data.frame(de), log2FoldChange, padj)
                    })
merged.m1a.df <- merge(mda.m1a.list[[1]], mda.m1a.list[[2]], by = 0, all = TRUE) %>%
                    remove_rownames() %>%
                    column_to_rownames("Row.names")


make_table_up <- function(deseq.one, deseq.two) {
    upreg.one <- rownames(splitDE(deseq.one)[[1]])
    upreg.two <- rownames(splitDE(deseq.two)[[1]])
    notupreg.one <- rownames(deseq.one[! (rownames(deseq.one) %in% upreg.one), ])
    notupreg.two <- rownames(deseq.two[! (rownames(deseq.two) %in% upreg.two), ])

    both.upreg <- calculate_overlap(deseq.one, deseq.two)[[1]]
    one.not.two <- intersect(upreg.one, notupreg.two)
    two.not.one <- intersect(upreg.two, notupreg.one)
    not.either <- intersect(notupreg.one, notupreg.two)

    fisher.df <- data.frame(
        upreg_one = c(length(both.upreg), length(one.not.two)),
        not_upreg_one = c(length(two.not.one), length(not.either)),
        row.names = c("upreg_two", "not_upreg_two")
    )
    
    return(fisher.df)
}

mda.m1a.fisher <- make_table_up(lung.vs.primary, m1a.o.comp)

