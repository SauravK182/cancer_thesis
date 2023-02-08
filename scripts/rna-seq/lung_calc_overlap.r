require(automateR)
setClass("Overlap", representation = representation(
        permutedf = "data.frame",
        observed = "numeric",
        expected = "numeric",
        pval = "numeric",
        type = "character"
    ))

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
        gene.sample <- sample(x = gene.pop, size = num.total.suc)
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

mda.m1a.overlap <- overlap_signif(lung.vs.primary, m1a.o.comp)
