# Calculate a significance level for the intersection of the upreg genes from ccRCC
intersect.genes <- intersect(rownames(splitDE(m1a.o.comp)[[1]]),
                             rownames(splitDE(lm.rc.comp)[[1]]))

mb.786.int <- intersect(rownames(signifDE(m1a.o.comp)), rownames(signifDE(lung.vs.primary)))
numTot = 16909
numSucState = 2145
geneListSize = 3039
observedOverlap = length(intersect.genes)

# Checking to see how many of these genes are near changes in H3K27ac in both samples
require(ChIPpeakAnno)
test.anno <- annotatePeakInBatch(dba.report(rod.chip[[3]][[1]]),        # results for M1A vs. O
                                 AnnotationData = TSS.human.GRCh37,     # seems to work with the TSS GRCh37 data
                                 output = "nearestBiDirectionalPromoter",
                                 bindingRegion = c(-2000, 2000))
sum(intersect.genes %in% test.anno$feature)     # only 51 genes proximal to change in H3K27ac in M1A

test.anno2 <- annotatePeakInBatch(dba.report(rod.chip[[3]][[2]]),        # results for LM1 vs. RC2
                                 AnnotationData = TSS.human.GRCh37,     # seems to work with the TSS GRCh37 data
                                 output = "nearestBiDirectionalPromoter",
                                 bindingRegion = c(-2000, 2000))
sum(intersect.genes %in% test.anno2$feature)    # 28 genes proximal to changes in H3K27ac in LM1

# Calculate p-value with hypergeometric dist
# pval <- 0
# for (i in 0:(observedOverlap - 1)) {
#     marginal <- (choose(numSucState, i) * choose(numTot - numSucState, geneListSize - i)) / choose(numTot, geneListSize)
#     pval <- pval + marginal
# }
# pval <- 1 - pval

# Calculate p-value with Fisher's exact test
fisher.df <- data.frame(genes.not.interest = c(2145 - 738, 16909 - 3039 - (2145 - 738)),
                        genes.in.interest = c(738, 3039 - 738))
rownames(fisher.df) <- c("In_category", "Not_in_category")
fisher.test(fisher.df, alternative = "greater")

# Calculate empirical p-value by randomly sampling gene sets of size 2145
numPerms <- 10000
permOverlap <- c()
for (i in 1:numPerms) {
    rand.set <- sample(rownames(lm.rc.comp), size = numSucState)
    num.intersect <- length(intersect(rand.set, rownames(splitDE(m1a.o.comp)[[1]])))
    permOverlap <- c(permOverlap, num.intersect)
}

# Plot histogram
overlap.df <- data.frame(size = permOverlap, greater = permOverlap > observedOverlap)
ggplot(data = overlap.df, aes(x = size, fill = greater)) +
    geom_histogram(color = "black")
