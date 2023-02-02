require(edgeR)
require(Rsamtools)

# Get RPKM values for genes
## We're going to see if genes near changes in H3K27ac occupancy
## Have higher average RPKM compared to in the primary tumor

# First, check with LFC values
require(ChIPpeakAnno)
require(EnsDb.Hsapiens.v75)
test.anno <- annotatePeakInBatch(dba.report(rod.chip[[3]][[1]]),        # results for M1A vs. O
                                 AnnotationData = TSS.human.GRCh37,     # seems to work with the TSS GRCh37 data
                                 output = "nearestBiDirectionalPromoter",
                                 bindingRegion = c(-2000, 2000))
length(unique(test.anno$feature))   # 1831 unique features here

# In total, there are 1852 regions within 2kb of a TSS out of 13205 total
test.anno.up <- test.anno[test.anno$Fold > 0, ]
test.anno.down <- test.anno[test.anno$Fold < 0, ]

# Isolate features of each
up.h3.features <- test.anno.up$feature
down.h3.features <- test.anno.down$feature

# Subset the DESeq2 results table to include only these features
# Apparently, some features are missing, so for preliminary analysis, let's use non-missing features
up.h3.features <- intersect(up.h3.features, rownames(m1a.o.comp))
down.h3.features <- intersect(down.h3.features, rownames(m1a.o.comp))

# Now subset
genes.prox.h3up <- m1a.o.comp[up.h3.features, ]
genes.prox.h3down <- m1a.o.comp[down.h3.features, ]

# Boxplot with significance
upreg.vec <- rep("upreg", length(up.h3.features))
downreg.vec <- rep("downreg", length(down.h3.features))
data.df <- data.frame(
                      LogFC = c(genes.prox.h3up$log2FoldChange,
                                genes.prox.h3down$log2FoldChange),
                      DiffExp = c(upreg.vec, downreg.vec)
)
require(ggsignif)
require(cowplot)
genes.bplot <- ggplot(data = data.df, aes(x = DiffExp, y = LogFC)) +
                geom_boxplot(color = "black", aes(fill = as.factor(DiffExp))) +
                geom_signif(comparisons = list(c("upreg", "downreg")),
                            map_signif_level = TRUE,
                            test = "wilcox.test",
                            test.args = list(alternative = "two.sided", paired = FALSE)) +
                theme_cowplot(font_size = 14, line_size = 1) +
                xlab("Change in H3K27ac occupancy proximal to TSS")
genes.bplot
