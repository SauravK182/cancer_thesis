# Test with ChIPseeker
require(ChIPseeker)
require(EnsDb.Hsapiens.v75)
# Note, an optional argument annoDb can be passed, which will be used to annotate for gene names
test.seek.anno <- annotatePeak(dba.report(rod.chip[[3]][[1]]),
                               tssRegion = c(-2000, 2000),
                               TxDb = EnsDb.Hsapiens.v75)   # TxDb takes TxDb or EnsDb object
as.GRanges(test.seek.anno)
# Note that ChIPseeker will report the gene with the nearest TSS for all peaks
# However, peaks outside of the desired range will be labeled as "distal intergenic"
sum(abs(as.GRanges(test.seek.anno)$distanceToTSS) < 2000)   # seems there are 3423 TSS-proximal peaks

# Make some visualizations
plotAnnoPie(test.seek.anno)
plotAnnoBar(test.seek.anno)
plotDistToTSS(test.seek.anno)
upsetplot(test.seek.anno, vennpie = TRUE)

# Test with ChIPpeakAnno
require(ChIPpeakAnno)
test.anno <- annotatePeakInBatch(dba.report(rod.chip[[3]][[2]]),        # results for M1A vs. O
                                 AnnotationData = TSS.human.GRCh37,     # seems to work with the TSS GRCh37 data
                                 output = "nearestBiDirectionalPromoter",
                                 bindingRegion = c(-2000, 2000))

# note that some peaks are duplicated - overlap multiple TSS's
# this is due to the fact that the output is nearestBiDirectionalPromoter
# here, strand is considered, and the function will report bidirectional promoters if there are promoters in both
# directions of the given region
require(org.Hs.eg.db)
test.anno.up <- mapIds(org.Hs.eg.db,
                       keys = test.anno[test.anno$Fold > 0, ]$feature,
                       keytype = "ENSEMBL",
                       column = "ENTREZID",
                       multiVals = "first")

test.anno.down <- mapIds(org.Hs.eg.db,
                         keys = test.anno[test.anno$Fold < 0, ]$feature,
                         keytype = "ENSEMBL",
                         column = "ENTREZID",
                         multiVals = "first")

# Use enrichKEGG to get pathway analysis results
# Note that we can perform ORA or GSEA with KEGG gene sets
require(clusterProfiler)
# Originally trying to run enrichKEGG gave errors about failure to download data
# Try setting download method to auto with R.utils - https://github.com/YuLab-SMU/clusterProfiler/issues/305
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
# keep getting no hits with test.anno.up
rod.786.kegg <- enrichKEGG(gene = unname(test.anno.up[!is.na(test.anno.up)]), # genes need to be in EntrezID!
                           pvalueCutoff = 0.05,
                           keyType = "ncbi-geneid",
                           pAdjustMethod = "BH")

rod.786.kdown <- enrichKEGG(gene = test.anno.down)
dotplot(rod.786.kdown)

# Try GO term analysis - also null
test.anno.up <- test.anno[test.anno$Fold < 0, ]$feature
rod.786.go <- clusterProfiler::enrichGO(gene = test.anno.up,
                       keyType = "ENSEMBL",
                       pvalueCutoff = 0.1,
                       OrgDb = org.Hs.eg.db,
                       ont = "MF")
