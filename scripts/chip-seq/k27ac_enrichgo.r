k27.genes <- lapply(seq_len(length(anno.chip.list.full)), function(i) {
    gr <- anno.chip.list.full[[i]]
    return(gr[gr$Fold > 0, ]$feature)
})
# Need to name list for compareCluster to run properly
names(k27.genes) <- names.comp

k27.entrez <- lapply(k27.genes, function(x) {
    mapIds(org.Hs.eg.db,
           keys = x,
           keytype = "ENSEMBL",
           column = "ENTREZID") %>% unname()
})

# Run compareCluster, use dotplot to visualize results
# Note: can pass ENSEMBL IDs to enrichGO, just need to indicate they are ENSEMBL with keyType
# and provide org.Hs.eg.db as OrgDb so it can be internally mapped to Entrez ID
k27.go <- automateR::compareCluster(geneClusters = k27.genes,
                                    fun = "enrichGO",
                                    keyType = "ENSEMBL",
                                    OrgDb = org.Hs.eg.db,
                                    pvalueCutoff = 0.1,
                                    ont = "MF")

cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/k27_go.pdf", height = 10, width = 8)
clusterProfiler::dotplot(k27.go, showCategory = 10) +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9))
dev.off()








# Originally trying to run enrichKEGG gave errors about failure to download data
# Try setting download method to auto with R.utils - https://github.com/YuLab-SMU/clusterProfiler/issues/305
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
k27.kegg <- automateR::compareCluster(geneClusters = k27.entrez,
                                      fun = "enrichKEGG",
                                      organism = "hsa",
                                      keyType = "ncbi-geneid",
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)