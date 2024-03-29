# Set up GSEA list
gsea.list <- lapply(dge.list.full, prepare_gsea)
names(gsea.list) <- names.comp

# HALLMARK GSEA RESULTS
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene)
hallmark.compare <- automateR::compareCluster(geneClusters = gsea.list,
                                   fun = "GSEA",
                                   TERM2GENE = hallmark,
                                   seed = TRUE)
dotplot_compareClusterResult_gsea(hallmark.compare, showCategory = 50)

#------Save hallmark gsea results--------
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/rna_gsea_hallmark.pdf",
          width = 9,
          height = 10)
dotplot_compareClusterResult_gsea(hallmark.compare, showCategory = 50) +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 9,
                                     angle = 30,
                                     hjust = 1,
                                     vjust = 1))
dev.off()

# CELL TYPE GSEA RESULTS
cell.type <- msigdbr(species = "Homo sapiens", category = "C8") %>%
    dplyr::select(gs_name, entrez_gene)
cell.type.compare <- automateR::compareCluster(geneClusters = gsea.list,
                                               fun = "GSEA",
                                               TERM2GENE = cell.type,
                                               seed = TRUE)
dotplot_compareClusterResult_gsea(cell.type.compare, showCategory = 10)


# IMMUNLOGICAL GSEA RESULTS
immuno <- msigdbr(species = "Homo sapiens", category = "C7") %>%
    dplyr::select(gs_name, entrez_gene)
immuno.compare <- automateR::compareCluster(geneClusters = gsea.list,
                                               fun = "GSEA",
                                               TERM2GENE = immuno,
                                               seed = TRUE)
dotplot_compareClusterResult_gsea(immuno.compare, showCategory = 20)


# ONCOGENIC SIGNATURE GSEA RESULTS
onco <- msigdbr(species = "Homo sapiens", category = "C6") %>%
    dplyr::select(gs_name, entrez_gene)
onco.compare <- automateR::compareCluster(geneClusters = gsea.list,
                                          fun = "GSEA",
                                          TERM2GENE = onco,
                                          seed = TRUE)
dotplot_compareClusterResult_gsea(onco.compare, showCategory = 20)