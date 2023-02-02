require(fgsea)
require(clusterProfiler)
require(cowplot)

# GSEA with the fgsea package
# First, isolate only Wald statistic for ranking and symbol
symbol.stat <- capan.panc %>%
                as_tibble() %>%
                dplyr::select(Symbol, log2FoldChange) %>%
                na.omit() %>%  # remove NA symbols since we can't use these
                distinct() %>%
                group_by(Symbol) %>%
                summarize(stat = mean(stat))
# In order to deal with the fact that one symbol may have multiple stats, average the stats

# Use deframe to convert 2-col df into named list - first col is names
ranked.list <- deframe(symbol.stat)

# Load in hallmark gene set
hallmark.sets <- gmtPathways("gsea/h.all.v2022.1.Hs.symbols.gmt")

# Run fgsea with 10000 permutations for initial p-value estimation
fgsea.res <- fgsea(pathways = hallmark.sets,
                   stats = ranked.list,
                   nPermSimple = 10000)

# Tidy up results
fgsea.res.tidy <- fgsea.res %>%
    as_tibble() %>%
    arrange(desc(NES)) # arranges in descending order by NES

# Plot NES
ggplot(fgsea.res.tidy, aes(x = reorder(pathway, NES), y = NES)) +
    geom_col(aes(fill = padj < 0.05), color = "black") +
    coord_flip() +  # flip so NES is on x, y is the pathways
    labs(x = "Pathway", y = "Normalized Enrichment Score") +
    theme_cowplot()


#----RANKING WITH LFC-----
symbol.log <- capan.panc %>%
                as_tibble() %>%
                dplyr::select(Symbol, log2FoldChange) %>%
                na.omit() %>%  # remove NA symbols since we can't use these
                distinct() %>%
                group_by(Symbol)

ranked.list.log <- deframe(symbol.log)
fgsea.res.log <- fgsea(pathways = hallmark.sets,
                       stats = ranked.list.log)

res.log.tidy <- fgsea.res.log %>%
    as_tibble() %>%
    arrange(desc(NES))

ggplot(res.log.tidy, aes(x = reorder(pathway, NES), y = NES)) +
    geom_col(aes(fill = padj < 0.05), color = "black") +
    coord_flip() +  # flip so NES is on x, y is the pathways
    labs(x = "Pathway", y = "Normalized Enrichment Score") +
    theme_cowplot()


#-----RUNNING GSEA WITH CLUSTERPROFILER-----
# Get hallmark gene set with msigdbr
require(msigdbr)
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene)

symbol.log <- capan.panc %>%
                as_tibble() %>%
                dplyr::select(Symbol, log2FoldChange) %>%
                na.omit() %>%  # remove NA symbols since we can't use these
                distinct() %>%
                group_by(Symbol)

require(org.Hs.eg.db)
# clusterProfiler requires Entrez ID for GSEA/GO/KEGG analysis
entrez.id <- mapIds(org.Hs.eg.db,
                    keys = symbol.log$Symbol,
                    keytype = "SYMBOL",
                    column = "ENTREZID")
# Note after conversion to Entrez, we have 14,513 genes from a start of 17,501

symbol.entrez.log <- add_column(symbol.log, entrez = unname(entrez.id), .after = 1)[, -1] %>%
    na.omit()

ranked.list.log <- deframe(symbol.entrez.log) %>%
    sort(decreasing = TRUE)

# By default, the GSEA function implements the fgsea algorithm
cp.gsea <- GSEA(ranked.list.log, TERM2GENE = hallmark)
head(cp.gsea)
dotplot(cp.gsea)    # dotplot is nice, but doesnt show direction of enrichment
