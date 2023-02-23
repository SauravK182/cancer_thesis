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

res.log.tidy$pathway <- gsub(pattern = "_", replacement = " ", res.log.tidy$pathway)

# Bar plot
ggplot(res.log.tidy, aes(x = reorder(pathway, NES), y = NES)) +
    geom_col(aes(fill = padj < 0.05), color = "black") +
    coord_flip() +  # flip so NES is on x, y is the pathways
    labs(x = "Pathway", y = "Normalized Enrichment Score") +
    theme_cowplot()

# Let's try custom barplot to show p-value and NES - dotplot looks kind of ugly for single sample
require(scales)
res.log.tidy <- filter(res.log.tidy, padj < 0.05)
ggplot(res.log.tidy, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
    geom_col(aes(fill = padj)) +
    scale_fill_gradient2(position = "bottom",
                         low = "red",
                         mid = muted("magenta"),
                         high = "blue",
                         midpoint = median(res.log.tidy$padj)) +
    coord_flip()

#-----RUNNING GSEA WITH CLUSTERPROFILER-----
# Get hallmark gene set with msigdbr
require(msigdbr)
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene)

symbol.log <- capan.panc %>%
                as.data.frame() %>%
                rownames_to_column("Ensembl") %>%
                as_tibble() %>%
                dplyr::select(Ensembl, log2FoldChange) %>%
                na.omit() %>%  # remove NA symbols since we can't use these
                distinct() %>%
                group_by(Ensembl)

require(org.Hs.eg.db)
# clusterProfiler requires Entrez ID for GSEA/GO/KEGG analysis
entrez.id <- mapIds(org.Hs.eg.db,
                    keys = symbol.log$Ensembl,
                    keytype = "ENSEMBL",
                    column = "ENTREZID")
# Note after conversion to Entrez from HGNC symbol, we have 14,513 genes from a start of 17,501

symbol.entrez.log <- add_column(symbol.log, entrez = unname(entrez.id), .after = 1)[, -1] %>%
    na.omit()
# Going directly from Ensembl to Entrez we have 15396 - better to do this way then.

ranked.list.log <- deframe(symbol.entrez.log) %>%
    sort(decreasing = TRUE)

# By default, the GSEA function implements the fgsea algorithm
cp.gsea <- GSEA(ranked.list.log, TERM2GENE = hallmark)
head(cp.gsea)
dotplot(cp.gsea, showCategory = 20)    # dotplot is nice, but doesnt show direction of enrichment

require(DOSE)
panc.gsea <- capan.panc %>%
    prepare_gsea() %>%
    run_gsea_cp(category = "H", nPermSimple = 10000, seed = FALSE)

dotplot_enrichResult_Col(panc.gsea, showCategory = 20)

#----TRY THE COMPARECLUSTERFUNCTION
# Not sure if this works with the default GSEA function
gsea.list <- lapply(dge.list.full, prepare_gsea)
require(plyr)
test.compare <- compareCluster(geneClusters = gsea.list,
                               fun = "GSEA",
                               TERM2GENE = hallmark,
                               seed = TRUE)
# dotplot(test.compare, showCategory = 20)   # while all the below is important, updating enrichplot allows me to make a dotplot!
dotplot_compareClusterResult_gsea(test.compare, showCategory = 20)
# however, wish to display NES instead of geneRatio


# Keep getting "no enrichment found" - must be an error at lines 120 - 122 from the compareCluster.R file in the GitHub
# https://github.com/YuLab-SMU/clusterProfiler/blob/master/R/compareCluster.R

# Using the development version of the function seems to work - i.e., using from the GitHub (see below) does not produce error

# Need to update to the latest version of DOSE for this to work - https://github.com/YuLab-SMU/clusterProfiler/issues/409
# Unfortunately, can't install latest GitHub version since it requires HDO.db
# Can't install HDO.db from Bioconductor because my R is too out of date
# Going to try from GitHub (https://github.com/YuLab-SMU/HDO.db), cant, HDO.db requires R >= 4.2.0

require(GOSemSim)
# Need to update clusterProfiler entirely likely in order to use compareCluster with GSEA
# Most recent commits on Yu lab repos are for support for GSEA with clusterProfiler
test.out <- enrichplot::pairwise_termsim(test.compare)  # updating to most recent enrichplot version from GitHub seems to work!
test <- setReadable2(test.out, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")  # still doesnt work after updating enrichplot...





setReadable2 <- function(x, OrgDb, keyType="auto") {
    OrgDb <- load_OrgDb(OrgDb)
    if (!'SYMBOL' %in% columns(OrgDb)) {
        warning("Fail to convert input geneID to SYMBOL since no SYMBOL information available in the provided OrgDb...")
    }

    if (!(is(x, "enrichResult") || is(x, "groupGOResult") || is(x, "gseaResult") || is(x,"compareClusterResult")))
        stop("input should be an 'enrichResult' , 'gseaResult' or 'compareClusterResult' object...")

    isGSEA <- FALSE
    isCompare <- FALSE
    if (is(x, 'gseaResult'))
        isGSEA <- TRUE

    if (is(x, 'compareClusterResult'))
        isCompare <- TRUE

    if (keyType == "auto") {
        keyType <- x@keytype
        if (keyType == 'UNKNOWN') {
            stop("can't determine keyType automatically; need to set 'keyType' explicitly...")
        }
    }

    if (x@readable)
        return(x)

    gc <- geneInCategory(x)
    if (isGSEA) {
        genes <- names(x@geneList)
    } else if (isCompare) {
        if ("core_enrichment" %in% colnames(as.data.frame(x))) {
            geneslist <- x@geneClusters
            names(geneslist) <- NULL
            genes <- unique(names(unlist(geneslist)))
        } else {
            genes <- unique(unlist(x@geneClusters))
        }     
    } else {
        genes <- x@gene
    }

    gn <- EXTID2NAME(OrgDb, genes, keyType)


    if(isCompare) {
        gc2 <- list()
        k <- 1
        for(i in seq_len(length(gc))) {
            for(j in seq_len(length(gc[[i]]))) {
                gc2[[k]] <- gc[[i]][[j]]
                names(gc2)[k] <- paste(names(gc)[[i]], names(gc[[i]])[j], sep="-")
                k <- k + 1
            }
        }
        gc <- gc2
        gc <- lapply(gc, function(i) gn[i])
        res <- x@compareClusterResult
        gc <- gc[paste(res$Cluster, res$ID, sep= "-")]
    } else {
        gc <- lapply(gc, function(i) gn[i])
        res <- x@result
        gc <- gc[as.character(res$ID)]
    }

    ## names(gc) should be identical to res$ID

    ## gc <- gc[as.character(res$ID)]


    geneID <- sapply(gc, paste0, collapse="/")
    # if (isGSEA) {
    if ("core_enrichment" %in% colnames(as.data.frame(x))) {
        res$core_enrichment <- unlist(geneID)
    } else {
        res$geneID <- unlist(geneID)
    }
    x@gene2Symbol <- gn
    x@keytype <- keyType
    x@readable <- TRUE
    if(isCompare){
        x@compareClusterResult <- res
    } else {
        x@result <- res
    }


    return(x)
}

#-----TRY TO REDEFINE DOTPLOT FUNCTION-----
# ep_str_wrap and default_labeller from utilities.R in enrichplot github: https://rdrr.io/bioc/enrichplot/src/R/utilities.R#sym-ep_str_wrap






#----TRY TO ADAPT TO SHOW MULTIPLE GSEA RESULTS AT ONCE

