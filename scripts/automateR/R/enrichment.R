#' @title Prepare `DESeq` results for GSEA
#'
#' @description
#' Will automatically create a ranked list of genes with proper identifiers from a `DESeq2 results` object for
#' either `clusterProfiler` or `fgsea`. Utilizes org.Hs.eg.db for any necessary gene identifier conversions.
#'
#' @param deseqObject A `results` object from `DESeq2`, or any other valid data-frame like object that possesses
#' a numerical ranking score and gene identifers
#' @param rankBy Column name in the data frame to rank genes by. Default = "log2FoldChange".
#' @param colSymbol Column name for the column in the data frame containing the gene identifiers.
#' Note that these identifiers must be a valid keytype in the org.Hs.eg database. Default = "rownames" (i.e., gene
#' identifiers are the rownames of the object).
#' @param useSymbol Keytype for the gene identifiers housed in the column `colSymbol`. Must be a valid
#' keytype that can be queried by org.Hs.eg.db. Default = "ENSEMBL"
#' @param method Either `"cluterProfiler"` or `"fgsea"`, indicating which package will be used for downstream analysis.
#' Default = "clusterProfiler"
#' @param convertSymbol Boolean value indicating whether to convert symbols to the appropriate type for `clusterProfiler`
#' or `fgsea`. Default = TRUE
#'
#' @return Named numeric vector sorted from highest to lowest `rankBy` stat. Names correspond to the gene IDs given the keytype
#' required by `clusterProfiler` (Entrez ID) or `fgsea` (HGNC Symbol), respectively.
#' @export
#'
#' @examples
prepare_gsea <- function(deseqObject,
                         rankBy = "log2FoldChange",
                         colSymbol = "rownames",
                         useSymbol = "ENSEMBL",
                         method = "clusterProfiler",
                         convertSymbol = TRUE) {
    # Ensure the rankBy variable is a column name of the df
    if (! (rankBy %in% colnames(as.data.frame(deseqObject))) | (class(as.data.frame(deseqObject)[, rankBy]) != "numeric")) {
        stop("Parameter rankBy must point to a numeric column in deseqObject for GSEA ranking.")
    }

    # If user passed "ensembl", make a column from the ensembl row; else, save the colname
    if (tolower(colSymbol) == "rownames") {
        symbol.col <- useSymbol
        deseqObject <- deseqObject %>%
                        as.data.frame() %>%
                        rownames_to_column(symbol.col)
    } else if (colSymbol %in% colnames(deseqObject)) {
        symbol.col <- useSymbol
    } else {
        stop("useSymbol must either be set to 'ensembl' or a column name containing gene symbols for annotation.")
    }

    # Organize deseqObject to obtain key information required for GSEA
    symbol.stat <- deseqObject %>%
                    as_tibble() %>%
                    dplyr::select(symbol.col, rankBy) %>%
                    na.omit() %>%
                    distinct() %>%
                    group_by(!!! syms(symbol.col))

    # If user desires to convert symbol, do so - clusterProfiler requires Entrez ID
    if (tolower(method) == "clusterprofiler") {
        if (tolower(useSymbol) == "ensembl") {
            keytype <- "ENSEMBL"
        } else if (useSymbol %in% keytypes(org.Hs.eg.db)) {
            keytype <- useSymbol
        } else {
            stop("Must provide a valid keytype for use in org.Hs.eg.db with parameter useSymbol.")
        }

        if (convertSymbol) {
            entrez.id <- mapIds(org.Hs.eg.db,
                                keys = symbol.stat %>% pull(all_of(symbol.col)),
                                keytype = keytype,
                                column = "ENTREZID")
            symbol.stat <- add_column(symbol.stat, entrez = unname(entrez.id), .after = 1)[, -1]
        }
    } else if (tolower(method) == "fgsea") {
        if (tolower(useSymbol) == "ensembl") {
            keytype <- "ENSEMBL"
        } else if (useSymbol %in% keytypes(org.Hs.eg.db)) {
            keytype <- useSymbol
        } else {
            stop("Must provide a valid keytype for use in org.Hs.eg.db with parameter useSymbol.")
        }

        if (convertSymbol) {
            entrez.id <- mapIds(org.Hs.eg.db,
                                keys = symbol.stat %>% pull(all_of(symbol.col)),
                                keytype = keytype,
                                column = "SYMBOL")
            symbol.stat <- add_column(symbol.stat, entrez = unname(entrez.id), .after = 1)[, -1]
        }
    }
    ranked.list <- symbol.stat %>%
                        na.omit() %>%
                        deframe() %>%
                        sort(decreasing = TRUE)
    return(ranked.list)
}

#' @title Run `clusterProfiler` GSEA
#' 
#' @description
#' Simple wrapper function for the `GSEA` function in the `clusterProfiler` package. Utilizes `msigdbr` to read in
#' a hallmark set and run `GSEA` with a pre-ranked list (which can be constructed with the `prepare_gsea()` function).
#'
#' @param ranked.list Ranked list in the form of a named numeric vector, sorted in descending order with Entrez ID identifiers.
#' @param category Category from MSigDB to use by `GSEA()`. Will be fetched using `msigdbr()`.
#' @param ... Additional arguments to pass to `GSEA()`.
#'
#' @return `gseaResult` object containing the results of the call to `GSEA()`.
#' @export
#'
#' @examples
run_gsea_cp <- function(ranked.list, category, ...) {
    # Load in desired gene set with msigdbr
    gene.set <- msigdbr(species = "Homo sapiens", category = category) %>%
        dplyr::select(gs_name, entrez_gene)

    # Run GSEA with desired arguments
    gsea <- GSEA(ranked.list, TERM2GENE = gene.set, ...)
    return(gsea)
}

# FROM YU ET AL.
# Function included here from Yu lab source code (https://github.com/YuLab-SMU/clusterProfiler/blob/master/R/compareCluster.R)
# Due to issues with package/R versions

#' @title Yu Lab's compareCluster
#' 
#' @description compareCluster function from the Yu Lab (https://github.com/YuLab-SMU/clusterProfiler/blob/master/R/compareCluster.R).
#' Imported in this package due to versioning issues. See the above link for more information. All credit to Yu et al.
#' 
#' @export
compareCluster <- function(geneClusters,
                           fun="enrichGO", data='',
                           source_from=NULL, ...) {

   if(is.character(fun)){
     if(fun %in% c("groupGO", "enrichGO", "enrichKEGG",
                   "gseGO", "gseKEGG", "GSEA", "gseWP")){
       fun <- utils::getFromNamespace(fun, "clusterProfiler")
     } else if(fun %in% c("enrichDO", "enrichDGN", "enrichDGNv",
                          "enrichNCG", "gseDO", "gseNCG", "gseDGN")){
       fun <- utils::getFromNamespace(fun , "DOSE")
     } else if(fun %in% c("enrichPathway", "gsePathway")){
        fun <- utils::getFromNamespace(fun , "ReactomePA")
     } else if(fun %in% c("enrichMeSH", "gseMeSH")){
        fun <- utils::getFromNamespace(fun , "meshes")
     } else {
       source_env <- .GlobalEnv
       if(!is.null(source_from)){
         source_env <- loadNamespace(source_from)
       }
       # If fun is in global or any loaded package, this will get it
       # This assumes that a user will actually load said package.
       fun <- get(fun, envir = source_env)
     }

   }


    # Use formula interface for compareCluster
    if (typeof(geneClusters) == 'language') {
        if (!is.data.frame(data)) {
            stop ('no data provided with formula for compareCluster')
        } else {
            genes.var = all.vars(geneClusters)[1]
            n.var = length(all.vars(geneClusters))
            # For formulas like x~y+z
            grouping.formula = gsub('^.*~', '~',
                       as.character(as.expression(geneClusters)))
            n.group.var = length(all.vars(formula(grouping.formula)))
            geneClusters = dlply(.data=data, formula(grouping.formula),
                                 .fun=function(x) {
                if ( (n.var - n.group.var) == 1 ) {
                    as.character(x[[genes.var]])
                } else if ( (n.var - n.group.var) == 2 ) {
                    fc.var = all.vars(geneClusters)[2]
                    geneList = structure(x[[fc.var]], names = x[[genes.var]])
                    sort(geneList, decreasing=TRUE)
                } else {
          stop('only Entrez~group or Entrez|logFC~group type formula is supported')
                }
            })
        }
    }
    clProf <- llply(geneClusters,
                    .fun=function(i) {
                        x=suppressMessages(fun(i, ...))

                        if (inherits(x, c("enrichResult",
                                          "groupGOResult", "gseaResult"))){
                            as.data.frame(x)
                        }
                    }
                    )
    clusters.levels = names(geneClusters)
    clProf.df <- ldply(clProf, rbind)

    if (nrow(clProf.df) == 0) {
        warning("No enrichment found in any of gene cluster, please check your input...")
        return(NULL)
    }

    #clProf.df <- dplyr::rename(clProf.df, c(.id="Cluster"))
    clProf.df <- plyr::rename(clProf.df, c(.id="Cluster"))
    clProf.df$Cluster = factor(clProf.df$Cluster, levels=clusters.levels)

    if (is.data.frame(data) && grepl('+', grouping.formula)) {
        groupVarName <- strsplit(grouping.formula, split="\\+") %>% unlist %>%
            gsub("~", "", .) %>% gsub("^\\s*", "", .) %>% gsub("\\s*$", "", .)
        groupVars <- sapply(as.character(clProf.df$Cluster),
                            strsplit, split="\\.") %>% do.call(rbind, .)
        for (i in seq_along(groupVarName)) {
            clProf.df[, groupVarName[i]] <- groupVars[,i]
        }
        i <- which(colnames(clProf.df) %in% groupVarName)
        j <- (1:ncol(clProf.df))[-c(1, i)]
        clProf.df <- clProf.df[, c(1, i, j)]
    }

    ##colnames(clProf.df)[1] <- "Cluster"
    res <- new("compareClusterResult",
               compareClusterResult = clProf.df,
               geneClusters = geneClusters,
               .call = match.call(expand.dots=TRUE)
               )

    params <- modifyList(extract_params_yu(args(fun)),
                         extract_params_yu(res@.call))

    keytype <- params[['keyType']]
    if (is.null(keytype)) keytype <- "UNKNOWN"
    readable <- params[['readable']]
    if (length(readable) == 0) readable <- FALSE

    res@keytype <- keytype
    res@readable <- as.logical(readable)
    ## work-around for bug in extract_parameters -- it doesn't match default args
    res@fun <- params[['fun']] %||% 'enrichGO'

    return(res)
}

#' @export
extract_params_yu <- function(x) {
    y <- rlang::quo_text(x)
    if (is.function(x)) y <- sub('\nNULL$', '', y)

    y <- gsub('"', '', y) %>%
        ## sub(".*\\(", "", .) %>%
        sub("[^\\(]+\\(", "", .) %>%
        sub("\\)$", "", .) %>%
        gsub("\\s+", "", .)

    y <- strsplit(y, ",")[[1]]
    params <- sub("=.*", "", y)
    vals <- sub(".*=", "", y)
    i <- params != vals
    params <- params[i]
    vals <- vals[i]
    names(vals) <- params
    return(as.list(vals))
}







# USE FOR SLIGHTLY MODIFIED DOSE PLOT TO SHOW GSEA RESULTS
# All credit due to Guanchuang Yu et al. - see https://rdrr.io/bioc/enrichplot/src/R/utilities.R#sym-ep_str_wrap
# for ep_str_wrap and default_labler,
#' @export
ep_str_wrap <- function(string, width) {
    x <- gregexpr(' ', string)
    vapply(seq_along(x),
           FUN = function(i) {
               y <- x[[i]]
               n <- nchar(string[i])
               len <- (c(y,n) - c(0, y)) ## length + 1
               idx <- len > width
               j <- which(!idx)
               if (length(j) && max(j) == length(len)) {
                   j <- j[-length(j)]
               }
               if (length(j)) {
                   idx[j] <- len[j] + len[j+1] > width
               }
               idx <- idx[-length(idx)] ## length - 1
               start <- c(1, y[idx] + 1)
               end <- c(y[idx] - 1, n)
               words <- substring(string[i], start, end)
               paste0(words, collapse="\n")
           },
           FUN.VALUE = character(1)
    )
}

#' @export
default_labeller <- function(n) {
    function(str){
        str <- gsub("_", " ", str)
        ep_str_wrap(str, n)
    }
}

# Modified version of Enrichplot dotplot. See https://github.com/YuLab-SMU/enrichplot/blob/master/R/dotplot.R for original.
# All credit to Guangchuang Yu et al.

#' @title Modified dotplot for `gseaResult` objects
#' 
#' @description Dotplot modified from Yu et al.'s `clusterProfiler' to made more amenable to GSEA results.
#' Instead of plotting `GeneRatio`, the `NES` is plotted as a horizontal bar plot and colored by `p-adjusted`.
#' See https://github.com/YuLab-SMU/enrichplot/blob/master/R/dotplot.R for the original code. All credit
#' to Yu et al.
#' 
#' @export
dotplot_enrichResult_Col <- function(object, x = "NES", color = "p.adjust",
                             showCategory=10, size=NULL, split = NULL,
                             font.size=12, title = "", orderBy="x",
                             label_format = 30, decreasing=TRUE) {

    colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))

    if (inherits(object, c("enrichResultList", "gseaResultList"))) {
        ldf <- lapply(object, fortify, showCategory=showCategory, split=split)
        df <- dplyr::bind_rows(ldf, .id="category")
        df$category <- factor(df$category, levels=names(object))
    } else {
        df <- fortify(object, showCategory = showCategory, split=split)
        ## already parsed in fortify
        ## df$GeneRatio <- parse_ratio(df$GeneRatio)
    }

    if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
        message('wrong orderBy parameter; set to default `orderBy = "x"`')
        orderBy <- "x"
    }

    if (orderBy == "x") {
        df <- dplyr::mutate(df, x = eval(parse(text=x)))
    }

    label_func <- default_labeller(label_format)
    if(is.function(label_format)) {
        label_func <- label_format
    }

    idx <- order(df[[orderBy]], decreasing = decreasing)

    df$Description <- factor(df$Description,
                          levels=rev(unique(df$Description[idx])))
    ggplot(df, aes_string(x=x, y="Description", size=size, fill=colorBy)) +
        geom_col() +
        scale_fill_continuous(low="red", high="blue", name = color,
            guide=guide_colorbar(reverse=TRUE)) +
        scale_y_discrete(labels = label_func) +
        ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
        scale_size(range=c(3, 8)) +
        guides(size  = guide_legend(order = 1),
               color = guide_colorbar(order = 2))
}

# Modified version of the dotplot for compareClusterResult objects
# Colors by NES and uses reverse scale for p-adjusted size
# Original at https://github.com/YuLab-SMU/enrichplot/blob/master/R/dotplot.R. All credit to Yu et al.

#' @title Modified dotplot for `compareClusterResult` objects.
#' 
#' @description Dotplot modified from the `clusterProfiler` package to be more amenable to viewing
#' GSEA results from the `compareCluster` package. Instead of the doplot being colored by `p-adjust` and
#' sized by `GeneRatio`, the plot is colored by `NES` and sized by `p-adjust` for visibility in which
#' processes are enriched vs. depleted in each phenotype. See https://github.com/YuLab-SMU/enrichplot/blob/master/R/dotplot.R
#' for the original code. All credit to Yu et al.
#' 
#' @export
dotplot_compareClusterResult_gsea <- function(object, x= "Cluster", colorBy="NES",
                                         showCategory=5, by="geneRatio", size="p.adjust",
                                         split=NULL, includeAll=TRUE,
                                         font.size=12, title="", label_format = 30,
                                         group = FALSE, shape = FALSE) {
    color <- NULL
    if (is.null(size)) size <- by ## by may deprecated in future release

    df <- fortify(object, showCategory=showCategory, by=size,
                  includeAll=includeAll, split=split)

    if (by != "geneRatio")
        df$GeneRatio <- parse_ratio(df$GeneRatio)
    label_func <- default_labeller(label_format)
    if(is.function(label_format)) {
        label_func <- label_format
    }
    if (size %in% c("rowPercentage", "count", "geneRatio")) {
        by2 <- switch(size, rowPercentage = "Percentage",
                            count         = "Count",
                            geneRatio     = "GeneRatio")
    } else {
        by2 <- size
    }

    p <- ggplot(df, aes_string(x = x, y = "Description", size = by2)) +
        scale_y_discrete(labels = label_func)

    ## show multiply GO enrichment result in separate panels
    #if ("ONTOLOGY" %in% colnames(df) && length(unique(df$ONTOLOGY)) > 1){
    #    p = p + facet_grid(
    #        ONTOLOGY ~ .,
    #        scales = "free",
    #        space = "free"
    #    )
    #}

    if (group) {
        p <- p + geom_line(aes_string(color = "Cluster", group = "Cluster"), size=.3) +
          ggnewscale::new_scale_colour()
    }

    if (shape) {
        ggstar <- "ggstar"
        require(ggstar, character.only=TRUE)
        # p <- p + ggsymbol::geom_symbol(aes_string(symbol = "Cluster", fill = colorBy)) +
        p <- p + ggstar::geom_star(aes_string(starshape="Cluster", fill=colorBy)) +
            scale_fill_continuous(high="red", low="blue", guide=guide_colorbar(reverse=TRUE))
    }  else {
        p <- p +  geom_point(aes_string(color = colorBy))
    }

    p + scale_color_continuous(low="blue", high="red",
                    guide=guide_colorbar(reverse=TRUE)) +
        ylab(NULL) + ggtitle(title) + DOSE::theme_dose(font.size) +
        scale_size_continuous(range=c(8, 3)) +
        guides(size  = guide_legend(order = 1),
                color = guide_colorbar(order = 2))
}
