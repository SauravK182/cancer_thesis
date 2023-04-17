#' @name Overlap-class
#'
#' @slot permutedf data frame representing the result of manual permutation tests - in particular,
#' should have one column for permuted number of overlaps, and optionally, a second column that indicates
#' whether the permuted overlap was greated than the observed overlap.
#' @slot observed Observed number of overlaps of two sets.
#' @slot expected Expectation from the permutation-based distribution for number of overlaps by chance.
#' @slot pval Manually computed p-value for observed vs. expected enrichment - i.e.,
#' P(X >= k) where X is the random variable representing the permuted overlaps, and k is the observed
#' overlap
#' @slot type I don't remember - probably for upregulated or downregulated gene expression. Optional.
#'
#' @return
#' @export
#'
#' @examples
setClass("Overlap", representation = representation(
  permutedf = "data.frame",
  observed = "numeric",
  expected = "numeric",
  pval = "numeric",
  type = "character"
))



#' @name ProximalGeneExp-class
#'
#' @slot upregProx data frame representing the log fold changes of genes annotated to be near significantly
#' increased enrichment for histone marks or chromatin accessibility
#' @slot downregProx data frame representing the log fold changes of genes annotated to be near significantly
#' decreased enrichment for histone marks/chromatin accessibility.
#' @slot upregIntPercent Number of genes tested for DE that coincided with the set of all genes that are
#' annotated with a LFC of differential binding > 0
#' @slot downregIntPercent Number of genes tested for DE that coincided with the set of all genes  that
#' were annotated with a LFC of differential binding < 0
#' @slot totalIntersect Number of genes tested for DE that coincided with the set of all genes that
#' were annotated with a significant difference in ChIP signal
#' @slot totalGenes Total number of genes tested for differential expression
#' @slot binomTestUp Number representing the p-value that the median LFC of the genes in `upregProx` are > 0
#' as computed by the exact binomial test
#' @slot binomTestDown Number representing the p-value that the median LFC of the genes in `downregProx` are
#' < 0 as computed by the exact binomial test.
#'
#' @return
#' @export
#'
#' @examples
setClass("ProximalGeneExp", representation = representation(
  upregProx = "data.frame",
  downregProx = "data.frame",
  upregIntPercent = "numeric",
  downregIntPercent = "numeric",
  totalIntersect = "numeric",
  totalGenes = "numeric",
  binomTestUp = "numeric",
  binomTestDown = "numeric"
))

setGeneric(name = "intersect_prox", def = function(obj1, obj2, type = "") standardGeneric("intersect_prox"))



#' @export
setMethod(
          f = "intersect_prox",
          signature = "ProximalGeneExp",
          definition = function(obj1, obj2, type = "up") {
            if (type == "up") {
              proxgenes.one <- obj1@upregProx
              proxgenes.two <- obj2@upregProx
            } else if (type == "down") {
              proxgenes.one <- obj1@downregProx
              proxgenes.two <- obj2@downregProx
            } else {
              stop("Variable 'type' should be 'up' or 'down' to indicate which intersection to take.")
            }

            genes.in.both <- intersect(rownames(proxgenes.one), rownames(proxgenes.two))
            mutual.lfc <- data.frame(log2FoldChange = proxgenes.one[genes.in.both, ],
                                     row.names = genes.in.both)
            prox.list <- list(prox_one = proxgenes.one, prox_two = proxgenes.two)
            prox.list <- lapply(prox.list, function(df) {
              df %>%
                rownames_to_column(var = "rowname") %>%
                filter(!(rowname %in% genes.in.both)) %>%
                column_to_rownames(var = "rowname")
            })

            prox.list[["mutual"]] <- mutual.lfc
            merged.df <- purrr::reduce(prox.list, .f = function(df1, df2) {
              merge(df1, df2, by = 0, all = TRUE) %>%
                remove_rownames() %>%
                column_to_rownames(var = "Row.names")
            })

            return(merged.df)
})

setGeneric(name = "get_proximal_genes", def = function(ranges, deseqObject, minlfc = 0, maxlfc = Inf) standardGeneric("get_proximal_genes"))


#' @export
setMethod(
          f = "get_proximal_genes",
          signature = c("GRanges", "DESeqResults"),
          definition = function(ranges, deseqObject, minlfc = 0, maxlfc = Inf) {
              # Get up/down changing features
              anno.up <- ranges[ranges$Fold > minlfc & ranges$Fold < maxlfc, ]
              anno.down <- ranges[ranges$Fold < -minlfc & ranges$Fold > -maxlfc, ]
              anno.up.features <- anno.up$feature
              anno.down.features <- anno.down$feature

              # Subset DESeq2 results object to isolate these features in particular
              anno.up.subset <- intersect(anno.up.features, rownames(deseqObject))
              anno.down.subset <- intersect(anno.down.features, rownames(deseqObject))
              gene.anno.up <- deseqObject[anno.up.subset, ] %>%
                as.data.frame() %>%
                dplyr::select(log2FoldChange)
              gene.anno.down <- deseqObject[anno.down.subset, ] %>%
                as.data.frame() %>%
                dplyr::select(log2FoldChange)

              # Calculate intersection stats
              upreg.int <- length(anno.up.subset) / length(anno.up.features)
              downreg.int <- length(anno.down.subset) / length(anno.down.features)
              total.int <- length(anno.up.subset) + length(anno.down.subset)

              # Calculate sign test for median > 0 or median < 0
              binomUp <- BSDA::SIGN.test(x = gene.anno.up$log2FoldChange, md = 0, alternative = "greater")$p.value
              binomDown <- BSDA::SIGN.test(x = gene.anno.down$log2FoldChange, md = 0, alternative = "less")$p.value

              # Instantiate ProximalGeneExp class and return
              pge.obj <- new("ProximalGeneExp",
                             upregProx = gene.anno.up,
                             downregProx = gene.anno.down,
                             upregIntPercent = upreg.int,
                             downregIntPercent = downreg.int,
                             totalIntersect = total.int,
                             totalGenes = length(rownames(deseqObject)),
                             binomTestUp = binomUp,
                             binomTestDown = binomDown)
              return(pge.obj)
          }
)
