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
