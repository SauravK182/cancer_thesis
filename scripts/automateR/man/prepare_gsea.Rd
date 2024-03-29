% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment.R
\name{prepare_gsea}
\alias{prepare_gsea}
\title{Prepare \code{DESeq} results for GSEA}
\usage{
prepare_gsea(
  deseqObject,
  rankBy = "log2FoldChange",
  colSymbol = "rownames",
  useSymbol = "ENSEMBL",
  method = "clusterProfiler",
  convertSymbol = TRUE
)
}
\arguments{
\item{deseqObject}{A \code{results} object from \code{DESeq2}, or any other valid data-frame like object that possesses
a numerical ranking score and gene identifers}

\item{rankBy}{Column name in the data frame to rank genes by. Default = "log2FoldChange".}

\item{colSymbol}{Column name for the column in the data frame containing the gene identifiers.
Note that these identifiers must be a valid keytype in the org.Hs.eg database. Default = "rownames" (i.e., gene
identifiers are the rownames of the object).}

\item{useSymbol}{Keytype for the gene identifiers housed in the column \code{colSymbol}. Must be a valid
keytype that can be queried by org.Hs.eg.db. Default = "ENSEMBL"}

\item{method}{Either \code{"cluterProfiler"} or \code{"fgsea"}, indicating which package will be used for downstream analysis.
Default = "clusterProfiler"}

\item{convertSymbol}{Boolean value indicating whether to convert symbols to the appropriate type for \code{clusterProfiler}
or \code{fgsea}. Default = TRUE}
}
\value{
Named numeric vector sorted from highest to lowest \code{rankBy} stat. Names correspond to the gene IDs given the keytype
required by \code{clusterProfiler} (Entrez ID) or \code{fgsea} (HGNC Symbol), respectively.
}
\description{
Will automatically create a ranked list of genes with proper identifiers from a \verb{DESeq2 results} object for
either \code{clusterProfiler} or \code{fgsea}. Utilizes org.Hs.eg.db for any necessary gene identifier conversions.
}
