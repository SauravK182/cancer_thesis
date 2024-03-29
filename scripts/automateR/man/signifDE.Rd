% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dge_analysis.R
\name{signifDE}
\alias{signifDE}
\title{Isolate DGE genes}
\usage{
signifDE(results, threshold = 0.05, lfc = 0)
}
\arguments{
\item{results}{Object produced by a call to \code{results} or \code{lfcShrink} in \code{DESeq2}}

\item{threshold}{FDR-adjusted p-value threshold. Only genes below this threshold will be kept. Default: 0.05.}

\item{lfc}{Log-fold change minimum to be imposed to extract DE genes. Default = 0 (i.e., no threshold)}
}
\value{
Subset of the \code{results} object, containing only genes below the user-defined threshold.
}
\description{
Takes in a \code{results} object produced by \code{DESeq2} and isolates only differentially expressed genes.
}
