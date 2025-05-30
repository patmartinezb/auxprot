% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pn_funs.R
\name{fisher_enrich}
\alias{fisher_enrich}
\title{Wrapper for Fisher's exact test function, using exact2x2}
\usage{
fisher_enrich(signature, geneset, N, method_p_val, label)
}
\arguments{
\item{signature}{Character vector with DE genes (symbols)}

\item{geneset}{Character vector with the genes (symbol) belonging to that
gene set}

\item{N}{Numeric vector of length 1, indicating the number of genes that
compose the background/universe}

\item{method_p_val}{Method to calculate de P value from the hypergeometric
distribution, either using the minimum-likelihood approach ('minlike',
default), or doubling approach (mid-P value, 'central')}

\item{label}{Character vector of length 1, with the name of the gene set}
}
\value{
A dataframe with the following variables: \code{label} (specific branch
or class), \code{pval} (P value), \code{signature} (number of DE genes), \code{geneset}
(number of genes in the geneset), \code{overlap} (number of DE genes found in
the geneset), \code{background} (number of genes in the experiment), \code{hits}
(readable version of overlap).
}
\description{
\code{fisher_enrich()} is a wrapper function of exact2x2::exact2x2(),
to perform Fisher's exact test on a 2x2 contingency matrix, using either
the minimum-likelihood approach or the doubling approach to extract the P
value or mid-P value, respectively, from a hypergeometric distribution.
}
