% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pn_funs.R
\name{create_pn_enrich_lists}
\alias{create_pn_enrich_lists}
\title{Transforms PN database into list}
\usage{
create_pn_enrich_lists(pn_db, org)
}
\arguments{
\item{pn_db}{PN annotation database as a dataframe}

\item{org}{Organism, either "human" or "mouse"}
}
\value{
A list of lists, where the last level is comprised of vectors with
the genes belonging to each category
}
\description{
\code{create_pn_enrich_lists()} transforms the PN database into a
list of vectors with two levels, branch and class. Each vector contains the
genes beloging to each category.
}
