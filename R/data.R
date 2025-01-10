#' Human Proteostasis Network Annotation
#'
#' A dataset containing a curated and annotated list of all (known) proteostasis
#' network genes. The variables are as follows:
#'
#' @format A data frame with 3774 rows and 29 variables:
#' \describe{
#'   \item{order}{order}
#'   \item{gene_symbol}{Gene symbols for each gene}
#'   \item{ortholog_name}{Gene symbols corresponding to mouse orthologs}
#'   \item{gene_name}{Gene names for each gene}
#'   \item{number_ann}{Number of branches the gene appears on}
#'   \item{in_branches}{Branches the gene appears on}
#'   \item{branch}{Top-most level of the PN annotation hierarchy}
#'   \item{class}{Second level of the PN annotation hierarchy}
#'   \item{group}{Third level of the PN annotation hierarchy}
#'   \item{type}{Fourth level of the PN annotation hierarchy}
#'   \item{subtype}{Fifth level of the PN annotation hierarchy}
#'   \item{principal_domains}{Sixth level of the PN annotation hierarchy}
#'   \item{auxiliary_domains}{Seventh and last level of the PN annotation hierarchy}
#'   \item{ensg}{Ensembl ID for each gene}
#'   \item{gene_id}{ENTREZID for each gene}
#'   \item{uni_prot_id}{UniProt ID for each gene}
#'   \item{gene_synonyms}{Gene symbol synonyms for each gene}
#'   \item{protein_name}{Protein names}
#'   \item{length_aa}{Protein length in aminoacids}
#'   \item{mass_da}{Mass of the protein in daltons}
#'   \item{notes}{Other comments}
#'   \item{references}{References}
#' }
#'
#' @source <https://www.proteostasisconsortium.com/pn-annotation/>
"pn_db"

#' PhosphoSitePlus kinase-substrate mouse database
#'
#' A dataset containing a paired list of kinases and their respective
#' substrates. The variables are as follows:
#'
#' @format A data frame with 2264 rows and 2 variables:
#' \describe{
#'   \item{kinase}{Kinase gene symbols}
#'   \item{substrate}{Aminoacid sequence (14aa) of the corresponding substrate}
#' }
#'
#' @source <https://www.phosphosite.org/homeAction.action>
"psp"