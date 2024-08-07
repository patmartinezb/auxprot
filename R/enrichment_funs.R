#' Wrapper function for GO enrichment
#'
#' @description `enrich_prot()` is a wrapper function for both
#' `{clusterProfile}`'s `enrichGO()` and `gseGO()` functions. It performs
#' enrichment analysis (either GO or GSEA) of all ontologies (CC, BP and MF).
#'
#'
#' @param de Input character vector. If GO enrichment, a vector of proteins that
#'   are going to be fed to the enrichment analysis, and UNIPROT ids are
#'   mandatory. If GSEA, a named numeric, ranked protein vector: fold change
#'   sorted in decreasing order, and names with ENTREZID.
#' @param universe Input character vector with the background proteins, only for
#'   GO analysis. UNIPROT ids are mandatory.
#' @param org Character vector of length 1. Either 'human' or 'mouse'.
#' @param type.enrich Character vector of length 1. Either 'go' (default) or
#'   'gsea'.
#'
#' @return A list of lists, with the results from the enrichent.
#' @export enrich_prot
enrich_prot <- function(de, universe, org, type.enrich = "go"){
  
  org <- switch(org,
                "human" = org.Hs.eg.db::org.Hs.eg.db,
                "mouse" = org.Mm.eg.db::org.Mm.eg.db)
  
  if (type.enrich == "go"){
    
    go.cc <- clusterProfiler::enrichGO(gene = de,
                                       OrgDb         = org,
                                       universe      = universe, 
                                       keyType       = 'UNIPROT',
                                       ont           = "CC",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.01,
                                       qvalueCutoff  = 0.05,
                                       readable = TRUE)
    
    go.bp <- clusterProfiler::enrichGO(gene = de,
                                       OrgDb         = org,
                                       universe      = universe, 
                                       keyType       = 'UNIPROT',
                                       ont           = "BP",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.01,
                                       qvalueCutoff  = 0.05,
                                       readable = TRUE)
    
    go.mf <- clusterProfiler::enrichGO(gene = de,
                                       OrgDb         = org,
                                       universe      = universe, 
                                       keyType       = 'UNIPROT',
                                       ont           = "MF",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.01,
                                       qvalueCutoff  = 0.05,
                                       readable = TRUE)
    
    go.list <- list(go.cc = go.cc,
                    go.bp = go.bp,
                    go.mf = go.mf)
    
    return(go.list)
    
  } else {
    
    gsea.cc <- clusterProfiler::gseGO(geneList = de,
                                      OrgDb        = org,
                                      ont          = "CC",
                                      minGSSize    = 100,
                                      maxGSSize    = 500,
                                      pvalueCutoff = 0.05,
                                      verbose      = FALSE)
    
    
    gsea.bp <- clusterProfiler::gseGO(geneList = de,
                                      OrgDb        = org,
                                      ont          = "BP",
                                      minGSSize    = 100,
                                      maxGSSize    = 500,
                                      pvalueCutoff = 0.05,
                                      verbose      = FALSE)
    
    gsea.mf <- clusterProfiler::gseGO(geneList = de,
                                      OrgDb        = org,
                                      ont          = "MF",
                                      minGSSize    = 100,
                                      maxGSSize    = 500,
                                      pvalueCutoff = 0.05,
                                      verbose      = FALSE)
    
    gsea.list <- list(gsea.cc = clusterProfiler::setReadable(gsea.cc, org),
                      gsea.bp = clusterProfiler::setReadable(gsea.bp, org),
                      gsea.mf = clusterProfiler::setReadable(gsea.mf, org))
    
    return(gsea.list)
    
  }
}

#' Creates dotplot for GO enrichment results
#'
#' @param ego An enrichResult instance, output from `enrich_prot()` function.
#'
#' @return A ggplot2 object.
#' 
#' @export dotplot_enrich
dotplot_enrich <- function(ego){
  
  fig <- ego@result %>% 
    dplyr::filter(p.adjust < .05) %>%
    tidyr::separate_wider_delim(GeneRatio, 
                                delim = "/", 
                                names = c("var1", "var2")) %>%
    type.convert(as.is = T) %>% 
    dplyr::mutate(GeneRatio = var1 / var2) %>% 
    ggplot2::ggplot(ggplot2::aes(GeneRatio, reorder(Description, GeneRatio))) +
    ggplot2::geom_point(ggplot2::aes(color = pvalue)) +
    ggplot2::theme_classic() +
    ggplot2::labs(y = "Description")
  
  return(fig)
  
}