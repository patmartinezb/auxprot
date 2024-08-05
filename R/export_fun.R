#' Export limma result(s)
#'
#' @description `export_results()` is a function that prepares the results from
#' differential expression - performed with `{limma}` - to be exported to excel.
#' It creates a single dataframe with annotation, log2-expression
#'
#'
#' @param limma_res Input list. Output from the function `limma_de()`.
#' @param df_pre_limma Input dataframe. Data that is fed to `limma_de()`, should
#'   be normalized and, optionally, imputed.
#' @param metadata Input dataframe. Metadata associated to the data.
#' @param phospho Input vector of length 1. Either 'yes' or 'no' (default).
#' @param org Input vector of length 1. Either 'human' or 'mouse'. Default NULL.
#'
#' @return A dataframe containing annotation, log2-abundances of each sample,
#'   logFC, p value, and adjusted p value for each comparison
#' @export export_results
export_results <- function(limma_res, 
                           df_pre_limma, 
                           metadata, 
                           phospho = "no", 
                           org = NULL){
  
  # choose organism
  org <- switch(org,
                "human" = "org.Hs.eg.db",
                "mouse" = "org.Mm.eg.db")
  
  # start empty list
  listy <- list()
  
  # loop through the result from limma to extrat the necessary info
  for (i in 1:length(limma_res)){
    
    nam <- names(limma_res)[i]
    
    temp <- limma_res[[i]] %>% 
      dplyr::select(Protein.IDs, # annotation
                    logFC, # fold change
                    P.Value,
                    adj.P.Val) %>% 
      dplyr::rename_with(.,
                         ~ paste0(nam, "_", # unique names for merge
                                  .x, recycle0 = TRUE)) %>% 
      dplyr::rename_with(.,
                         ~ paste0("Protein.IDs"),
                         ends_with("Protein.IDs"))
    
    listy[[i]] <- temp
    names(listy)[i] <- names(limma_res)[i]
    
  }
  
  if (phospho == "yes"){ # alternate processing for phospho annotation
    
    dd <- listy %>% # make one dataframe from list of dataframes
      purrr::reduce(dplyr::left_join, by = c("Protein.IDs")) %>% 
      dplyr::mutate(Protein = stringr::str_extract(Protein.IDs, "[^_]+")) %>% 
      dplyr::relocate(Protein, .after = Protein.IDs)
    
    ids_dd <- clusterProfiler::bitr(unique(dd$Protein), # extra annotation
                                    fromType = "UNIPROT", 
                                    toType = c("ENTREZID",
                                             "SYMBOL",
                                             "GENENAME"), 
                                    OrgDb = org) %>% 
      dplyr::distinct(UNIPROT, .keep_all = TRUE) %>% 
      dplyr::rename(Protein = UNIPROT)
    
    dd2 <- ids_dd %>% 
      dplyr::right_join(dd, by = "Protein") %>% 
      dplyr::relocate(Protein.IDs)
    
    
    res_exp <- prot_log2(df_pre_limma, metadata) %>% # expression in log2-abundance
      dplyr::right_join(dd2, by = "Protein.IDs") %>% 
      dplyr::relocate(dplyr::any_of(c("Protein.IDs",
                                      "Protein",
                                      "ENTREZID",
                                      "SYMBOL",
                                      "GENENAME")))
    
  } else {
    
    dd <- listy %>% # make one dataframe from list of dataframes
      purrr::reduce(dplyr::left_join, by = c("Protein.IDs"))
    
    ids_dd <- clusterProfiler::bitr(unique(dd$Protein.IDs), # extra annotation
                                    fromType = "UNIPROT", 
                                    toType = c("ENTREZID",
                                               "SYMBOL",
                                               "GENENAME"), 
                                    OrgDb = org) %>% 
      dplyr::distinct(UNIPROT, .keep_all = TRUE) %>% 
      dplyr::rename(Protein.IDs = UNIPROT)
    
    dd2 <- ids_dd %>% 
      dplyr::right_join(dd, by = "Protein.IDs")
    
    
    res_exp <- prot_log2(df_pre_limma, metadata) %>% # expression in log2-abundance
      dplyr::right_join(dd2, by = "Protein.IDs") %>% 
      dplyr::relocate(dplyr::any_of(c("Protein.IDs",
                                      "ENTREZID",
                                      "SYMBOL",
                                      "GENENAME")))
    
  }
  
  return(res_exp) # for excel export
  
}