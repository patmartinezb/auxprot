#'Proteostasis Network (PN) annotation database enrichment analysis
#'
#'@description `pn_enrich()` performs an enrichment analysis using the
#'  hypergeometric test, using the PN annotation database, for both the branch
#'  and class levels.
#'
#'@param pn_db PN annotation database as a dataframe.
#'@param bio_list List populated by vectors containing the gene symbols of
#'  interest
#'@param df Original dataframe from where `bio_list` comes from. It serves as
#'  the background for creating the enrichment list.
#'@param org Organism, either "human" or "mouse".
#'@param method_p_adjust P-value adjustment method to filter by the PN terms
#'  (0.05 cutoff), either "fdr" (default) or "bonferroni"
#'
#'@returns A list with the enrichment result tables and dotplots.
#'@export pn_enrich
pn_enrich <- function(pn_db, bio_list, df, org, method_p_adjust = "fdr"){
  
  enrich_lists <- create_pn_enrich_lists(pn_db, df, org)
  
  hyp_objs <- lapply(enrich_lists, create_hyp_obj, bio_list = bio_list, df = df)
  
  enrich_plot_dfs <- sapply(names(hyp_objs), function(x){
    foo <- hyp_objs[[x]]
    
    prep_hyp_df(foo, x, method_p_adjust) 
  })
  
  enrich_plots <- lapply(enrich_plot_dfs, enrich_pn_plot, method_p_adjust)
  
  foo_list <- list("tables" = enrich_plot_dfs,
                   "figures" = enrich_plots)
  
  return(foo_list)
  
}


create_pn_enrich_lists <- function(pn_db, df, org){
  
  # Built PN lists for enrichment
  
  pn_hyper <- switch(org, # Filter by those detected in the experiment
                     "human" = pn_db %>% 
                       dplyr::filter(gene_symbol %in% df$symbol),
                     "mouse" = pn_db %>% # Use orthologs if mouse
                       dplyr::filter(ortholog_name %in% df$symbol) %>% 
                       dplyr::select(- gene_symbol) %>% 
                       dplyr::rename(gene_symbol = ortholog_name)
  )
  
  pn_hyper <- pn_hyper %>%
    dplyr::select(gene_symbol,
                  branch,
                  class) %>% 
    dplyr::mutate(class = paste(branch, class, sep = "_")) # So that the parent branch is known
  
  
  ## Branch list
  pn_branch <- unique(pn_hyper$branch)
  
  
  foo.b <- vector("list", length(pn_branch))
  names(foo.b) <- pn_branch
  
  for (i in 1:nrow(pn_hyper)){
    
    for (j in 1:length(foo.b)){
      
      if (pn_hyper$branch[i] == names(foo.b)[j]){
        
        foo.b[[j]] <- c(foo.b[[j]], pn_hyper$gene_symbol[i])
        
      }
      
    }
  }
  
  
  ## Class list
  pn_class <- unique(pn_hyper$class)
  
  foo.c <- vector("list", length(pn_class))
  names(foo.c) <- pn_class
  
  for (i in 1:nrow(pn_hyper)){
    
    for (j in 1:length(foo.c)){
      
      if (pn_hyper$class[i] == names(foo.c)[j]){
        
        foo.c[[j]] <- c(foo.c[[j]], pn_hyper$gene_symbol[i])
        
      }
      
    }
  }
  
  pn_db_list <- list("branch" = foo.b,
                     "class" = foo.c)
  
  return(pn_db_list)
  
}


create_hyp_obj <- function(bio_list, enrich_list, df){
  
  hyp_obj <- list()
  
  for (i in 1:length(bio_list)){
    
    temp <- hypeR::hypeR(bio_list[[i]], # Run hypergeometric test
                         enrich_list,
                         background = nrow(df))
    
    hyp_obj[[i]] <- temp$data
    names(hyp_obj)[i] <- names(bio_list)[i]
    
  }
  
  return(hyp_obj)
  
}


prep_hyp_df <- function(hyp_obj, level, method_p_adjust){
  
  enrich_plot <- list()
  for (i in 1:length(hyp_obj)){
    
    df_hyp <- hyp_obj[[i]] %>%
      dplyr::mutate(bonferroni = p.adjust(pval, method = "bonferroni")) %>%
      dplyr::filter(.data[[method_p_adjust]] < .05) %>%
      dplyr::mutate(cluster = names(hyp_obj)[i],
                    BgRatio = geneset / background,
                    GeneRatio = overlap / signature,
                    FoldEnrichment = GeneRatio / BgRatio
      )
    
    if (level == "class"){
      
      df_hyp <- df_hyp %>% 
        tidyr::separate_wider_delim(label, 
                                    delim = "_", 
                                    names = c("branch", "label")) %>%
        dplyr::mutate(branch = stringr::str_to_title(branch),
                      branch = stringr::str_squish(stringr::str_remove_all(branch, 
                                                                           "[^[:UPPER:]]")),
                      label = paste0(label, " (", branch, ")"))
      
    }
    
    enrich_plot[[i]] <- df_hyp
    names(enrich_plot)[i] <- names(hyp_obj)[i]
    
  }
  
  enrich_plot <- dplyr::bind_rows(enrich_plot)
  
  return(enrich_plot)
  
}


enrich_pn_plot <- function(hyp_df, method_p_adjust){
  
  foo_fig <- hyp_df %>%
    ggplot2::ggplot(ggplot2::aes(x = as.factor(cluster), 
                                 y = label, 
                                 color = .data[[method_p_adjust]])) +
    ggplot2::geom_point(ggplot2::aes(size = FoldEnrichment)) +
    ggplot2::geom_point() +
    hrbrthemes::theme_ipsum(base_size = 18) +
    ggplot2::labs(title = "",
                  y = "",
                  x = "Cluster",
                  color = "Adjusted P value",
                  size = "FoldEnrichment") +
    ggplot2::theme(plot.title.position = "plot") +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45)) +
    ggplot2::scale_y_discrete(labels = scales::label_wrap(30))
  
  return(foo_fig)
  
}