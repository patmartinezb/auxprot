#'Proteostasis Network (PN) annotation database enrichment analysis
#'
#'@description `pn_enrich()` performs an enrichment analysis using the
#'  hypergeometric test, using the PN annotation database, for both the branch
#'  and class levels.
#'
#'@param pn_db PN annotation database as a dataframe
#'@param bio_list List populated by vectors containing the gene symbols of
#'  interest
#'@param df Original dataframe from where `bio_list` comes from. It serves as
#'  the background for creating the enrichment list
#'@param org Organism, either "human" or "mouse"
#'@param method_p_val Method to calculate de P value from the hypergeometric
#'  distribution, either using the minimum-likelihood approach ('minlike',
#'  default), or doubling approach (mid-P value, 'central')
#'@param method_p_adjust P-value adjustment method to filter by the PN terms
#'  (0.05 cutoff), either "fdr" (default) or "bonferroni"
#'
#'@returns A list with the enrichment result tables and dotplots.
#'@export pn_enrich
pn_enrich <- function(pn_db, bio_list, df, org, method_p_val = "minlike", method_p_adjust = "fdr"){
  
  enrich_lists <- create_pn_enrich_lists(pn_db, df, org)
  
  hyp_objs <- create_hyp_obj(bio_list, 
                             enrich_lists, 
                             df, 
                             method_p_val)
  
  enrich_plot_dfs <- sapply(names(hyp_objs), function(x){
    foo <- hyp_objs[[x]]
    
    prep_hyp_df(foo, x) 
  })
  
  enrich_plots <- lapply(enrich_plot_dfs, enrich_pn_plot, method_p_adjust)
  
  foo_list <- list("tables" = enrich_plot_dfs,
                   "figures" = enrich_plots)
  
  return(foo_list)
  
}


#' Transforms PN database into list
#'
#' @description `create_pn_enrich_lists()` transforms the PN database into a
#'   list of vectors with two levels, branch and class. Each vector contains the
#'   genes beloging to each category.
#'
#' @param pn_db PN annotation database as a dataframe
#' @param org Organism, either "human" or "mouse"
#'
#' @returns A list of lists, where the last level is comprised of vectors with
#'   the genes belonging to each category
#'
#' @export create_pn_enrich_lists
create_pn_enrich_lists <- function(pn_db, org){
  
  # Built PN lists for enrichment
  
  pn_hyper <- switch(org,
                     "human" = pn_db,
                     "mouse" = pn_db %>% # Use orthologs if mouse
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


#' Perform Fisher's exact test
#'
#' @description `create_hyp_obj()` performs Fisher's exact test (using either
#'   the minimum-likelihood or doubling approach) and creates result dataframe
#'
#' @param bio_list List populated by vectors containing the gene symbols of
#'   interest
#' @param enrich_list PN database as a list, results from
#'   `create_pn_enrich_lists()`
#' @param df Original dataframe from where `bio_list` comes from. It serves as
#'   the background for creating the enrichment list
#' @param method_p_val Method to calculate de P value from the hypergeometric
#'   distribution, either using the minimum-likelihood approach ('minlike',
#'   default), or doubling approach (mid-P value, 'central')
#'
#' @returns A list of lists, with dataframes for each category (branch or class)
#'   and any other extra clustering (e.g., up or down), with the following
#'   variables: `label` (specific branch or class), `pval` (P value),
#'   `signature` (number of DE genes), `geneset` (number of genes in the
#'   geneset), `overlap` (number of DE genes found in the geneset), `background`
#'   (number of genes in the experiment), `hits` (readable version of overlap).
#'
#' @export create_hyp_obj
create_hyp_obj <- function(bio_list, enrich_list, df, method_p_val){
  
  hyp_obj <- lapply(enrich_list, function(x){
    
    bio_list_flat <- purrr::list_flatten(bio_list)
    
    temp <- list()
    
    for (i in 1:length(x)){
      
      res <- bio_list_flat %>%
        purrr::map(fisher_enrich, 
                   geneset = x[[i]],
                   N = nrow(df),
                   method_p_val = method_p_val,
                   label = names(x)[i])
      
      temp[[i]] <- res
      
    }
    
    res_flat <- purrr::list_flatten(temp)
    
    res_flat_comb <- purrr::map(split(res_flat, names(res_flat)), dplyr::bind_rows)
    
  })
  
  return(hyp_obj)
  
}


#' Wrapper for Fisher's exact test function, using exact2x2
#'
#' @description `fisher_enrich()` is a wrapper function of exact2x2::exact2x2(),
#'   to perform Fisher's exact test on a 2x2 contingency matrix, using either
#'   the minimum-likelihood approach or the doubling approach to extract the P
#'   value or mid-P value, respectively, from a hypergeometric distribution.
#'
#' @param signature Character vector with DE genes (symbols)
#' @param geneset Character vector with the genes (symbol) belonging to that
#'   gene set
#' @param N Numeric vector of length 1, indicating the number of genes that
#'   compose the background/universe
#' @param method_p_val Method to calculate de P value from the hypergeometric
#'   distribution, either using the minimum-likelihood approach ('minlike',
#'   default), or doubling approach (mid-P value, 'central')
#' @param label Character vector of length 1, with the name of the gene set
#'
#' @returns A dataframe with the following variables: `label` (specific branch
#'   or class), `pval` (P value), `signature` (number of DE genes), `geneset`
#'   (number of genes in the geneset), `overlap` (number of DE genes found in
#'   the geneset), `background` (number of genes in the experiment), `hits`
#'   (readable version of overlap).
#'
#' @export fisher_enrich
fisher_enrich <- function(signature, geneset, N, method_p_val, label){
  
  geneset = unique(geneset) # Remove duplicates
  overlap.hits = intersect(signature, geneset)
  
  
  n = length(signature) # number of DE genes
  m = length(geneset) # number of genes in the gene set
  k = length(overlap.hits) # number of DE genes in the gene set
  
  # Create 2x2 contingency matrix
  con_table = data.frame(de = c(k, n-k), 
                         non_de = c(m-k, N+k-n-m))
  
  # Either Fisher's exact test (two-sided), with P value calculated using the minimum-likelihood ('minlike') approach or mid-P value using the doubling approach ('central')
  pval <- switch(method_p_val,
                 "minlike" = exact2x2::exact2x2(con_table, 
                                                tsmethod = "minlike")$p.value,
                 "central" = exact2x2::exact2x2(con_table, 
                                                tsmethod = "central",
                                                midp = TRUE)$p.value)
  
  res.df <- data.frame(label = label,
                       pval = pval,
                       signature = m,
                       geneset = n,
                       overlap = k,
                       background = N,
                       hits = paste(overlap.hits, collapse = ", "))
  
  return(res.df)
  
}


#' Modify the enrichment dataframe
#'
#' @description `prep_hyp_df()` takes the resulting dataframe from the
#'   enrichment analysis, calculates the adjusted P value, as well as the
#'   BgRatio, GeneRatio and the FoldEnrichment. It also makes the class labels
#'   readable
#'
#' @param hyp_obj list resulting from `create_hyp_obj()`
#' @param level either 'class' or 'branch'
#'
#' @returns A dataframe as the one from `fisher_enrich()`, but with the
#'   following extra variables: `fdr` (adjusted P value by BH method),
#'   `bonferroni` (adjusted P value by Bonferroni method), `cluster` (extra
#'   grouping, if any), `BgRatio` (number of all genes in specific term / number
#'   of universal genes), `GeneRatio` (number of genes enriched in specific term
#'   / number of input genes), `FoldEnrichment` (GeneRatio / BgRatio)
#'
#' @export prep_hyp_df
prep_hyp_df <- function(hyp_obj, level){
  
  enrich_plot <- list()
  for (i in 1:length(hyp_obj)){
    
    df_hyp <- hyp_obj[[i]] %>%
      dplyr::mutate(bonferroni = p.adjust(pval, method = "bonferroni"),
                    fdr = p.adjust(pval, method = "BH")) %>%
      # dplyr::filter(.data[[method_p_adjust]] < .05) %>%
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


#' Enrichment dot plot
#'
#' @description `enrich_pn_plot()` creates a dot plot showcasing the PN
#'   enrichment
#'
#' @param hyp_df List resulting from `prep_hyp_df()`
#' @param method_p_adjust Type of adjusted P value to use in the dot plot
#'
#' @returns A ggplot2 figure, a dot plot, showing the branch/class labels, the
#'   fold enrichment, adjusted P value and faceted (or not) depending on the
#'   cluster/grouping.
#'
#' @export enrich_pn_plot
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