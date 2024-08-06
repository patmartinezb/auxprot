#' Performs hierarchical clustering and plots number of clusters
#'
#' @description `dend_calc()` performs hierarchical clustering and prints a
#' silhouette plot to help choose number of clusters.
#'
#' @param df Input dataframe. Raw data matrix, with protein ids and sample
#'   abundances.
#' @param meta Input dataframe. Metadata associated to `df`.
#' @param plot Boolean indicating if silhouette plot should be displayed or not.
#'   Default is TRUE.
#'
#' @return An object of class "hclust".
#'
#' @export dend_calc
dend_calc <- function(df, meta, plot = TRUE){
  
  df.dend <- dplyr::select_if(prot_log2(df, meta), is.numeric)
  rownames(df.dend) <- df$Protein.IDs
  
  df.dend <- t(apply(df.dend, 1, cal_z_score))
  
  if (plot == TRUE){
    print(factoextra::fviz_nbclust(df.dend,
                                   factoextra::hcut, 
                                   method = 'silhouette'))
  }
  
  hc <- hclust(dist(df.dend), method = "complete")
  
  return(hc)
  
}

#' Assignment of phosphopeptides to respective clusters.
#'
#' @description `select_cluster()` assigns each phosphopeptide to its cluster,
#' based on the number of clusters selected by the user.
#'
#'
#' @param hc An object of class "hclust" created with the `dend_calc()`
#'   function.
#' @param k Number of clusters.
#' @param df Input dataframe. Raw data matrix, with protein ids and sample
#'   abundances.
#'
#' @return A list of dataframes with as many levels as specified in `k`.
#'
#' @export select_cluster
select_cluster <- function(hc, k = NULL, df){
  
  if (is.null(k) | (k <= 0)){
    stop("k should be a valid, positive number, and different from zero.")
  }
  
  hr <- cutree(hc, k)
  
  cl.list <- list()
  
  for (i in 1:k){
    
    cluster.n <- names(hr)[hr == i]
    cl <- df %>% 
      dplyr::filter(Protein.IDs %in% cluster.n) %>% 
      tibble::as_tibble()
    
    cl.list[[i]] <- cl
    names(cl.list)[i] <- paste0("cluster", i)
    
  }
  
  return(cl.list)
  
}

#' Plots cluster profile.
#' 
#' @description
#' `cluster_plot()` plots a profile of the specified cluster of phosphopeptides. Values a z-score transformed previous to plotting.
#'
#' @param df Input dataframe. Raw data matrix, with protein ids and sample
#'   abundances.
#' @param meta Input dataframe. Metadata associated to `df`.
#' @param hc An object of class "hclust" created with the `dend_calc()`
#'   function. 
#' @param k Number of clusters.
#'
#' @return A ggplot2 object.
#' 
#' @export cluster_plot
cluster_plot <- function(df, meta, hc, k){
  
  if ((length(k) > 1) | (length(k) < 1)){
    stop("`k` should be a numeric vector of length 1.")
  }
  
  if (is.null(k) | (k <= 0)){
    stop("k should be a valid, positive number, and different from zero.")
  }
  
  hr <- cutree(hc, k)
  
  plot.ca <- t(apply(dplyr::select_if(prot_log2(df, meta), 
                                      is.numeric), 
                     1, 
                     cal_z_score)) %>%
    as.data.frame() %>% 
    dplyr::mutate(cluster = hr) # hr is in the same order as df
  
  rownames(plot.ca) <- df$Protein.IDs
  
  plot.ca$cluster <- as.factor(plot.ca$cluster)
  
  plot.ca.ord <- colnames(plot.ca)[order(sub("_.*", "", colnames(plot.ca)))]
  names.ord <- meta$key[match(plot.ca.ord, meta$key)]
  
  sum_meta <- meta %>% 
    dplyr::select(Condition,
                  key) %>% 
    dplyr::arrange(factor(key, levels = names.ord)) %>%
    dplyr::mutate(Condition = forcats::fct_inorder(Condition)) %>% 
    dplyr::group_by(Condition) %>% 
    dplyr::count()
  
  n_lines <- dplyr::pull(sum_meta, n)
  
  # calculate/generalize where the dashed lines should go
  vlines <- tibble::tibble(lines = cumsum(n_lines[1:length(n_lines)-1]) + 0.5)
  
  g.col = ncol(plot.ca) # column where with the cluster information
  
  fig.list <- list()
  
  for (i in 1:k){
    
    fig <- plot.ca %>%
      dplyr::filter(cluster == i) %>%
      GGally::ggparcoord(columns = 1:g.col-1, 
                         groupColumn = g.col, 
                         order = match(plot.ca.ord, colnames(plot.ca)[-g.col]),
                         scale = "globalminmax",
                         showPoints = FALSE, 
                         title = paste("Cluster", i),
                         alphaLines = 0.3
      ) + 
      ggplot2::scale_color_manual(values = c("#a29eb8")) +
      hrbrthemes::theme_ipsum(grid = "",
                              axis_title_size = 14,
                              base_size = 12) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 16),
        axis.text.x = ggplot2::element_text(angle = 90, 
                                            vjust = 0.5, 
                                            hjust=1, 
                                            size = 10),
        axis.line = ggplot2::element_line(colour = "grey50")
      ) +
      ggplot2::scale_x_discrete(breaks = plot.ca.ord,
                                labels = names.ord) +
      ggplot2::xlab("Samples") +
      ggplot2::ylab("Standarized Log2-intensities") +
      ggplot2::stat_summary(ggplot2::aes(group = cluster), 
                            fun = mean, 
                            geom = "line", 
                            colour = "#8a3223", 
                            size = 1.1) +
      ggplot2::geom_vline(data = vlines, ggplot2::aes(xintercept = lines), 
                          linetype = "dashed", 
                          color = "grey50") +
      ggplot2::annotate("text", 
                        x = cumsum(n_lines) - n_lines/2 + 0.5,
                        y = 2.8, 
                        label = dplyr::pull(sum_meta, Condition), 
                        size = 3)
    
    fig.list[[i]] <- fig
    names(fig.list)[i] <- paste0("cluster", i)
    
  }
  
  return(fig.list)
  
}