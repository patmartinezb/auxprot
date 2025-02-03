#' Plots volcano plot
#'
#' @description `volcano_prot()` creates a volcano plot using the output
#'   dataframe from `limma_de()`.
#'
#'
#' @param df Input dataframe. Output dataframe from `limma_de()`.
#' @param var Input character vector of lenght 1. Either `P.Value` or
#'   `adj.P.Val`.
#' @param anno Input character vector of lenght equal to the number of rows of
#'   `df`. Determines the type of labels (uniprot id, symbol...) marking the
#'   points of the volcano plot.
#'
#' @return A ggplot2 object.
#' @export volcano_prot
volcano_prot <- function(df, var, anno){
  
  if (!var %in% c("P.Value", "adj.P.Val")){
    stop("`var` should be either `P.Value` or `adj.P.Val`.")
  }
  
  fig <- df %>% 
    dplyr::mutate(dir = dplyr::case_when((.data[[var]] < 0.05) & (logFC > 0.3) ~ "Up",
                                         (.data[[var]] < 0.05) & (logFC < -0.3) ~ "Down",
                                         .default = "Not significant"),
                  lab = dplyr::case_when((.data[[var]] < 0.05) & (abs(logFC) > 0.3) ~ .data[[anno]],
                                         .default = NA_character_)) %>% 
    ggplot2::ggplot(ggplot2::aes(x = logFC, y = -log10(.data[[var]]))) +
    ggplot2::geom_vline(xintercept = c(-0.3, 0.3), 
                        col = "gray", 
                        linetype = 'dashed') +
    ggplot2::geom_hline(yintercept = -log10(0.05), 
                        col = "gray", 
                        linetype = 'dashed') +
    ggplot2::geom_point(ggplot2::aes(color = dir)) +
    ggplot2::scale_color_manual(values = c("firebrick", "grey", "#00AFBB")) +
    ggplot2::theme_classic(base_size = 14) +
    ggrepel::geom_text_repel(ggplot2::aes(label = lab))
  
  fig <- switch(var,
                "P.Value" = fig + ggplot2::labs(color = "",
                                                y = expression(-Log[10]*"(P value)"),
                                                x = "LogFC"),
                "adj.P.Val" = fig + ggplot2::labs(color = "",
                                                  y = expression(-Log[10]*"(adjusted P value)"),
                                                  x = "LogFC"))
  
  return(fig)
  
}

#' Plots heatmap with color-coded annotation.
#'
#' @description `heatmap_prot()` plots a heatmap with extra annotation in the
#'   form of colors. Wrapper for `pheatmap()`.
#'
#'
#' @param df Input dataframe. Data matrix.
#' @param meta Input dataframe. Metadata associated to `df`.
#' @param vars Input character vector of length >= 1. Indicates the variables
#'   that are going to be
#' @param z.score Input boolean vector. Default is FALSE. Performs z score.
#' @param cutrow Input numeric vector. Indicates the k groups or clusters that
#'   the row dendrogram should be divided into.
#' @param cutcol Input numeric vector. Indicates the k groups or clusters that
#'   the column dendrogram should be divided into.
#'
#' @return A pheatmap object.
#' @export heatmap_prot
heatmap_prot <- function(df, 
                         meta, 
                         vars, 
                         z.score = FALSE, 
                         # pal,
                         cutrow = NA, 
                         cutcol = NA){
  
  
  df <- as.data.frame(df)
  meta <- as.data.frame(meta)
  
  iterations = nrow(meta)
  variables = length(vars)
  
  my_sample_col <- matrix(ncol = variables, 
                          nrow = iterations)
  
  for(i in 1:variables){
    my_sample_col[,i] <- meta[[vars[i]]]
  }
  
  colnames(my_sample_col) <- vars
  my_sample_col <- as.data.frame(my_sample_col)
  
  rownames(my_sample_col) <- meta$key
  
  
  
  df_heat <- dplyr::select_if(df, is.numeric)
  
  # colnames(df_heat) <- meta$key
  rownames(df_heat) <- df$Protein.IDs
  
  if (z.score == TRUE){
    
    df_heat <- t(apply(df_heat, 1, cal_z_score))
    
  }
  
  # TODO: generalize pal param
  # if (variables == 1){
  #   pal = list(vars = pal)
  # } else {
  #   
  # }
  
  
  p <- pheatmap::pheatmap(na.omit(df_heat),
                          annotation_col = my_sample_col,
                          cutree_rows = cutrow,
                          cutree_cols = cutcol,
                          # annotation_colors = pal,
                          show_rownames = FALSE)
  
  return(p)
  
}

#' Computes PCA
#'
#' @description `pca_prot()` calculates principal component analysis (PCA) using
#'   stats::prcomp(). This function removes any missing values and then prepares
#'   the input dataframe so that prcomp() can use it.
#'
#'
#' @param df input dataframe
#' @param meta dataframe with the metadata
#' @param scale logical. Whether scale the data or not (default FALSE)
#'
#' @returns a list containing the results from computing the PCA, as well as any
#'   confounder variables and the transformed dataframe.
#'
#' @export pca_prot
pca_prot <- function(df, meta, scale = FALSE){

  df <- na.omit(df)
  
  meta <- na.omit(meta)
  
  # Selects cofounder names
  vars <- meta %>%
    dplyr::select(-BioReplicate,
                  -key) %>%
    colnames(.)
  
  out_vars <- colnames(meta)
  
  out_vars <- dplyr::setdiff(out_vars, vars)
  
  tt <- as.data.frame(t(dplyr::select(df, tidyselect::any_of(meta$key))))
  
  tt <- tt %>%
    dplyr::mutate(key = rownames(.)) %>%
    dplyr::inner_join(meta, by = "key") %>%
    dplyr::select(-all_of(out_vars)) %>%
    dplyr::select(all_of(vars), everything())
  
  rownames(tt) <- as.vector(na.omit(meta$BioReplicate[match(colnames(df), meta$key)]))
  
  res.pca <- stats::prcomp(dplyr::select(tt, -all_of(vars)), scale = scale)
  
  # List
  res.pca <- list(pca = res.pca,
                  cofounders = vars,
                  tt = tt)
  
  return(res.pca)
  
}

#' Plots PCA
#'
#' @param res.pca list as output from pca_prot()
#' @param var_p string. Variable from metadata to color.
#'
#' @returns a ggplot2 object (list)
#' 
#' @export plot_pca
plot_pca <- function(res.pca, var_p){
  
  pca <- res.pca$pca
  tt <- res.pca$tt
  
  fig <- factoextra::fviz_pca_ind(pca,
                                  geom.ind = "point",
                                  habillage = tt[[var_p]],
                                  legend.title = "", 
                                  title = "",
                                  invisible = "quali",
                                  pointsize = 4,
                                  pointshape = 19,
                                  addEllipses = TRUE,
                                  # ellipse.type = "euclid",
                                  repel = TRUE,
                                  col.ind = "black"
  ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggsci::scale_color_npg() +
    ggrepel::geom_text_repel(aes(label = rownames(tt)),
                             min.segment.length = 0, 
                             seed = 42, 
                             box.padding = 0.5)
  
  return(fig)
  
}