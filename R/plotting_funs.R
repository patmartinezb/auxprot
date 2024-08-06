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
