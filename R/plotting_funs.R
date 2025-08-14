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
                                  # addEllipses = TRUE,
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


#' Boxplots for proteomics data
#'
#' @description `prot_boxplot()` creates boxplots for all samples - can be
#'   reordered or not, based on the intensity/abundance variable
#'
#' @param df Dataframe to be plotted
#' @param metadata Dataframe with associated metadata
#' @param fill String with name of variable to fill the boxplots
#' @param reorder TRUE or FALSE (default). Reorder boxplots based on intensity
#' @param var_reorder String with name of the variable to use for reorder
#'
#' @returns A ggplot2 object
#' @export prot_boxplot
prot_boxplot <- function(df, metadata, fill, reorder = FALSE, var_reorder = "counts"){
  
  if (reorder == TRUE){
    
    df <- df %>%
      dplyr::select(Protein.IDs,
                    which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
      tidyr::pivot_longer(!Protein.IDs, names_to = "key", values_to = "counts") %>%
      dplyr::left_join(metadata, by = "key")
    
    df$BioReplicate <- reorder(as.factor(df$BioReplicate), df[[var_reorder]], na.rm = TRUE)
    
    plot <- ggplot2::ggplot(df) +
      ggplot2::geom_boxplot(
        ggplot2::aes(x = BioReplicate, y = counts, fill = .data[[fill]])
      ) +
      ggplot2::labs(y = "Log2-intensity", x = "Samples") +
      ggsci::scale_fill_npg() +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, vjust = 0.5, size = 12))
    
    return(plot)
    
  } else {
    
    plot <- df %>%
      dplyr::select(Protein.IDs, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
      tidyr::pivot_longer(!Protein.IDs, names_to = "key", values_to = "counts") %>%
      dplyr::left_join(metadata, by = "key") %>%
      ggplot2::ggplot() +
      ggplot2::geom_boxplot(
        ggplot2::aes(x = BioReplicate, y = counts, fill = .data[[fill]])
      ) +
      ggplot2::labs(y = "Log2-intensity", x = "Samples") +
      ggsci::scale_fill_npg() +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, vjust = 0.5, size = 12))
    
    return(plot)
    
  }
  
}


#' Density plot for proteomics data
#'
#' @description `prot_denplot()` creates density plots of proteomics data
#'
#' @param df Dataframe to be plotted
#' @param metadata Dataframe with associated metadata
#' @param color String with variable name to be used as grouping variable
#'
#' @returns A ggplot2 object
#' @export prot_denplot
prot_denplot <- function(df, metadata, color = "BioReplicate"){
  
  plot <- df %>%
    dplyr::select(Protein.IDs, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    tidyr::pivot_longer(!Protein.IDs, names_to = "key", values_to = "counts") %>%
    dplyr::left_join(metadata, by = "key") %>%
    ggplot2::ggplot() +
    ggplot2::geom_density(
      ggplot2::aes(x = counts, color = .data[[color]])
    ) +
    ggplot2::labs(y = "Density", x = "Samples") +
    ggsci::scale_fill_npg() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, vjust = 0.5, size = 12))
  
  return(plot)
  
}


#' Principal component analysis (or PCA) regression plot
#'
#' @description `cofounders_plot()` uses the information provided by calculating
#'   a PCA, extracts the top principal components (PCs), which are orthogonal
#'   variables explaining the information contained within it. Then, the
#'   different clinical or technical variables indicated in the metadata are
#'   regressed against the PCs. The heatmap shows the values of the R² and
#'   significance (FDR), indicating the strength of the relationship between the
#'   variables and the PCs.
#'
#' @param res.pca Output list from `pca_prot()`
#'
#' @returns A ggplot2 object
#' @export cofounders_plot
cofounders_plot <- function(res.pca){

  # Define variables
  pca <- res.pca$pca
  cofounders2 <- res.pca$cofounders
  tt <- res.pca$tt
  
  cofounders <- c() # Check that the confounders have more than one level
  for (i in 1:length(cofounders2)){
    
    if (!is.numeric(tt[,i])){
      
      levs <- nlevels(as.factor(tt[, i]))
      
      if (levs > 1){
        
        cofounders[i] <- cofounders2[i]
        cofounders <- cofounders[!is.na(cofounders)]
      }
    } else{
      
      cofounders[i] <- cofounders2[i]
      
    }
  }
  
  # Number of PCs to interrogate
  if (dim(pca$x)[2] > 15){
    ndims <- 15
  } else {
    ndims <- dim(pca$x)[2]
  }
  
  # We may have both numerical and categorical variables, so we will do linear modeling:
  confound <- calc_confound(pca$x, tt, cofounders, ndims)
  
  # Get p-values and adjust
  ps <- sapply(confound, function(x){sapply(x, function(y){y$p})})
  ps.adj <- matrix(p.adjust(c(as.matrix(ps)), method = 'fdr'), ncol = ncol(ps), byrow = F)
  
  # Get R2 values and filter by FDR
  r2s <- sapply(confound, function(x){sapply(x, function(y){y$R2_adj})})
  r2s <- matrix(as.character(round(r2s, 2)),ncol = ncol(ps.adj))
  
  r2s[ps.adj >= 0.05] <- ""
  
  rownames(r2s) <- rownames(ps)
  colnames(r2s) <- colnames(ps)
  
  # Also set as NA non-significant ps.adj
  ps.adj[ps.adj >= 0.05] <- NA
  ps.adj.log <- -log10(ps.adj)
  
  rownames(ps.adj.log) <- rownames(ps)
  colnames(ps.adj.log) <- colnames(ps)
  
  
  # Transform for plotting
  ps.adj.log <- ps.adj.log %>%
    tibble::as_tibble() %>%
    tibble::rowid_to_column(var = "PC") %>%
    dplyr::mutate(PC = paste0("PC", PC)) %>%
    tidyr::pivot_longer(!PC, names_to = "vars", values_to = "fdr")
  
  
  # Transform for plotting
  r2s <- r2s %>%
    tibble::as_tibble() %>%
    tibble::rowid_to_column(var = "PC") %>%
    dplyr::mutate(PC = paste0("PC", PC)) %>%
    tidyr::pivot_longer(!PC, names_to = "vars", values_to = "r2s")
  
  
  # Merge both FDR and R2 to make plot
  ff <- r2s %>%
    dplyr::left_join(ps.adj.log, by = c("PC", "vars")) %>%
    dplyr::mutate(PC = forcats::fct_reorder(PC, readr::parse_number(PC), .desc = TRUE))
  
  # Plot
  
  fig <- ggplot2::ggplot(ff, ggplot2::aes(x = vars, y = PC, fill = fdr)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = r2s)) +
    ggplot2::scale_fill_continuous(low = "thistle2",
                                   high = "darkred",
                                   guide = "colorbar",
                                   na.value="white") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90)) +
    ggplot2::labs(x = "",
                  y = "",
                  fill = "-log10 (FDR)",
                  caption = "Variables and PCs (FDR < 0.05)")
  
  return(fig)
  
}


#' Barplot for proteomics data
#'
#' @description `prot_barplot()` creates barplots of proteomics data
#'
#' @param df Dataframe to be plotted
#' @param metadata Dataframe with associated metadata
#' @param fill String with variable name to be used as filling/grouping variable
#'
#' @returns A ggplot2 object
#' @export prot_barplot
prot_barplot <- function(df, metadata, fill = "Condition"){
  
  reporter_names_clean <- gsub(" ", ".", metadata$key)
  
  df_plot <- df %>%
    dplyr::select(
      Protein.IDs,
      dplyr::any_of(reporter_names_clean)
    ) %>%
    tidyr::pivot_longer(!Protein.IDs, names_to = "key", values_to = "vals") %>%
    dplyr::left_join(metadata, by = "key") %>%
    dplyr::select(
      key,
      vals,
      Condition,
      BioReplicate,
      Mixture,
      tidyselect::all_of(fill)
    ) %>%
    dplyr::filter(!is.na(vals)) %>%
    dplyr::group_by(BioReplicate, .data[[fill]]) %>%
    dplyr::count() %>%
    ggplot2::ggplot() +
    ggplot2::geom_bar(ggplot2::aes(reorder(BioReplicate, n), n, fill = .data[[fill]]),
                      stat = "identity", color = "black") +
    ggplot2::geom_text(ggplot2::aes(reorder(BioReplicate, n), n, label = n),
                       hjust = 1.1, angle = 90) +
    ggplot2::xlab("Sample") +
    ggplot2::ylab("Count") +
    ggplot2::ggtitle("Number of proteins per sample") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
    ggsci::scale_fill_npg()
  
  return(df_plot)
  
}