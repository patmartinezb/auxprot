#' WGCNA wrapper function
#'
#' @description `wgcna_prot()` is a wrapper function of several functions from
#' the well-known `{WGCNA}` package, that are run in a pipeline to conduct the
#' WGCNA analysis, start to end. It outputs the figures and list with clusters.
#'
#'
#' @param df Input dataframe. Data matrix with normalized (log2) abundances.
#' @param meta Input dataframe. Metadata associated to `df`.
#' @param Group Input dataframe. Dataframe with the grouping, should have two
#'   columns: Sample and Group.
#' @param traits Input dataframe. For module-trait correlation, should have
#'   Sample column, plus the rest of clinical variables, dummy-encoded (0-1),
#'   and controls are not included.
#' @param var Input character vector of length 1, optional. Extra variable to
#'   take into account. Default is NULL.
#' @param org Input character vector of length 1. Either "human" or "mouse".
#' @param RCutoff Cut-off value for the correlation R square, used to adjust the
#'   soft threshold selection, default is 0.85.
#' @param MCutheight Height in the dendrogram, establishes if modules are going
#'   to be merged, default is 0.1.
#' @param PowerUpper Maximum soft-thresholding power to choose from, default is
#'   30.
#' @param minModuleSize Controls the minimum number of proteins in each cluster
#'   - if a few hundred/a thousand proteins, should be between 20-30, more than
#'   that, 50. Default is 20.
#' @param verbose Print intermediate steps. Default is TRUE.
#'
#' @return A list of lists, with figures and the proteins that conform each
#'   cluster.
#' @export wgcna_prot
wgcna_prot <- function(df, 
                       meta, 
                       Group, 
                       traits, 
                       var = NULL,
                       org,
                       RCutoff = 0.85, 
                       MCutheight = 0.1, 
                       PowerUpper = 30, 
                       minModuleSize = 20,
                       verbose = TRUE){
  
  org <- switch(org,
                "human" = "org.Hs.eg.db",
                "mouse" = "org.Mm.eg.db")
  
  if ("Protein.IDs" %in% colnames(df)){
    allDat <- df %>% 
      as.data.frame()
  } else if (!"Protein.IDs" %in% colnames(df)){
    allDat <- df %>%
      tibble::rownames_to_column(var = "Protein.IDs")%>% 
      as.data.frame()
  } else {
    stop("There is no Protein.IDs column in the data")
  }
  
  rownames(allDat) <- allDat$Protein.IDs
  
  # Only abundance data - samples as rows and proteins as columns
  IgnoreCols = 1 	
  datExpr <- as.data.frame(t(allDat[,-c(IgnoreCols)]) )
  
  rownames(datExpr) <- colnames(allDat)[-c(IgnoreCols)]
  colnames(datExpr) <- rownames(allDat)
  
  Group <- Group$Group
  
  WGCNA::enableWGCNAThreads()
  
  # Build WGCNA network
  # Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 12, to = PowerUpper, by = 2))
  
  # Call the network topology analysis function
  sft <- WGCNA::pickSoftThreshold(datExpr, 
                                  powerVector = powers, 
                                  RsquaredCut = RCutoff, 
                                  verbose = 0)
  
  # Plot the results
  sft_pwr_fig <- sft_power_fig(sft, RCutoff)
  
  # Get the soft power
  softPower <- sft$powerEstimate
  
  # Build the adjacency table - use "signed" for proteomics data
  adjacency <- WGCNA::adjacency(datExpr, power = softPower, type = "signed")
  
  # Turn adjacency into topological overlap distance
  TOM <- WGCNA::TOMsimilarity(adjacency)
  dissTOM <- 1 - TOM
  
  
  # Clustering using TOM-based dissimilarity
  proTree <- hclust(as.dist(dissTOM), method = "average")
  
  
  # Module identification using dynamic tree cut
  dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = proTree, 
                                               distM = dissTOM, 
                                               deepSplit = 2, 
                                               pamRespectsDendro = FALSE,
                                               minClusterSize = minModuleSize)
  
  if (verbose == TRUE){
  print("Dynamic tree cut results:")
  print(table(dynamicMods))
  }
  
  # Convert numeric labels into colors
  dynamicColors <- WGCNA::labels2colors(dynamicMods)
  if (verbose == TRUE){
  table(dynamicColors)
  }
  
  # Plot the dendrogram and colors underneath
  # plotDendroAndColors(proTree, dynamicColors, "Dynamic Tree Cut",
  #                     dendroLabels = FALSE, hang = 0.03,
  #                     addGuide = TRUE, guideHang = 0.05,
  #                     main = "Protein dendrogram and module colors")
  
  # Merge clusters
  mergedClust <- WGCNA::mergeCloseModules(datExpr, 
                                          dynamicColors, 
                                          cutHeight = MCutheight, 
                                          verbose = 3)
  
  mergedColors <- mergedClust$colors
  
  mergedMEs <- mergedClust$newMEs
  
  # Rename to moduleColors
  moduleColors <- mergedColors
  
  if (verbose == TRUE){
  print("Modules after merging:")
  print(table(moduleColors))
  }
  
  # dendro_fig <- hush(plotDendroAndColors(proTree, cbind(dynamicColors, moduleColors),
  #                                   c("Dynamic Tree Cut", "Merged dynamic"),
  #                                   dendroLabels = FALSE, hang = 0.03,
  #                                   addGuide = TRUE, guideHang = 0.05))
  
  # Get the module eigenproteins
  MEs <- mergedMEs
  
  rownames(MEs) <- rownames(datExpr)
  
  
  # reorder MEs by color names of modules
  MEs <- MEs[,order(colnames(MEs))]
  
  # Plot module profiles with eigenproteins overlaid
  WGCNAClusterID <- moduleColors
  
  cluster.data = t(datExpr)
  group = Group
  
  all_clust_figs <- all_clusters_fig(WGCNAClusterID, 
                                     cluster.data, 
                                     moduleColors, 
                                     meta, 
                                     var)
  
  
  # Eigenprotein network
  
  # eigen_network_fig <- plotEigengeneNetworks(MEs, "Eigenprotein Network",
  #                                            marHeatmap = c(3,4,2,2),
  #                                            marDendro = c(3,4,2,5),
  #                                            plotDendrograms = TRUE,
  #                                            xLabelsAngle = 90,
  #                                            heatmapColors=blueWhiteRed(50),
  #                                            excludeGrey = TRUE)
  
  #The grey module comprises the set of genes which have not been clustered in any module
  MEs <- WGCNA::removeGreyME(MEs)
  
  MEs2 <- MEs
  
  # Heatmap
  heat_fig <- heatmap_wgcna(MEs2, meta, var)
  
  
  # Signed KMEs
  kmes <- WGCNA::signedKME(datExpr, MEs)
  
  # Separate results by modules, order by kME, hub proteins on top
  dat.res <- data.frame(allDat, moduleColors , kmes)
  
  if ("grey" %in% moduleColors){
    cols <- moduleColors[-which(moduleColors == "grey")]
  } else {
    cols <- moduleColors
  }
  
  
  list.clusters <- list()
  for (i in 1:length(levels(as.factor(cols)))){
    
    col <- levels(as.factor(cols))[i]
    
    id <- paste0("kME", col)
    
    temp <- dat.res %>% 
      dplyr::select(Protein.IDs:moduleColors,
                    all_of(id)) %>% 
      dplyr::filter(moduleColors == col) %>% 
      dplyr::arrange(desc(.data[[id]]))
    
    # prot.ids <- clusterProfiler::bitr(temp$Protein.IDs, 
    #                                   fromType = "UNIPROT", 
    #                                   toType = "SYMBOL", 
    #                                   OrgDb = org) %>% 
    #   dplyr::distinct(UNIPROT, .keep_all = TRUE) %>% 
    #   dplyr::rename(Protein.IDs = UNIPROT)
    # 
    # temp <- temp %>%
    #   dplyr::left_join(prot.ids)
    
    list.clusters[[i]] <- temp
    names(list.clusters)[i] <- col
  }
  
  
  # Boxplot for eigenproteins
  all_boxplot_eign_fig <- boxplots_eigenprot(MEs, 
                                             meta, 
                                             var, 
                                             moduleColors)
  
  
  # Boxplot for top 6 hub proteins
  boxplot_hub_prots_fig <- boxplot_hub_prots(MEs, 
                                             list.clusters, 
                                             meta, 
                                             var, 
                                             moduleColors)
  
  # Module-trait correlation analysis
  # Traits
  datTraits <- traits[, -1]
  rownames(datTraits) <- traits$Sample
  
  # Define numbers of genes and samples
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)
  
  
  MEs0 <- WGCNA::orderMEs(MEs)
  modTraitCor <- WGCNA::cor(MEs0, datTraits, use = "p")
  modTraitP <- WGCNA::corPvalueStudent(modTraitCor, nSamples)
  
  
  # Plot
  corr_fig <- cor_plot_wgcna(modTraitP, modTraitCor)
  
  
  # Return both figures and list of proteins per cluster
  figs <- mget(ls(pattern = "_fig"))
  
  all.dat <- list(figs = figs,
                  list.clusters = list.clusters)
  
  return(all.dat)
}



## Plotting functions for wgcna

sft_power_fig <- function(sft, RCutoff){
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  p1 <- sft$fitIndices %>%
    as.data.frame %>% 
    dplyr::mutate(signo = -sign(slope) * SFT.R.sq) %>% 
    ggplot2::ggplot(ggplot2::aes(Power, signo)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = RCutoff,
                        linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(label = Power),
                       nudge_y = .03) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Soft Threshold (power)") +
    ggplot2::ylab("Scale Free Topology Model Fit, signed R^2") +
    ggplot2::labs(title = "Scale independence")
  
  
  # Mean connectivity as a function of the soft-thresholding power
  p2 <- sft$fitIndices %>%
    as.data.frame %>%
    ggplot2::ggplot(ggplot2::aes(Power, mean.k.)) +
    ggplot2::geom_point() +
    ggplot2::geom_text(ggplot2::aes(label = Power),
                       nudge_y = 5) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Soft Threshold (power)") +
    ggplot2::ylab("Mean Connectivity") +
    ggplot2::labs(title = "Mean Connectivity")
  
  # arrange the two figures together
  sft_power_fig <- ggpubr::ggarrange(p1, p2,
                                     align = "h")
  
  return(sft_power_fig)
  
}


all_clusters_fig <- function(WGCNAClusterID, cluster.data, moduleColors, meta, var){
  
  clust_fis <- list()
  
  for (i in 1:length(unique(WGCNAClusterID))){
    
    foo = mean(cluster.data, na.rm = TRUE)
    
    col = unique(WGCNAClusterID)[i]
    
    
    stat_clust <- cluster.data %>%
      tibble::as_tibble() %>% 
      dplyr::mutate(cols = moduleColors,
                    Prots = rownames(cluster.data)) %>% 
      dplyr::filter(cols == col) %>% 
      tidyr::pivot_longer(!c(cols, Prots), 
                          names_to = "key", 
                          values_to = "vals") %>% 
      dplyr::left_join(dplyr::select(meta, 
                                     key, 
                                     dplyr::all_of(var))) %>% 
      dplyr::group_by(.data[[var]]) %>% 
      dplyr::summarise(mean2 = mean(vals, na.rm = TRUE),
                       sd2 = sd(vals, na.rm = TRUE)) %>% 
      dplyr::mutate(grupo = "grupo")
    
    
    fig <- cluster.data %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(cols = moduleColors,
                    Prots = rownames(cluster.data)) %>% 
      dplyr::filter(cols == col) %>% 
      tidyr::pivot_longer(!c(cols, Prots), 
                          names_to = "key", 
                          values_to = "vals") %>% 
      dplyr::left_join(dplyr::select(meta, 
                                     key, 
                                     dplyr::all_of(var))) %>% 
      dplyr::group_by(Prots, .data[[var]]) %>% 
      dplyr::summarise(mean = mean(vals, na.rm = TRUE)) %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[[var]], 
                                   y = mean)) +
      ggplot2::geom_line(ggplot2::aes(group = Prots), 
                         color = "grey80") +
      ggplot2::geom_line(data = stat_clust, 
                         ggplot2::aes(y = mean2, 
                                      x = .data[[var]], 
                                      group = grupo), 
                         color = col, 
                         linewidth = 1.5) +
      ggplot2::geom_linerange(data = stat_clust, 
                              ggplot2::aes(ymin = mean2-sd2, 
                                           ymax = mean2+sd2, 
                                           y = mean2), 
                              linewidth = 1.2) +
      ggplot2::geom_point(data = stat_clust, 
                          ggplot2::aes(y = mean2, 
                                       x = .data[[var]], 
                                       group = grupo), 
                          color = col, 
                          size = 3) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = foo), 
                          linetype = 2, 
                          linewidth = 1) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = NA,
                                                              colour = "black")) +
      ggplot2::scale_x_discrete(expand = c(0.01,0.01)) +
      ggplot2::labs(y = "Average log2-intensities",
                    x = "",
                    title = paste("ME", 
                                  col, 
                                  " (", 
                                  length(moduleColors[moduleColors == col]), 
                                  " proteins)", 
                                  sep = ""))
    
    clust_fis[[i]] <- fig
    names(clust_fis)[i] <- col
    
  }
  
  all_clust_figs <- ggpubr::ggarrange(plotlist = clust_fis)
  
  return(all_clust_figs)
  
}


heatmap_wgcna <- function(MEs2, meta, var){
  
  # Coloring legend
  
  if (length(unique(meta$Condition)) > length(feathers::get_pal("oriole"))){
    pal <- rainbow(length(coefs))
  } else {
    pal <- sample(feathers::get_pal("oriole"))
  }
  
  # create named vector: colors - conditions
  col1 = c()
  for (i in 1:length(unique(meta$Condition))){
    
    col1[i] = pal[i]
    names(col1)[i] <- unique(meta$Condition)[i]
    
  }
  
  
  if (var != "Condition"){
    
    # coloring for extra variable
    col2 <- c()
    
    for (i in 1:length(unique(meta[[var]]))){
      
      col2[i] <- feathers::get_pal("cassowary")[i]
      names(col2)[i] <- unique(meta[[var]])[i]
      
    }
    
    col.ha <- list(Condition = col1,
                   col2)
    names(col.ha)[2] <- eval(var)
    
    # Annotation columns
    df_anno <- data.frame(Condition = meta$Condition,
                          meta[[var]])
    colnames(df_anno)[2] <- eval(var)
    
    # Legend title
    anno_legend <- list(
      Condition = list(title_position = "lefttop-rot",
                       legend_direction = "vertical")
      ,
      list(title_position = "lefttop-rot",
           legend_direction = "vertical"
      ))
    
    names(anno_legend)[2] <- eval(var)
    
  } else {
    
    col.ha <- list(Condition = col1)
    
    # Annotation columns
    df_anno <- data.frame(Condition = meta$Condition)
    
    # Legend title
    anno_legend <- list(
      Condition = list(title_position = "lefttop-rot",
                       legend_direction = "vertical")
    )
    
  }
  
  
  # annotation object
  ha = ComplexHeatmap::HeatmapAnnotation(df = df_anno,
                                         which = 'col',
                                         col = col.ha,
                                         annotation_legend_param = anno_legend)
  
  
  # Heatmap
  heat_fig <- ComplexHeatmap::Heatmap(t(MEs2), 
                                      heatmap_legend_param = list(
                                        title = "MEs",
                                        title_position = "lefttop-rot",
                                        legend_direction = "vertical"
                                      ),
                                      top_annotation = ha)
  
  return(heat_fig)
  
}


boxplots_eigenprot <- function(MEs, meta, var, moduleColors){
  
  # Boxplot for eigenproteins
  boxplot_eigen <- list()
  
  for (i in 1:length(colnames(MEs))){
    
    col_me = colnames(MEs)[i]
    
    fig <- MEs %>% 
      tibble::rownames_to_column("key") %>% 
      tibble::as_tibble() %>% 
      dplyr::left_join(dplyr::select(meta, 
                                     key, 
                                     dplyr::all_of(var))) %>% 
      dplyr::select(key,
                    dplyr::all_of(c(var,
                                    col_me))) %>% 
      ggplot2::ggplot(ggplot2::aes(.data[[var]], 
                                   .data[[col_me]])) +
      ggplot2::geom_boxplot(fill = gsub("ME", "", col_me),
                            alpha = 0.6) +
      ggplot2::theme_classic() +
      ggplot2::labs(y = "Log ratio",
                    x = "",
                    title = paste(col_me, 
                                  " (", 
                                  length(moduleColors[moduleColors == gsub("ME", 
                                                                           "", 
                                                                           col_me)]),
                                  ")", 
                                  sep = ""))
    
    boxplot_eigen[[i]] <- fig
    names(boxplot_eigen)[i] <- col_me
    
  }
  
  all_boxplot_eign_fig <- ggpubr::ggarrange(plotlist = boxplot_eigen)
  
  return(all_boxplot_eign_fig)
  
}


boxplot_hub_prots <- function(MEs, list.clusters, meta, var, moduleColors){
  
  # Boxplot for top 6 hub proteins
  boxplot_hub_prots <- list()
  
  for (i in 1:length(colnames(MEs))){
    
    col_me = paste0("k", colnames(MEs)[i])
    
    col = gsub("kME", "", col_me)
    
    fig <- list.clusters[[col]] %>%
      tibble::as_tibble() %>% 
      dplyr::slice(1:6) %>% 
      tidyr::pivot_longer(!c(Protein.IDs, 
                             dplyr::all_of(col_me), 
                             moduleColors), 
                          names_to = "key", 
                          values_to = "vals") %>%
      dplyr::mutate(prot_kme = paste0(Protein.IDs, 
                                      " (kME = ", 
                                      round(.data[[col_me]], 2), 
                                      ")")) %>% 
      dplyr::left_join(dplyr::select(meta, 
                                     key, 
                                     dplyr::all_of(var))) %>% 
      ggplot2::ggplot(ggplot2::aes(.data[[var]], vals)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::facet_grid(. ~ prot_kme) + 
      ggplot2::labs(y = "Log ratio",
                    x = "",
                    title = col_me)
    
    boxplot_hub_prots[[i]] <- fig
    names(boxplot_hub_prots)[i] <- col_me
    
  }
  
  return(boxplot_hub_prots)
  
}


cor_plot_wgcna <- function(modTraitP, modTraitCor){
  
  modTraitP2 <- modTraitP %>% 
    as.data.frame %>% 
    tibble::rownames_to_column("module") %>% 
    tidyr::pivot_longer(!module, 
                        names_to = "traits", 
                        values_to = "pvals")
  
  
  corr_fig <- modTraitCor %>%
    as.data.frame %>% 
    tibble::rownames_to_column("module") %>% 
    tibble::as_tibble() %>% 
    tidyr::pivot_longer(!module, 
                        names_to = "traits", 
                        values_to = "cor") %>% 
    dplyr::left_join(modTraitP2) %>% 
    ggplot2::ggplot(ggplot2::aes(traits, module, fill = cor)) +
    ggplot2::geom_tile() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(color = "black",
                                                        fill = NA),
                   legend.key.width = ggplot2::unit(0.4, "line"),
                   legend.key.height = ggplot2::unit(2, "line"),
                   legend.direction = "vertical",
                   legend.title = ggplot2::element_text(angle = 90,
                                                        hjust = 1)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::scale_fill_gradient2(low = "blue", 
                                  mid = "white", 
                                  high ="red") +
    ggplot2::geom_text(aes(label = paste0(round(cor,2), 
                                          "\n", 
                                          "(", 
                                          round(pvals, 3),
                                          ")")), 
                       size = 4) +
    ggplot2::labs(x = "",
                  y = "",
                  title = "Module-trait relationships",
                  fill = "Pearson correlation") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(title.position ="left",
                                                    title.hjust = 0.5))
  
  return(corr_fig)
  
}