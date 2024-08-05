#' Main KinSwing function
#'
#' @description
#' `kinswing()` runs the four main functions of KinSwing (`buildPWM()`,
#' `scoreSequences()`, `scoreSequences()` and `swing()`) in sequential order,
#' starting with importing the PSP database, creating the annotation dataframe,
#' and running the `topTable()` function from `{limma}` to generate the logFC
#' and p values from all phosphopeptides and comparisons of interest.
#'
#' @param phospho_raw Input dataframe. Raw phosphoproteomic dataframe.
#' @param coefs Input vector. A character vector created with `transform_comp()`
#'   from ProteoMarker.
#' @param fit2 Input list. A list resulting from the last step of `{limma}`
#'   pipeline.
#' @param org Input vector. Character vector of length 1. Refers to the
#'   organism, either "human" or "mouse".
#' @param Protein.IDs Input vector. Character vector of length equivalent to the
#'   number of rows of phospho_raw.
#'
#' @return A list containing the annotated dataframe, the scores list resulting
#'   from the `scoreSequences()` function, and the swing scores dataframe
#'   resulting from the `swing()` function.
#' @export kinswing
kinswing <- function(phospho_raw, coefs, fit2, org, Protein.IDs){
  
  # call the PSP database according to organism
  if(org == "human"){
    
    psp <- KinSwingR::phosphositeplus_human
    
  } else {
    
    psp <- as.matrix(readr::read_rds("D:/PhD/Proteomics/MS_data_PLATONIC/phospho_DBs/PPS_KIN_SUB_mouse_kinswing.rds"))

    psp[,1] <- make.names(psp[,1]) # some kinases are composed of two names separated by a blank space, and if not corrected, the process crashes
    
  }
  
  pwms <- KinSwingR::buildPWM(psp)  
  
  # construct annotation based on phospho raw df
  annotation <- phospho_raw %>% 
    dplyr::rename(Protein.IDs = Index,
                  Gene.names = Gene) %>% 
    dplyr::select(Protein.IDs,
                  Gene.names,
                  SequenceWindow) %>% 
    dplyr::mutate(SequenceWindow = toupper(SequenceWindow),
                  Positions.within.proteins = stringr::str_split(Protein.IDs, 
                                                                 '_', 
                                                                 simplify = TRUE)[,2],
                  Positions.within.proteins = readr::parse_number(Positions.within.proteins),
                  Protein = stringr::str_extract(Protein.IDs, "[^_]+"))
  
  listy <- list()
  n = 100
  
  for (i in 1:length(coefs)){
    
    # get the comparison of interest
    coeff <- coefs[i]
    
    # extract logfc and pval from limma fit2
    annotated_data <- limma::topTable(fit2,
                                      coef = coeff, 
                                      number = Inf, 
                                      genelist = as.data.frame(Protein.IDs)) %>% 
      dplyr::filter(P.Value < .05) %>%
      # dplyr::rename(PPS = anno) %>% 
      dplyr::left_join(annotation, by = "Protein.IDs") %>% 
      dplyr::select(Protein.IDs,
                    Protein,
                    Gene.names,
                    Positions.within.proteins,
                    SequenceWindow,
                    logFC,
                    P.Value) %>%
      dplyr::rename(peptide = SequenceWindow) %>% 
      # KinSwingR only takes windows of 14aa, as found in PSP - take the central 14aa of the original sequence window
      dplyr::mutate(annotation = paste(Protein, 
                                       Gene.names, 
                                       Positions.within.proteins, 
                                       peptide, 
                                       sep = "|")) %>% 
      dplyr::rename(fc = logFC,
                    pval = P.Value) %>% 
      dplyr::select(annotation,
                    peptide,
                    fc,
                    pval)
    
    # Score PWM matches against peptide sequences
    BiocParallel::register(BiocParallel::SnowParam(workers = 4))
    
    set.seed(1234)
    scores <- KinSwingR::scoreSequences(input_data = annotated_data, 
                                        pwm_in = pwms,
                                        n = n,
                                        force_trim = TRUE)
    
    # Predict kinase activity using swing()
    BiocParallel::register(BiocParallel::SnowParam(workers = 4))
    
    set.seed(1234)
    swing_out <- KinSwingR::swing(input_data = annotated_data, 
                                  pwm_in = pwms, 
                                  pwm_scores = scores,
                                  return_network = TRUE)
    
    # Create list for further plotting
    list_temp <- list(scores = scores,
                      annotated_data = annotated_data,
                      swing_scores = swing_out$scores %>% 
                        dplyr::mutate(coef = coeff))
    
    listy[[i]] <- list_temp
    names(listy)[i] <- coeff
    
  }
  
  return(listy)
  
}

#' Volcano plot of predicted kinases
#'
#' @description `kinswing_plot()` returns a ggplot2 object, a volcano plot of
#' the output of the `kinswing()` function.
#'
#'
#' @param ks_list Input list. Output list of `kinswing()` function.
#' @param coefs Input vector. A character vector created with `transform_comp()`
#'   from ProteoMarker.
#'
#' @return A ggplot2 object.
#' @export kinswing_plot
kinswing_plot <- function(ks_list, coefs){
  
  # get df to plot depending on starting list
  if (length(ks_list) > 1){
    
    df_plot <- purrr::list_rbind(lapply(ks_list, `[[`, 3))
    
  } else {
    
    df_plot <- purrr::flatten(ks_list)$swing_scores
    
  }
  
  # extract color pallete for plotting - depends on the number of coefficients
  if (length(coefs) > length(feathers::get_pal("oriole"))){
    pal <- rainbow(length(coefs))
  } else {
    pal <- sample(feathers::get_pal("oriole"))
  }
  
  # create named vector: colors - coefs
  col = c()
  for (i in 1:length(coefs)){
    
    col[i] = pal[i]
    names(col)[i] <- coefs[i]
    
  }
  
  # transform named vector col into dataframe so that it can be merged with df_plot
  pal.df <- tibble::enframe(col, name = "coef", value = "col")
  
  # plot figure
  fig <- df_plot %>% 
    dplyr::left_join(pal.df, by = "coef") %>% 
    dplyr::mutate(p.val = dplyr::case_when(swing < 0 ~ p_less,
                                           .default = p_greater),
                 col = dplyr::case_when((p.val < .05) & (swing > 0) ~ col,
                                        (p.val < .05) & (swing < 0) ~ col,
                                        .default = "#6c6b75"),
                 c.label = dplyr::case_when((p.val < .05) & (swing > 0) ~ kinase,
                                            (p.val < .05) & (swing < 0) ~ kinase,
                                            .default = "")) %>% 
    ggplot2::ggplot(ggplot2::aes(swing, -log10(p.val))) +
    ggplot2::geom_point(ggplot2::aes(color = col, 
                                     size = all), 
                        alpha = .6) +
    ggplot2::scale_color_identity() +
    ggplot2::geom_hline(yintercept = -log10(.05),
                        linetype = 2) +
    hrbrthemes::theme_ipsum(grid = "none",
                            axis_title_size = 12) +
    ggplot2::theme(axis.line.x = ggplot2::element_line(color="black", 
                                                       size = .7),
                   axis.line.y = ggplot2::element_line(color="black", 
                                                       size = .7),
                   legend.position = "bottom") +
    ggplot2::labs(x = "Predicted Ser/Thr/Tyr PK activity",
                  y = "-Log10(P value)",
                  size = "# substrates") +
    ggrepel::geom_label_repel(ggplot2::aes(label = c.label),
                              size = 2.9, # font size in the text labels
                              point.padding = 0, # additional padding around each point
                              min.segment.length = 0, # draw all line segments
                              max.time = 1, max.iter = 1e3, # stop after 1 second, or after 100,000 iterations
                              box.padding = 0.3,
                              max.overlaps = Inf) +
    ggplot2::facet_wrap(.~ coef, 
                        labeller = ggplot2::labeller(ag = coefs)) # if it's only one, you only get one figure
  
  return(fig)
  
}

#' Summarized results from KinSwing analysis
#'
#' @description `kinswing_datable()` generates a list with the results of the
#' KinSwing analysis, that is, it reports the predicted phosphopeptide-kinase
#' relationships.
#'
#'
#' @param ks_list Input list. Output list of `kinswing()` function.
#'
#' @return A list of dataframe(s) with the predicted phosphopeptide-kinase
#'   relationships per comparison.
#' @export kinswing_datable
kinswing_datable <- function(ks_list){
  # extracts kinswing result table
  
  # get information on the predicted kinases
  kinases <- purrr::list_rbind(lapply(ks_list, `[[`, 3)) %>%
    dplyr::mutate(p.val = dplyr::case_when(swing < 0 ~ p_less,
                                           .default = p_greater)) %>% 
    # dplyr::filter(p.val < .05) %>% 
    dplyr::select(kinase,
                  coef,
                  p.val)
  
  
  # merge info on phosphopeptides (logfc and pval of the comparison) with the ks result (given by the kinases df)
  res_kin <- list()
  
  for (i in 1:length(ks_list)){
    
    # pivot kinase info
    foo <- ks_list[[i]]$scores$peptide_p %>% 
      tidyr::pivot_longer(!c("annotation", "peptide"), 
                          names_to = "kinase", 
                          values_to = "swing_pval") %>% 
      dplyr::filter(swing_pval < .05)
    
    # join with peptide scores, parse annotation, order by kinase
    foo2 <- ks_list[[i]]$scores$peptide_scores %>% 
      tidyr::pivot_longer(!c("annotation", "peptide"), 
                          names_to = "kinase", 
                          values_to = "swing_score") %>% 
      dplyr::right_join(foo) %>% 
      dplyr::right_join(ks_list[[i]]$annotated_data) %>% 
      dplyr::filter(pval < .05) %>%
      tidyr::separate_wider_delim(annotation, 
                                  delim = "|",
                                  names = c("uniprot_id",
                                            "symbol",
                                            "phosphosite",
                                            "phosphopeptide")) %>% 
      dplyr::mutate(phosphosite = paste0(phosphosite,
                                         stringr::str_sub(phosphopeptide, 8, 8))) %>% 
      dplyr::select(-peptide,
                    -swing_score,
                    -swing_pval) %>% 
      dplyr::relocate(uniprot_id:phosphopeptide) %>% 
      na.omit() %>% 
      dplyr::arrange(kinase) %>% 
      dplyr::left_join(dplyr::filter(kinases, 
                                     coef == names(ks_list)[i])) %>% 
      dplyr::rename(pval_kinase = p.val) %>% 
      dplyr::select(-coef)
    # dplyr::filter(kinase %in% dplyr::filter(kinases,
    #                                           coef == names(ks_list)[i])$kinase)
    
    res_kin[[i]] <- foo2
    names(res_kin)[i] <- names(ks_list)[i]
    
  }
  
  return(res_kin)
  
}
