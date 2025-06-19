#' limma wrapper function
#'
#' @description `limma_de()` is a wrapper function of limma, and performs the
#'   whole pipeline.
#'
#' @param df.norm input dataframe or matrix. Ideally it is normalized, and it
#'   will be subject to log2 transformation
#' @param meta dataframe with the metadata
#' @param comp dataframe with the comparisons of interest. Must have a Condition
#'   and Control variables
#' @param org organism. Either "human", "mouse" or "yeast"
#' @param phospho is it a phosphoproteomics analysis or not (default "no")
#' @param covars input character vector indicating blocking/confounding
#'   variable(s), default is NULL
#' @param anova input logical to indicate whether topTable should output results
#'   from each comparison or as an ANOVA
#'
#' @returns A list of dataframes, one for each comparison (or ANOVA) with the
#'   results yielded by topTable
#' @export limma_de
limma_de <- function(df.norm, meta, comp, org, phospho = "no", covars = NULL, anova){

  # Transform files
  comp.files <- prep_dfs_limma(df.norm, meta)
  comp.v <- transform_comp(comp)
  
  # Make sure that there are no proteins with 100% NAs
  comp.files$df.comp <- comp.files$df.comp %>%
    dplyr::select(dplyr::any_of(meta$key)) %>%
    janitor::remove_empty(., which = "rows")
  
  # Create model matrix
  mod <- complex_model_matrix(comp.files$df.comp,
                              comp,
                              comp.v,
                              comp.files$meta,
                              covars)
  
  
  # We use limma::lmFit to fit the model Fit the model
  fit <- limma::lmFit(dplyr::select(comp.files$df.comp,
                                    tidyselect::any_of(comp.files$meta$key)),
                      mod)
  
  # Comparisons - only when there is no intercept
  if (!"Intercept" %in% colnames(mod)){
    
    # Create contrast matrix
    contrast.matrix <- limma::makeContrasts(contrasts = comp.v,
                                            levels = mod)
    
    # Fit the comparisons
    fit <- limma::contrasts.fit(fit, contrast.matrix)
  }
  
  
  # Test these coefficients being != 0 with limma::eBayes
  fit2 <- limma::eBayes(fit, trend = TRUE)
  
  # Get results table:
  if (!"Intercept" %in% colnames(mod)){ # So that topTable gets relevant coefs only
    coefs <- colnames(contrast.matrix)
  } else {
    coefs <- colnames(mod)[which(colnames(mod) %in% as.vector(as.matrix(comp)))]
  }
  
  # For annotation
  org <- switch(org,
                "human" = "org.Hs.eg.db",
                "mouse" = "org.Mm.eg.db",
                "yeast" = "org.Sc.sgd.db")
  
  listy <- get_limma_res(fit2, 
                         comp.v, 
                         coefs, 
                         comp.files$df.comp, 
                         org, 
                         phospho, 
                         anova)
  
  # Calculate standard error of effect sizes / logFC
  SE_df <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
  
  SE_df <- SE_df %>% 
    as.data.frame %>% 
    tibble::rownames_to_column("Protein.IDs")
  
  # Extract degrees of freedom (total)
  d_freedom <- data.frame(Protein.IDs = SE_df$Protein.IDs,
                          dfs = fit2$df.total)
  
  
  # Incorporate everything
  for (i in 1:length(listy)){
    
    nam_l <- names(listy)[i]
    
    nam <- colnames(SE_df)[sapply(colnames(SE_df), grepl, nam_l)]
    
    listy[[i]] <- listy[[i]] %>% 
      dplyr::left_join(dplyr::select(SE_df,
                                     Protein.IDs,
                                     dplyr::any_of(nam)),
                       by = "Protein.IDs") %>% 
      dplyr::rename(logFC_SE = dplyr::any_of(nam)) %>%
      dplyr::relocate(logFC_SE, .after = logFC) %>% 
      dplyr::left_join(d_freedom,
                       by = "Protein.IDs")
    
  }
  
  # checkpoint - if differences between calculated logFCs and those yielded by limma - manually check
  
  if (!isTRUE(anova)){
    check_log2(df.norm, listy, meta)
  }
  
  return(listy)
  
}

#' Transform comparison file
#'
#' @description `transform_comp()` transforms the comparison file so that it is ready
#'   to use by limma
#'
#' @param comp dataframe with the comparisons of interest. Must have a Condition
#'   and Control variables
#'
#' @return A character vector with all the comparisons
#'
#' @noRd
transform_comp <- function(comp){
  
  if (all(colnames(comp) == c("Condition", "Control"))){
    
    comp.v <- comp %>%
      dplyr::mutate(comp = stringr::str_c(Condition, Control, sep = "-")) %>%
      dplyr::select(comp) %>%
      dplyr::pull()
    
  } else {
    
    stop("Variable names of comparison file should be 'Condition' and 'Control'.")
    
  }
  
  return(comp.v)
  
}

#' Transform data and metadata
#'
#' @description Function that transforms the data and medatada files so that it
#'   is ready to use by limma
#'
#' @param df input dataframe or matrix. Ideally it is normalized, and it will be
#'   subject to log2 transformation
#' @param meta dataframe with the metadata
#'
#' @return A list with the transformed data and metadata
#'
#' @noRd
prep_dfs_limma <- function(df, meta){
  
  # Data
  df.comp <- prot_log2(df, meta) %>%
    dplyr::select(!tidyselect::any_of(meta$key[meta$Condition == "Norm"])) %>%
    as.data.frame()
  
  rownames(df.comp) <- df.comp$Protein.IDs
  
  # Metadata vars
  meta.vars <- meta %>%
    dplyr::select(-BioReplicate,
                  -key) %>%
    colnames()
  
  # Metadata
  meta <- meta %>%
    dplyr::filter(!Condition == "Norm") %>%
    dplyr::mutate(across({{ meta.vars }}, ~ as.factor(.x)))
  
  
  listy <- list(df.comp = df.comp,
                meta = meta)
  
  return(listy)
  
}

#' Create model matrix
#'
#' @description `complex_model_matrix()` creates a model matrix for limma
#'
#' @param df.comp dataframe with transformed data, as output of prep_dfs_limma()
#' @param comp dataframe with the comparisons of interest. Must have a Condition
#'   and Control variables
#' @param comp.v input vector with comparisons, as output of transform_comp()
#' @param meta dataframe with the metadata
#' @param covars input character vector indicating blocking/confounding
#'   variable(s), default is NULL
#'
#' @return A model matrix
#'
#' @noRd
complex_model_matrix <- function(df.comp, comp, comp.v, meta, covars = NULL){
  
  # TODO: Implement more complex experimental designs

  # Get name of the grouping variable
  grouping.var <- grouping_var(meta, comp)
  
  # Get the variables that compose the formula
  # Decision tree for whether there are covariables (blocking) / paired analysis
  if (is.null(covars)){
    vars <- grouping.var
  } else {
    vars <- paste(grouping.var, " + ", paste0(covars, collapse = " + "))
  }
  
  
  # Decision tree for the intercept
  if (length(unique(comp$Control)) == 1){
    xnam <- paste0("~ ", vars) # If the comparison is against a control/baseline: intercept
    
    # Make sure that the control in the first level factor
    levs <- c(unique(comp$Control), unique(comp$Condition))
    
    if (length(setdiff(unique(meta$Condition), levs)) > 0){
      # Not all levels of Condition might want to be compared, and thus be present in the comparison file, so we have to make sure that they are there
      levs <- c(levs, setdiff(unique(meta$Condition), levs))
      
    }
    
    meta <- meta %>%
      dplyr::mutate(Condition = factor(Condition,
                                       levels = levs))
    
  } else {
    xnam <- paste0("~0 + ", vars) # To get all comparisons
  }
  
  
  # Transform to class formula
  fmla <- as.formula(xnam)
  
  mod <- model.matrix(fmla, data = meta) # data is always meta
  
  # Set rownames
  rownames(mod) <- meta$key
  
  # Clean colnames
  colnames(mod) <- gsub("[()]", "", colnames(mod))
  colnames(mod) <- gsub("Condition", "", colnames(mod))
  
  colnames(mod) <- gsub('\\B([[:upper:]])', ' \\1', colnames(mod))
  colnames(mod) <- stringr::word(colnames(mod), -1) # safer to get the last word, instead of removing the first (there could be more than one)
  # colnames(mod) <- stringr::str_remove(colnames(mod), '(\\w+\\s+){1}') # removes the first word
  
  return(mod)
  
}

#' topTable wrapper and annotation for phospho data
#'
#' @description `get_limma_res()` performs the last step of limma_de() depending
#'   on whether it's phospho data or not
#'
#' @param fit2 fitted model from limma_de()
#' @param comp.v input vector with comparisons, as output of transform_comp()
#' @param coefs comparisons as stipulated in the model matrix
#' @param df dataframe as input to limma_de(), yielded by prep_dfs_limma()
#' @param org organism. Either "human", "mouse" or "yeast"
#' @param phospho is it a phosphoproteomics analysis or not (default "no")
#' @param anova input logical to indicate whether topTable should output results
#'   from each comparison or as an ANOVA
#'
#' @return A list of dataframes, one for each comparison (or ANOVA)
#'
#' @noRd
get_limma_res <- function(fit2, comp.v, coefs, df, org, phospho = "no", anova){

  listy <- list()
  Protein.IDs <- rownames(df)
  
  if (phospho == "yes"){
    
    ids <- sub("(^[^-]+)-.*", "\\1", Protein.IDs)
    ids <- sub("(^[^_]+)_.*", "\\1", Protein.IDs)
    
    id <- clusterProfiler::bitr(ids,
                                OrgDb = org,
                                fromType = "UNIPROT",
                                toType = c("SYMBOL", "GENENAME", "ENTREZID")) %>%
      dplyr::distinct(UNIPROT, .keep_all = TRUE)
    
    if (isTRUE(anova)){
      
      listy <- get_res(fit2, anova, coefs[i], Protein.IDs, id, phospho = "yes")
      
    } else {
      
      for (i in 1:length(coefs)){
        
        res <- get_res(fit2, anova, coefs[i], Protein.IDs, id, phospho = "yes")
        
        listy[[i]] <- res
        names(listy)[i] <- comp.v[i]
        
      }
    }
    
  } else {
    
    if (isTRUE(anova)){
      
      listy <- get_res(fit2, anova, coefs[i], Protein.IDs, id)
      
    } else {
      
      for (i in 1:length(coefs)){
        
        res <- get_res(fit2, anova, coefs[i], Protein.IDs, id)
        
        listy[[i]] <- res
        names(listy)[i] <- comp.v[i]
        
      }
    }
  }
  
  # Remove NULL elements from listy
  listy <- purrr::compact(listy)
  
  return(listy)
  
}

#' topTable wrapper
#'
#' @description `get_res()` performs topTable() with some modifications to the
#'   foldchange and adjusted P value
#'
#' @param fit2 fitted model from limma_de()
#' @param anova input logical to indicate whether topTable should output results
#'   from each comparison or as an ANOVA
#' @param coefs comparisons as stipulated in the model matrix
#' @param Protein.IDs input vector with protein annotation
#' @param id dataframe with extra annotation for phospho analysis, default is
#'   NULL
#' @param phospho is it a phosphoproteomics analysis or not (default "no")
#'
#'
#' @return A list of dataframes, one for each comparison (or ANOVA)
#'
#' @noRd
get_res <- function(fit2, anova, coefs, Protein.IDs, id = NULL, phospho = "no"){
  
  if (isTRUE(anova)){
    
    res <- limma::topTable(fit2,
                           number = Inf,
                           genelist = as.data.frame(Protein.IDs),
                           confint = T) %>%
      dplyr::mutate(direction = dplyr::case_when(adj.P.Val < 0.05 ~ "DE",
                                                 TRUE ~ "No DE"),
                    log.adj.P.Val = -1 * log10(adj.P.Val))
    
  } else {
    
    res <- limma::topTable(fit2,
                           coef = coefs,
                           number = Inf,
                           genelist = as.data.frame(Protein.IDs),
                           confint = T) %>%
      dplyr::mutate(direction = dplyr::case_when((logFC > 0) & (adj.P.Val < 0.05) ~ "Up",
                                                 (logFC < 0) & (adj.P.Val < 0.05) ~ "Down",
                                                 TRUE ~ "No DE"),
                    log.adj.P.Val = -1 * log10(adj.P.Val))
    
  }
  
  
  if (phospho == "yes"){
    
    res <- res %>%
      dplyr::mutate(UNIPROT = sub("(^[^_]+)_.*", "\\1", Protein.IDs)) %>%
      dplyr::left_join(id, by = "UNIPROT")
    
  }
  
  return(res)
  
}


#' Formatting of results from limma
#'
#' @description `result_table()` generates a dataframe with a summary of the
#'   number of significantly changing proteins, both upregulated and
#'   downregulated, along with the number of those that do not change. The
#'   cutoff is based on the adjusted P value (< 0.05)
#'
#' @param res input list resulting from limma_de()
#'
#' @returns A dataframe
#' @export result_table
result_table <- function(res){

  if (is.list(res)){
    res <- dplyr::bind_rows(res, .id = "comp")
  }
  
  res.table <- res %>%
    dplyr::select(comp,
                  direction) %>%
    dplyr::group_by(comp) %>%
    dplyr::count(direction) %>%
    dplyr::mutate(direction = factor(direction,
                                     levels = c("Down", "No DE", "Up"))) %>%
    tidyr::pivot_wider(names_from = comp, values_from = n) %>%
    dplyr::mutate_if(is.numeric, ~replace(., is.na(.), 0))
  
  
  if (nrow(res.table) < 3){ # If there are no up and/or down regulated, add and make them 0
    
    levs.out <- setdiff(levels(res.table$direction), res.table$direction)
    
    foo <- list()
    
    for (i in 1:length(levs.out)){
      temp <- c(levs.out[i], rep(0, ncol(res.table)-1))
      foo[[i]] <- temp
    }
    
    res.table <- rbind(setNames(as.data.frame(do.call(rbind, foo)), colnames(res.table)), res.table)
    
    res.table <- plyr::ddply(res.table, c('direction')) # order based on the levels of direction
    
  }
  
  res.table <- res.table %>% t %>% as.data.frame() %>% janitor::row_to_names(1)
  
  return(res.table)
}

#' Checkpoint for limma
#'
#' @description `check_log2()` makes sure that the differential expression
#'   analysis took place with samples allocated to the correct groups (i.e.,
#'   makes sure the order of the samples is correct). For this, it calculates
#'   the fold-change separately, and then compares it to the one calculated with
#'   limma_de(), to see if they match.
#'
#' @param df input dataframe or matrix. Ideally it is normalized, and it will be
#'   subject to log2 transformation
#' @param res_limma output list from get_limma_res()
#' @param metadata dataframe with the metadata
#'
#' @return A string
#'
#' @noRd
check_log2 <- function(df, res_limma, metadata){
  
  for (i in 1:length(res_limma)){
    
    if (grepl("\\)", names(res_limma)[i])){
      
      foo1 <- gsub("\\(", "[", names(res_limma)[i])
      foo2 <- gsub("\\)", "]", foo1)
      
      foo3 <- unlist(strsplit(foo2, "-(?![^[]*\\])", perl=TRUE))
      
      foo3 <- gsub("\\[|\\]", "", foo3)
      
      name_comp1 <- strsplit(foo3, "-", perl = T)
      
      list_lfc <- list()
      
      for (j in 1:length(name_comp1)){
        
        name_comp <- name_comp1[[j]]
        
        meta_comp <- metadata %>%
          dplyr::filter(Condition %in% name_comp)
        
        
        df_org <- df %>%
          prot_log2(., metadata) %>%
          dplyr::select(Protein.IDs,
                        meta_comp$key) %>%
          tidyr::pivot_longer(-Protein.IDs,
                              names_to = "key",
                              values_to = "vals") %>%
          dplyr::left_join(dplyr::select(metadata,
                                         Condition,
                                         key),
                           by = "key") %>%
          dplyr::group_by(Protein.IDs,
                          Condition) %>%
          dplyr::summarise(mean_v = mean(vals, na.rm = TRUE),
                           .groups = "drop") %>%
          tidyr::pivot_wider(names_from = Condition,
                             values_from = mean_v) %>%
          dplyr::mutate(logFC = .data[[name_comp[1]]] - .data[[name_comp[2]]]) %>% # log2FC
          dplyr::select(Protein.IDs,
                        logFC)
        
        
        list_lfc[[j]] <- df_org
        
      }
      
      
      df_org <- purrr::reduce(list_lfc, dplyr::full_join, by = "Protein.IDs") %>% 
        dplyr::mutate(logFC = logFC.x - logFC.y) %>% 
        dplyr::select(Protein.IDs,
                      logFC)
      
      
      temp <- res_limma[[i]] %>%
        dplyr::arrange(match(Protein.IDs, df_org$Protein.IDs))
      
      df_org <- df_org %>%
        dplyr::filter(Protein.IDs %in% temp$Protein.IDs) # make sure dims are the same for the all.equal
      
      
      if (all(df_org$Protein.IDs == temp$Protein.IDs)){ # ensures the Protein.ID order is the same
        
        if (!isTRUE(all.equal(df_org$logFC, temp$logFC))){
          
          print(paste(names(res_limma)[i], "-",
                      all.equal(df_org$logFC, temp$logFC)))
        } else {
          print(paste(names(res_limma)[i], "-> log2FC are OK"))
        }
      }
      
    } else {
      
      # extract the 2 condition names to make sure the logFC is calculated correctly
      name_comp <- c(stringr::str_split(names(res_limma)[i],
                                        "-",
                                        simplify = TRUE))
      
      meta_comp <- metadata %>%
        dplyr::filter(Condition %in% name_comp)
      
      df_org <- df %>%
        prot_log2(., metadata) %>%
        dplyr::select(Protein.IDs,
                      meta_comp$key) %>%
        tidyr::pivot_longer(-Protein.IDs,
                            names_to = "key",
                            values_to = "vals") %>%
        dplyr::left_join(dplyr::select(metadata,
                                       Condition,
                                       key),
                         by = "key") %>%
        dplyr::group_by(Protein.IDs,
                        Condition) %>%
        dplyr::summarise(mean_v = mean(vals, na.rm = TRUE),
                         .groups = "drop") %>%
        tidyr::pivot_wider(names_from = Condition,
                           values_from = mean_v) %>%
        dplyr::mutate(logFC = .data[[name_comp[1]]] - .data[[name_comp[2]]]) %>% # log2FC
        dplyr::select(Protein.IDs,
                      logFC)
      
      temp <- res_limma[[i]] %>%
        dplyr::arrange(match(Protein.IDs, df_org$Protein.IDs))
      
      df_org <- df_org %>%
        dplyr::filter(Protein.IDs %in% temp$Protein.IDs) # make sure dims are the same for the all.equal
      
      if (all(df_org$Protein.IDs == temp$Protein.IDs)){ # ensures the Protein.ID order is the same
        
        if (!isTRUE(all.equal(df_org$logFC, temp$logFC))){
          
          print(paste(paste0(name_comp, collapse = "-"), "-",
                      all.equal(df_org$logFC, temp$logFC)))
        } else {
          print(paste(paste0(name_comp, collapse = "-"), "-> log2FC are OK"))
        }
      }
    }
    
  }
  
}

#' Model matrix function helper
#'
#' @description `grouping_var()` identifies which metadata variable contains the
#'   grouping information (usually Condition)
#'
#' @param metadata dataframe with the metadata
#' @param comp dataframe with the comparisons of interest. Must have a Condition
#'   and Control variables
#'
#' @return A string
#'
#' @noRd
grouping_var <- function(meta, comp){
  
  var_m <- meta %>%
    dplyr::select(dplyr::where(~select_grouping_var(.x, comp))) %>% 
    colnames
  
  if (any(c("BioReplicate", "key") %in% var_m)){
    
    var_m <- var_m[-which(var_m %in% c("BioReplicate", "key"))]
    
  }
  
  return(var_m)
  
}

select_grouping_var <- function(x, comp){
  
  x = paste0("\\b", x, "\\b") # forces exact match
  
  any(sapply(x, grepl, comp$Condition[1]))
  
}