#' Sample loading normalization
#'
#' @description `SL_norm()` normalizes the sum of each channel to the average
#'   grand total
#'
#' @param df dataframe with the raw data (not log2-transformed)
#' @param meta dataframe with metadata
#'
#' @returns a normalized dataframe
#'
#' @export SL_norm
SL_norm <- function(df, meta) {
  # Computes normalization factors and scales columns
  
  # Transform df to long format and add metadata
  df_long <- df %>% 
    tidyr::pivot_longer(!Protein.IDs, 
                        names_to = "key", 
                        values_to = "vals") %>% 
    dplyr::left_join(dplyr::select(meta,
                                   Mixture,
                                   key),
                     by = "key") 
  
  # Calculate the sum of the intensities of each sample
  df_sums <- df_long %>% 
    dplyr::group_by(key) %>% 
    dplyr::summarise(sums = sum(vals, na.rm = T)) %>% 
    dplyr::ungroup()
  
  # Calculate the mean of those sums, by plex
  mean_v <- df_sums %>%
    dplyr::left_join(meta, by = "key") %>% 
    dplyr::group_by(Mixture) %>% 
    dplyr::summarise(mean = mean(sums, na.rm = T))
  
  # For each sample, calculate the normalization factors
  df_norm_fact <- df_sums %>% 
    dplyr::left_join(meta, by = "key") %>%
    dplyr::left_join(mean_v,
                     by = "Mixture") %>% 
    dplyr::mutate(norm_fact = mean / sums) %>% 
    dplyr::select(key,
                  Mixture,
                  norm_fact)
  
  # Multiply each protein in each sample by the corresponding normalization factor
  df_sl <- df_long %>% 
    dplyr::left_join(df_norm_fact,
                     by = c("key", "Mixture")) %>% 
    dplyr::mutate(vals_norm = vals * norm_fact) %>% 
    dplyr::select(- vals,
                  - Mixture,
                  - norm_fact) %>% 
    tidyr::pivot_wider(names_from = key, 
                       values_from = vals_norm)
  
  return(df_sl) # return normalized data
  
}

#' calcNormFactors wrapper
#'
#' @description `tmm_norm()` is a wrapper function for edgeR::calcNormFactors().
#'   It accounts for compositional bias in the data.
#'
#' @param df dataframe with the raw data (not log2-transformed), usually as
#'   output of SL_norm()
#' @param meta dataframe with metadata
#'
#' @returns a normalized dataframe
#'
#' @export tmm_norm
tmm_norm <- function(df, meta){
  
  Protein.IDs <- df$Protein.IDs
  
  # Calculates factors
  tmm_factors <- edgeR::calcNormFactors(dplyr::select(na.omit(df), 
                                                      tidyselect::any_of(meta$key)))
  df_tmm <- sweep(dplyr::select(df, tidyselect::any_of(meta$key)), 
                  2, 
                  tmm_factors, 
                  FUN = "/")
  
  df_tmm <- cbind(Protein.IDs, df_tmm) # Recover annotation
  
  return(df_tmm)
}


#' Internal reference scaling (IRS) normalization
#'
#' @description `irs_norm()` performs internal reference scaling (IRS)
#'   normalization (after sample loading (SL) and TMM (optional) normalizations)
#'   in the presence of 2 or more TMT experiments. It uses one (or several)
#'   reference/pool channels, and in the absence of one, calculates a 'virtual'
#'   one using the information of all proteins of the plexes.
#'
#'
#' @param df Dataframe that is going to be IRS-normalized
#' @param meta Dataframe of associated metadata. Pool channels should be flagged
#'   as 'Norm', but if a single sample is used repeatedly across all TMTs as a
#'   sort of technical replicate that acts as a reference channel, and these
#'   samples should be then collapsed into a single one after normalization,
#'   they should be flagged as 'Ref'.
#' @param tmm TRUE (default) or FALSE, whether the function should perform TMM
#'   normalization after SL normalization and before IRS normalization
#'
#' @returns A normalized dataset
#' @export irs_norm
irs_norm <- function(df, meta, tmm = TRUE){

  # Previous normalizations
  df <- SL_norm(df, meta)
  
  if (tmm == TRUE){
    df <- tmm_norm(df, meta)
  }
  
  # IRS normalization
  # Decision making based on the existence (or not) of a reference channel
  
  if (("Norm" %in% meta$Condition) | (any(grepl("Ref|Norm", meta$key)))){ # There is a ref channel
    
    # Define variable with norm channels so that it's not dependent on the existence of Norm in Condition
    ref_names <- meta %>%
      dplyr::filter(Condition == "Norm" | grepl("Ref|Norm", key)) %>%
      dplyr::pull(key)
    
    print("Ref names:")
    print(paste(ref_names, sep = "\n"))
    
    
    if (exists("ref_names")){ # exists() argument, even if it is a variable, should be a string
      if (!"Norm" %in% meta$Condition){
        
        type_ref = "sample"
        
        name_ref = unique(gsub("(Ref|Norm)[0-9]+_", "", ref_names))
        # Note: this works okay if it is always the same sample (which it sould be)
      }
    }
    
    
    ## Decision tree based on the number of ref channels
    if (length(ref_names) > length(unique(meta$Mixture))){ # >1 ref channel / plex
      
      irs_channels <- df %>%
        dplyr::select(Protein.IDs,
                      tidyselect::any_of(ref_names)) %>% # First select ref channels
        # na.omit() %>%
        tidyr::pivot_longer(!Protein.IDs,
                            names_to = "key",
                            values_to = "vals") %>%
        dplyr::left_join(meta, by = "key") %>%
        dplyr::group_by(Mixture, Protein.IDs) %>%
        dplyr::summarise(mean = exp(mean(log(vals))),
                         .groups = "drop") %>% # Calculate the geometric mean of each group of refs
        tidyr::pivot_wider(names_from = Mixture,
                           values_from = mean)
      
      
      if(exists("ref_names") & ("Norm" %in% meta$Condition)){
        # Remove Norm (pool) samples
        
        df <- df %>%
          dplyr::select(-tidyselect::any_of(ref_names))
        
        meta <- transform_meta_irs(meta)
      }
      
      
    } else { # One ref channel / plex
      
      irs_channels <- df %>%
        dplyr::select(Protein.IDs,
                      tidyselect::any_of(ref_names)) %>%
        tidyr::pivot_longer(!Protein.IDs,
                            names_to = "key",
                            values_to = "vals") %>%
        dplyr::left_join(meta, by = "key") %>%
        dplyr::select(Protein.IDs,
                      Mixture,
                      vals) %>%
        tidyr::pivot_wider(names_from = Mixture,
                           values_from = vals)
      # na.omit()
      
      
      if(exists("ref_names") & ("Norm" %in% meta$Condition)){
        # Remove Norm (pool) samples
        
        df <- df %>%
          dplyr::select(-tidyselect::any_of(ref_names))
        
        meta <- transform_meta_irs(meta)
        
      }
      
    }
    
  } else { # There isn't a ref channel - create one per plex // virtual
    
    # Sum the intensities for each protein, across samples, within each plex
    
    irs_channels <- df %>%
      tidyr::pivot_longer(!Protein.IDs,
                          names_to = "key",
                          values_to = "vals") %>%
      dplyr::left_join(meta, by = "key") %>%
      dplyr::group_by(Mixture, Protein.IDs) %>%
      dplyr::summarise(sum = hablar::sum_(vals, ignore_na = TRUE),
                       .groups = "drop") %>%
      tidyr::pivot_wider(names_from = Mixture,
                         values_from = sum)
    # na.omit()
    
    # View(irs_channels)
    
  }
  
  # Common to everything: calculate IRS factors
  irs_factors <- irs_channels %>%
    dplyr::rowwise() %>%
    dplyr::mutate(geomean = exp(mean(log(dplyr::c_across(dplyr::starts_with("Mixture"))),
                                     na.rm = TRUE))) %>% # Across proteins calculate geometric mean
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ geomean / .x)) %>%
    dplyr::select(-geomean) %>%
    tidyr::pivot_longer(!Protein.IDs,
                        names_to = "mix",
                        values_to = "factors")
  
  # View(irs_factors)
  
  # Apply correction to all df
  df_irs <- df %>%
    tidyr::pivot_longer(!Protein.IDs,
                        names_to = "key",
                        values_to = "vals") %>%
    dplyr::left_join(meta, by = "key") %>%
    dplyr::right_join(irs_factors,
                      by = "Protein.IDs",
                      relationship = "many-to-many") %>%
    dplyr::filter(Mixture == mix) %>% # Make sure each factor matches its corresponding mixture
    dplyr::mutate(irs_vals = factors * vals) %>% # Normalize
    dplyr::select(Protein.IDs,
                  key,
                  irs_vals) %>%
    tidyr::pivot_wider(names_from = key, values_from = irs_vals)
  
  
  if(exists("type_ref")){ # if one sample is repeated in every tmt as ref channel
    
    #   # instead of calculating the median of ref channels before the actual irs normalization, I apply the normalization to all of them, as if they were regular sample channels, and then calculate the median
    
    df_irs <- df_irs %>%
      dplyr::rowwise() %>% # calculate median of all ref samples
      dplyr::mutate("{name_ref}" := median(dplyr::c_across(dplyr::any_of(ref_names)),
                                           na.rm = TRUE)) %>%
      dplyr::select(-tidyselect::any_of(ref_names))
    
  }
  
  return(df_irs)
  
}