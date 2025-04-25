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