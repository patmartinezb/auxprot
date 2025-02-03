#' Filters dataframe based on missing values
#'
#' @param df dataframe to be filtered
#' @param threshold numeric input indicating the minimum percentage of observed values across all samples. Goes from 0 to 1, default is 0.7 (i.e., 70%)
#'
#' @returns a filtered dataframe
#' 
#' @export filt_na
filt_na <- function(df, threshold = 0.7){
  
  # Filters proteins based on the number of NAs across all samples - threshold refers to non-missing values
  
  filt <- df %>%
    dplyr::mutate(filt = dplyr::case_when(
      rowSums(!is.na(dplyr::across(tidyselect::where(is.numeric)))) / ncol(dplyr::across(tidyselect::where(is.numeric))) >= threshold ~ TRUE,
      TRUE ~ FALSE
    )
    ) %>%
    dplyr::filter(filt == TRUE) %>%
    dplyr::select(-filt)
  
  return(filt)
  
}

prot_log2 <- function(df, meta){
  ## Log2-transform only numeric data, for wrangling or plotting purposes ##
  # df: inputa dataframe. data matrix
  # meta: metadata associated to df
  
  log2_df <- df %>%
    dplyr::mutate(dplyr::across(tidyselect::any_of(meta$key), ~ log2(.x)))
  
  return(log2_df) # log2-transformed data matrix
  
}

cal_z_score <- function(x) {
  ## Performs Z score ## 
  # x: numeric vector
  
  (x - mean(x)) / sd(x)
}