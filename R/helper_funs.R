prot_log2 <- function(df, meta){
  ## Log2-transform only numeric data, for wrangling or plotting purposes ##
  # df: inputa dataframe. data matrix
  # meta: metadata associated to df
  
  log2_df <- df %>%
    dplyr::mutate(dplyr::across(tidyselect::all_of(meta$key), ~ log2(.x)))
  
  return(log2_df) # log2-transformed data matrix
  
}

cal_z_score <- function(x) {
  ## Performs Z score ## 
  # x: numeric vector
  
  (x - mean(x)) / sd(x)
}