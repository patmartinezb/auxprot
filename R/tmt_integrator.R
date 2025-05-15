#' FragPipe's TMT Integrator wrapper
#'
#' @description `tmt_integrator()` mimics the FragPipe's TMT Integrator module
#'   filters for TMT data, minus de PSM normalization. This includes purity and
#'   peptide probability filters, as well as a minimum summed MS2 intensity,
#'   presence of TMT and, optionally, PTM labeling. It summarises PSM
#'   information to the peptide and subsequent protein level, for further
#'   analysis.
#'
#'
#' @param psm Input FragPipe data matrix at the PSM level.
#' @param metadata Input dataframe. Metadata associated to `psm`.
#' @param org Input vector of length 1. Indicates the species the data belongs
#'   to.
#' @param tmt Input vector of length 1. Indicates the type of TMT used, either
#'   `tmt`, `tmtpro` or `tmtpro0`.
#' @param purity_thr Numeric vector of length 1. Indicates purity threshold.
#'   Default is 0.5 (50%).
#' @param pep_prob_thr Numeric vector of length 1. Indicates peptide probability
#'   threshold. Default is 0.9.
#' @param sum_ms2_int_thr Numeric vector of length 1. Indicates the percentage
#'   of MS2 intensity below which PSMs are going to be filtered out. Default is
#'   0.05 (5%). For PTM data, such as phosphorylation, it should be 0.025
#'   (2.5%).
#'
#' @return A list containing both peptide and protein level information
#' @export tmt_integrator
tmt_integrator <- function(psm, 
                           metadata, 
                           org, 
                           tmt,
                           purity_thr = 0.5,
                           pep_prob_thr = 0.9,
                           sum_ms2_int_thr = 0.05){
  
  # TODO: if phosphoproteomics, sum_ms2_int_thr is 0.025 (2.5%), and must contain phosphorylation
  # TODO: use stringr::str_glue() and stringr::str_glue_data() instead of cat()
  # https://glue.tidyverse.org/
  # NOTE: The janitor::make_clean_names(Assigned.Modifications) takes forever when the dataset is large
  
  # 0.Initial transformation and filtering for contaminants of PSM data
  psm <- psm %>% 
    tibble::as_tibble(.name_repair = make.names) %>% # make.names() preserves the sample names as in metadata
    dplyr::rename_with(., # name fix for fragpipe v23
                       ~ gsub("Intensity.", "", .x),
                       dplyr::starts_with("Intensity.")) %>%
    dplyr::mutate(sum_int = rowSums(dplyr::across(dplyr::any_of(metadata$key)), 
                                    na.rm = T)) %>%  # sum all intensitites, rowwise
    dplyr::mutate(Assigned.Modifications = tolower(Assigned.Modifications),
                  Assigned.Modifications = gsub("\\-|\\(|\\)", 
                                                ".", 
                                                Assigned.Modifications),
                  Assigned.Modifications = gsub("\\.", "_", Assigned.Modifications),
                  organism = sub('.*\\_', '', Protein)) %>%  
    dplyr::filter(organism == toupper(org)) # remove contaminants from other species
  
  
  # 1. If there's a ref/pool channel (+1 plexes), no 0 in them
  
  if (any(grepl("Ref", metadata$key))){
    ref_names <- metadata$key[grepl("Ref", metadata$key)]
  }
  
  if(exists("ref_names")){
    
    psm <- psm %>% 
      dplyr::mutate(dplyr::across(dplyr::any_of(metadata$key), 
                                  ~ dplyr::na_if(.x, 0))) %>% 
      tidyr::drop_na(., dplyr::any_of(ref_names))
    
  }
  
  # 2. Precursor ion purity filtering
  # cat("Distribution of the PSMs purity:\n")
  # cat(summary(psm$Purity))
  
  cat("There are", nrow(psm[psm$Purity <= purity_thr,]), "PSMs below the",
      purity_thr, "threshold.\n")
  
  psm1 <- psm %>%
    dplyr::filter(Purity >= purity_thr)
  
  # 3. Filtering by minimum peptide probability
  # cat("Distribution of the PSMs probability:\n")
  # cat(summary(psm1$Probability))
  
  cat("There are", nrow(psm1[psm1$Probability <= pep_prob_thr,]), "PSMs below the",
      pep_prob_thr, "threshold.\n")
  
  psm2 <- psm1 %>%
    dplyr::filter(Probability >= pep_prob_thr)
  
  # 4. Summed MS2 intensity filtering
  five.pc <- quantile(psm$sum_int, probs = sum_ms2_int_thr)
  
  stringr::str_glue("The summed MS2 intensity {sum_ms2_int_thr*100}% threshold is {five.pc}.")
  
  cat("There are", nrow(psm2[psm2$sum_int < five.pc,]), "PSMs below the",
      sum_ms2_int_thr*100, "% threshold.\n")
  
  psm3 <- psm2 %>%
    dplyr::filter(sum_int >= five.pc)
  
  
  # Must be TMT labeled
  # the weight is different if tmt or tmtpro
  
  tmt <- switch(tmt,
                "tmt" = "n_term_229_1629",
                "tmtpro" = "n_term_304_207",
                "tmtpro0" = "n_term_295_1896")
  
  psm4 <- psm3 %>%
    dplyr::mutate(is.labeled = dplyr::case_when(
      stringr::str_detect(Assigned.Modifications, tmt) ~ TRUE,
      .default = FALSE)) %>% 
    dplyr::filter(is.labeled == TRUE)
  
  cat(nrow(psm3) - nrow(psm4), "PSMs were not TMT labeled.\n")
  
  cat("Number of peptides after filtering:",
      length(unique(psm4$Peptide)), "peptides.\n")
  
  cat("Number of proteins after filtering:",
      length(unique(psm4$Protein.ID)), "proteins.\n")
  
  # Calculate grouping charges
  charges <- psm4 %>%
    dplyr::select(Protein.ID,
                  Peptide,
                  Charge,
                  dplyr::any_of(metadata$key)) %>%
    tidyr::pivot_longer(cols = dplyr::any_of(metadata$key),
                        names_to = "key",
                        values_to = "vals") %>%
    dplyr::select(-vals) %>%  
    dplyr::group_by(Protein.ID,
                    Peptide,
                    key) %>% 
    dplyr::summarise(dplyr::across(dplyr::everything(), 
                                   ~ toString(unique(.))),
                     .groups = "drop")
  
  
  # Aggregate PSM information to protein information
  # First do it to the peptide level (non-redundant)
  pep <- psm4 %>%
    dplyr::select(Protein.ID,
                  Peptide,
                  dplyr::any_of(metadata$key)) %>%
    tidyr::pivot_longer(cols = dplyr::any_of(metadata$key),
                        names_to = "key",
                        values_to = "vals") %>%
    dplyr::group_by(key,
                    Protein.ID,
                    Peptide) %>%
    dplyr::summarise(sum_int = sum(vals),
                     .groups = "drop") %>%
    dplyr::left_join(charges,
                     by = c("Protein.ID", "Peptide", "key")) %>%
    tidyr::pivot_wider(names_from = "key",
                       values_from = "sum_int") %>% 
    dplyr::left_join(dplyr::select(psm3,
                                   Protein.ID,
                                   Gene,
                                   Protein.Description,
                                   Peptide,
                                   Peptide.Length,
                                   Protein.Start,
                                   Protein.End),
                     by = c("Protein.ID",
                            "Peptide")) %>% 
    dplyr::relocate(dplyr::any_of(metadata$key),
                    .after = dplyr::last_col()) %>% 
    dplyr::distinct()
  
  # Then, from the peptide level, aggregate (sum) the peptides to proteins
  protein <- pep %>% 
    dplyr::select(Protein.ID,
                  Peptide,
                  dplyr::any_of(metadata$key)) %>%
    tidyr::pivot_longer(cols = dplyr::any_of(metadata$key),
                        names_to = "key",
                        values_to = "vals") %>%
    dplyr::group_by(key,
                    Protein.ID) %>%
    dplyr::summarise(sum_int = sum(vals),
                     .groups = "drop") %>%
    tidyr::pivot_wider(names_from = "key",
                       values_from = "sum_int") %>% 
    dplyr::left_join(dplyr::select(psm4,
                                   Protein.ID,
                                   Gene,
                                   Protein.Description),
                     by = "Protein.ID") %>% 
    dplyr::relocate(dplyr::any_of(metadata$key),
                    .after = dplyr::last_col()) %>% 
    dplyr::distinct() %>% 
    dplyr::rename(Index = Protein.ID)
  
  listy <- list(peptide = dplyr::rename(pep, Index = Peptide),
                protein = protein)
  
  return(listy)
  
}