#' FragPipe's TMT Integrator wrapper
#'
#' @description `tmt_integrator()` mimics the FragPipe's TMT Integrator module
#'   filters for TMT data, minus de PSM normalization. This includes purity and
#'   peptide probability filters, as well as a minimum summed MS2 intensity,
#'   presence of TMT and, optionally, PTM labeling. If an Astral mass spec is
#'   used, it also performs minimum summed SNR and resolution filtering (but not
#'   minimum MS2 intensity). It summarises PSM information to the peptide and
#'   subsequent protein level, for further analysis.
#'
#'
#' @param psm Input FragPipe data matrix at the PSM level.
#' @param metadata Input dataframe. Metadata associated to `psm`.
#' @param org Input vector of length >= 1. Indicates the species (could be
#'   several) the data belongs to.
#' @param tmt Input vector of length 1. Indicates the type of TMT used, either
#'   `tmt`, `tmtpro` or `tmtpro0`.
#' @param astral Whether the mass spec used was an Astral (default FALSE)
#' @param use.unique Whether to use only unique peptides for quantification.
#'   Default is FALSE.
#' @param purity_thr Numeric vector of length 1. Indicates purity threshold.
#'   Default is 0.5 (50%).
#' @param pep_prob_thr Numeric vector of length 1. Indicates peptide probability
#'   threshold. Default is 0.9.
#' @param sum_ms2_int_thr Numeric vector of length 1. Indicates the percentage
#'   of MS2 intensity below which PSMs are going to be filtered out. Default is
#'   0.05 (5%). For PTM data, such as phosphorylation, it should be 0.025
#'   (2.5%).
#' @param min_resolution Numeric vector of length 1. Indicates the minimum
#'   resolution below which PSMs are going to be filtered out. Default is
#'   45,000.
#' @param min_snr Numeric vector of length 1. Indicates the minimum
#'   signal-to-noise ratio (SNR) below which PSMs (sum of all channels) are
#'   going to be filtered out. Default is 1,000.
#'
#' @return A list containing both peptide and protein level information.
#' @export tmt_integrator
tmt_integrator <- function(psm, 
                           metadata, 
                           org, 
                           tmt,
                           astral = FALSE,
                           use.unique = FALSE,
                           purity_thr = 0.5,
                           pep_prob_thr = 0.9,
                           sum_ms2_int_thr = 0.05,
                           min_resolution = 45000,
                           min_snr = 1000){
  
  # TODO: if phosphoproteomics, sum_ms2_int_thr is 0.025 (2.5%), and must contain phosphorylation
  # NOTE: The janitor::make_clean_names(Assigned.Modifications) takes forever when the dataset is large
  
  # 0. Initial transformation and filtering for contaminants of PSM data
  psm <- psm %>% 
    tibble::as_tibble(.name_repair = make.names) %>% # make.names() preserves the sample names as in metadata
    dplyr::rename_with(., # name fix for fragpipe v23
                       ~ gsub("Intensity.", "", .x),
                       dplyr::starts_with("Intensity.")) %>%
    dplyr::mutate(sum_int = rowSums(dplyr::across(dplyr::any_of(metadata$key)), 
                                    na.rm = TRUE)) %>%  # sum all intensitites, rowwise
    dplyr::mutate(Assigned.Modifications = tolower(Assigned.Modifications),
                  Assigned.Modifications = gsub("\\-|\\(|\\)", 
                                                ".", 
                                                Assigned.Modifications),
                  Assigned.Modifications = gsub("\\.", "_", Assigned.Modifications),
                  organism = sub('.*\\_', '', Protein)) %>%  
    dplyr::filter(organism %in% toupper(org)) # remove contaminants from other species
  
  if (isTRUE(astral)){
    
    psm <- psm %>% 
      dplyr::filter(Is.Contaminant == FALSE,
                    Is.Decoy == FALSE)
  }
  
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
  
  message("There are ", nrow(psm[psm$Purity <= purity_thr,]), " PSMs below the ",
          purity_thr, " purity threshold.\n")
  
  psm1 <- psm %>%
    dplyr::filter(Purity >= purity_thr)
  
  # 3. Filtering by minimum peptide probability
  
  message("There are ", nrow(psm1[psm1$Probability <= pep_prob_thr,]), 
          " PSMs below the ", pep_prob_thr, " probability threshold.\n")
  
  psm2 <- psm1 %>%
    dplyr::filter(Probability >= pep_prob_thr)
  
  # 4.1 Summed MS2 intensity filtering OR filtering by minimum resolution (e.g., 45,000)
  
  if (isTRUE(astral)){
    
    resol_names <- paste0("Resolution.", metadata$key)
    
    resol_out <- psm2 %>% 
      dplyr::select(Spectrum,
                    dplyr::any_of(resol_names)) %>% 
      tidyr::pivot_longer(cols = dplyr::any_of(resol_names),
                          names_to = "names",
                          values_to = "resolution") %>% 
      dplyr::filter(resolution < min_resolution) %>% 
      dplyr::pull(Spectrum) %>% 
      unique
    
    message("There are ", length(resol_out), " PSMs below the ",
            min_resolution, " minimum resolution threshold.\n")
    
    psm3 <- psm2 %>%
      dplyr::filter(!Spectrum %in% resol_out)
    
  } else {
    
    five.pc <- quantile(psm$sum_int, probs = sum_ms2_int_thr)
    
    # stringr::str_glue("The summed MS2 intensity {sum_ms2_int_thr*100}% threshold is {five.pc}.")
    
    message("There are ", nrow(psm2[psm2$sum_int < five.pc,]), " PSMs below the ",
            sum_ms2_int_thr*100, " % summed MS2 intensity threshold.\n")
    
    psm3 <- psm2 %>%
      dplyr::filter(sum_int >= five.pc)
    
  }
  
  # 4.2 if ASTRAL: Requires minimum SNR of 1,000
  # Removes the PSM if all channels' summed SNR is less than the min SNR
  
  if (isTRUE(astral)){
    
    snr_names <- paste0("SNR.", metadata$key)
    
    snr_out <- psm3 %>% 
      dplyr::select(Spectrum,
                    dplyr::any_of(snr_names)) %>% 
      tidyr::pivot_longer(cols = dplyr::any_of(snr_names),
                          names_to = "names",
                          values_to = "snr") %>% 
      dplyr::group_by(Spectrum) %>% 
      dplyr::summarise(sum_snr = sum(snr, na.rm = TRUE)) %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(sum_snr < min_snr) %>% 
      dplyr::pull(Spectrum) %>% 
      unique
    
    
    message("There are ", length(snr_out), " PSMs below the ",
            min_snr, " minimum SNR.\n")
    
    psm3 <- psm3 %>%
      dplyr::filter(!Spectrum %in% snr_out)
    
  }
  
  
  # 5. Must be TMT labeled
  # the weight is different if tmt or tmtpro
  
  tmt <- switch(tmt,
                "tmt" = "229_1629",
                "tmtpro" = "304_207",
                "tmtpro0" = "295_1896")
  
  psm4 <- psm3 %>%
    dplyr::mutate(is.labeled = dplyr::case_when(
      stringr::str_detect(Assigned.Modifications, tmt) ~ TRUE,
      .default = FALSE)) %>% 
    dplyr::filter(is.labeled == TRUE)
  
  message(nrow(psm3) - nrow(psm4), " PSMs were not TMT labeled.\n")
  
  # 6. Use only unique peptides (optional)
  
  if (use.unique == TRUE){
    
    message(nrow(psm4) - nrow(filter(psm4, Is.Unique == TRUE)),
            "PSMs are not unique. \n")
    
    psm4 <- psm4 %>% 
      filter(Is.Unique == TRUE)
    
  }
  
  # Final numbers
  message("Number of peptides after filtering: ",
          length(unique(psm4$Peptide)), " peptides.\n")
  
  message("Number of proteins after filtering: ",
          length(unique(psm4$Protein.ID)), " proteins.\n")
  
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
  
  # Calculate number of total PSMs per peptide and protein
  n_psm <- psm4 %>% 
    dplyr::group_by(Peptide, Protein.ID) %>% 
    dplyr::count(name = "n_psm")
  
  
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
                                   Protein.End,
                                   Is.Unique),
                     by = c("Protein.ID",
                            "Peptide")) %>% 
    dplyr::left_join(n_psm, 
                     by = c("Protein.ID",
                            "Peptide")) %>% 
    dplyr::relocate(dplyr::any_of(metadata$key),
                    .after = dplyr::last_col()) %>%
    dplyr::relocate(c(Charge, n_psm), 
                    .after = Protein.Description) %>% 
    dplyr::distinct()
  
  # Calculate total number of PSMs (including redundant) and number of peptides (unique PSMs)
  pep_n <- pep %>% 
    dplyr::group_by(Protein.ID) %>% 
    dplyr::count(name = "n_peptides") %>% 
    dplyr::ungroup()
  
  n_psm_total <- pep %>% 
    dplyr::group_by(Protein.ID) %>% 
    dplyr::summarise(n_psm = sum(n_psm), 
                     .groups = "drop")
  
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
    dplyr::left_join(pep_n, 
                     by = "Protein.ID") %>% 
    dplyr::left_join(n_psm_total, 
                     by = "Protein.ID") %>% 
    dplyr::relocate(dplyr::any_of(metadata$key),
                    .after = dplyr::last_col()) %>% 
    dplyr::relocate(n_peptides:n_psm, 
                    .after = Protein.Description) %>%
    dplyr::distinct() %>% 
    dplyr::rename(Index = Protein.ID)
  
  listy <- list(peptide = dplyr::rename(pep, Index = Peptide),
                protein = protein)
  
  return(listy)
  
}