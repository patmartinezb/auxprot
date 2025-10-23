#' Import TMT proteomic data
#'
#' @description `import()` imports one or several files of TMT proteomic data,
#'   from several formats (txt, csv, tsv). If several files are indicated, and
#'   depending on the complexity level (peptide spectrum matches (PSM) or
#'   protein) it can collapse the data to the desired level (either peptide or
#'   protein). It's also TMT-sensitive (TMT, TMT-pro or TMT-0) if input data is
#'   at the PSM-level.
#'
#' @param input_file_name String vector with name(s) of the files, including
#'   their extension
#' @param metadata Dataframe with associated metadata
#' @param type Either "psm" or "protein", indicating the level of the input data
#'   - only relevant if the input comprises several files, if only one file, then NULL
#' @param org Organism - only relevant if the input is several files, if only
#'   one file, then NULL Input vector of length 1.
#' @param tmt Indicates the type of TMT used, either `tmt`, `tmtpro` or
#'   `tmtpro0`.
#' @param use.unique Whether to use only unique peptides for quantification.
#'   Default is FALSE.
#' @param level Either "peptide" or "protein", indicating the desired
#'   output-level data - only relevant if the input comprises several files, if
#'   only one file, then NULL
#'
#' @returns A dataframe at the complexity level specified by the user
#' @export import
import <- function(input_file_name,
                   metadata = NULL,
                   type = NULL,
                   org = NULL,
                   tmt = NULL,
                   use.unique = FALSE,
                   level = NULL){
  
  # Main import function
  
  if (length(input_file_name) == 1){
    
    # Imports all data, checking the extension
    
    ext <- tools::file_ext(input_file_name)
    switch(ext,
           csv = vroom::vroom(input_file_name, delim = ","),
           txt = vroom::vroom(input_file_name, delim = "\t"),
           tsv = vroom::vroom(input_file_name, delim = "\t"),
           validate("Invalid file; Please upload a .csv, .tsv or .txt file")
    )
    
  } else {
    
    # Imports several files - .txt or .tsv, as yielded by FragPipe
    
    ext <- import_several(input_file_name, metadata, type, org, tmt, use.unique, level)
    
  }
}


#' Clean metadata
#'
#' @description `clean_key()` prepares the metadata so that all names match
#'   names in the data matrix, and clean the names so that there are no
#'   conflicts.
#'
#' @param metadata Dataframe with metadata
#'
#' @returns A cleaned dataframe
#' @export clean_key
clean_key <- function(metadata){
  
  # Makes sure that the metadata mandatory variables are written correctly and/or are all present
  
  true_names <- c("Run", "Mixture", "TechRepMixture", "Fraction", "Channel",
                  "Condition", "BioReplicate", "key")
  
  if (length(true_names %in% colnames(metadata)) != length(true_names)){
    
    stop("Either some mandatory metadata variables are missing, or they are spelled wrong.
         Check the 'Help' tab.")
    
  }
  
  metadata$key <- make.names(metadata$key)
  
  two_mand <- c("Mixture", "Condition")
  
  vars <- c(two_mand, colnames(metadata)[-which(colnames(metadata) %in% true_names)])
  
  metadata[, vars] <- apply(metadata[, vars], 2, function(x){
    stringr::str_to_title(tolower(x))
  })
  
  return(metadata)
  
}


#' Clean comparison matrix
#'
#' @description `clean_comp()` prepares the comparison file so that all names
#'   match names in the metadata, and clean the names so that there are no
#'   conflicts.
#'
#' @param comp Dataframe with comparisons
#'
#' @returns A cleaned dataframe
#' @export clean_comp
clean_comp <- function(comp){
  
  # Makes sure that the comparison mandatory variables are written correctly and/or are all present
  
  true_names <- c("Condition", "Control")
  
  if (length(true_names %in% colnames(comp)) != length(true_names)){
    
    stop("Either some mandatory comparison variables are missing, or they are spelled wrong.
         Check the 'Help' tab.")
    
  }
  
  # So that names are easier to parse - The first word has to be in capital letters
  
  if (nrow(comp) > 1){
    
    comp <- apply(comp, 2, function(x){
      stringr::str_to_title(tolower(x))
    }) |> tibble::as_tibble()
    
  } else {
    
    comp <- apply(comp, 2, function(x){
      stringr::str_to_title(tolower(x))
    })
    
    comp <- tibble::as_tibble(data.frame(as.list(comp)))
  }
  
  return(comp)
  
}


extract_protein_tsv <- function(input_file_paths,
                                names,
                                type = c("psm", "protein")){
  
  # Reads file and selects variables of interest
  
  df <- data.table::fread(input_file_paths, 
                          sep = "\t", 
                          header = T, 
                          na.strings = "NA") %>%
    dplyr::rename_with(.,
                       ~ make.names(.),
                       dplyr::everything()) %>%
    dplyr::filter(!grepl("contam", Protein)) %>% # Removes contaminants
    type.convert(as.is = TRUE)
  
  if (all(sapply(df[, which(colnames(df) %in% names)], is.numeric))){
    
    df <- df
    
  } else { # In case there is any problem with the locale and type.convert() does not work
    
    if (any(grepl("\\.", df[, which(colnames(df) %in% names)]))){
      
      df[,-1] <- apply(df[,-1], 2, function(x){
        
        temp <- gsub("\\.", "", x)
        
        as.numeric(temp)
        
      })
      
    } else { # If it doesn't coerce to numeric, then I have no clue
      
      stop("The intensity values cannot be read as numeric values. Please check format.")
      
    }
    
  }
  
  if (type != "psm"){ # if df at protein level, only the following info is necessary
    
    df <- df %>%
      dplyr::select(Protein.ID,
                    tidyselect::any_of(names))
    
  }
  
  return(df)
  
}


import_several <- function(input_file_paths,
                           metadata,
                           type = c("psm", "protein"),
                           org,
                           tmt,
                           use.unique,
                           level = c("peptide", "protein")){
  
  # Recursively imports files and outputs a single dataset
  
  if ("key" %in% colnames(metadata)){
    
    if (type == "psm"){
      
      ext0 <- input_file_paths |>
        purrr::set_names() |>
        purrr::map(extract_protein_tsv, names = metadata$key, type = type) |>
        purrr::map(tmt_integrator, metadata, org, tmt, use.unique)
      
      if (level == "peptide"){
        
        ext <- ext0 |>
          purrr::modify_depth( 1, level) |>
          purrr::map(dplyr::select, !c(Charge, n_psm, Is.Unique)) |>
          purrr::reduce(dplyr::full_join, by = c("Protein.ID",
                                                 "Index",
                                                 "Gene",
                                                 "Protein.Description",
                                                 "Peptide.Length",
                                                 "Protein.Start",
                                                 "Protein.End"))
        
      } else {
        
        ext <- ext0 |>
          purrr::modify_depth(1, level) |>
          purrr::map(dplyr::select, !c(n_peptides:n_psm)) |>
          purrr::reduce(dplyr::full_join, by = c("Index",
                                                 "Gene",
                                                 "Protein.Description"))
        
      }
      
    } else {
      
      ext <- input_file_paths |>
        purrr::map(extract_protein_tsv, names = metadata$key, type = type) |>
        purrr::reduce(dplyr::full_join, by = "Protein.ID")
      
    }
    
    ext <- ext %>%
      dplyr::rename_with(.,
                         ~ make.names(.),
                         dplyr::everything())
    
    return(ext)
    
  } else {
    
    stop("There is no 'key' variable in the metadata file")
    
  }
}


select_mq_vars <- function(raw, metadata, reporter_names, phospho = "no"){
  # raw: raw data
  # metadata
  # reporter names: key names
  # phospho: whether it's phosphoproteomic data or not, default is no
  
  # Selects variables of interests from MQ output
  
  raw <- raw %>%
    dplyr::rename_with(.,
                       ~ make.names(.),
                       dplyr::everything()) # clean colnames so that they match the metadata
  
  if (phospho == "yes"){
    # Phospho
    
    raw.names <- colnames(raw)[grep(paste(reporter_names, collapse = "|"), colnames(raw))]
    foo <- unique(stringr::str_sub(raw.names, end = -5))
    
    if (all(reporter_names %in% foo)){
      
      df <- raw %>%
        dplyr::select(
          dplyr::starts_with(reporter_names),
          Protein,
          Gene.names,
          Localization.prob,
          Positions.within.proteins,
          Amino.acid,
          Sequence.window,
          Protein.names,
          Reverse,
          Potential.contaminant)
      
    } else {
      
      stop("Names in metadata 'key' do not match names in data matrix")
      
    }
    
  } else {
    
    if (all(reporter_names %in% colnames(raw))){
      
      df <- raw %>%
        dplyr::select(
          Protein.IDs,
          Q.value,
          tidyselect::matches(reporter_names),
          Only.identified.by.site,
          Reverse,
          Potential.contaminant,
          Protein.names,
          Gene.names,
          Fasta.headers)
      
    } else {
      
      stop("Names in metadata 'key' do not match names in data matrix")
      
    }
  }
  return(df)
}


clean_mq_vars <- function(raw, reporter_names, phospho = "no"){
  # raw: data matrix
  # reporter_names: key names
  # phospho: whether it's phosphoproteomic data or not, default is no
  
  if (phospho == "yes"){
    
    rev.prots <- raw$Protein[grep("*REV_", raw$Protein)]
    
    df <- raw %>%
      dplyr::filter(!Reverse %in% "+",
                    !Potential.contaminant %in% "+", # Remove contaminants
                    Localization.prob >= .75, # Keep only high quality hits - localization probability over 0.75
                    !Protein %in% rev.prots) %>% # In case some decoy is not cleaned
      dplyr::select(-Reverse,
                    -Potential.contaminant,
                    -Localization.prob) %>%
      dplyr::mutate(Protein = sub(";.*", "", Protein), # clean variables
                    Gene.names = sub(";.*", "", Gene.names),
                    Protein.names = sub(";.*", "", Protein.names),
                    Positions.within.proteins = sub(";.*", "", Positions.within.proteins),
                    Sequence.window = sub(";.*", "", Sequence.window),
                    Phosphosite = paste0(Amino.acid, Positions.within.proteins),
                    PPS = paste(Protein, Phosphosite, sep = "-")) %>% # Complete phosphosite name
      tidyr::pivot_longer(!Protein:PPS, names_to = "samples", values_to = "vals") %>%
      dplyr::mutate(multiplicity = stringr::str_sub(samples, start = -4), # Get the multiplicity for each phosphopeptide
                    biorep = stringr::str_sub(samples, end = -5)) %>% # Clean bioreplicate names
      dplyr::select(- samples) %>%
      tidyr::pivot_wider(names_from = "biorep", values_from = "vals") %>%
      dplyr::mutate(Protein.IDs = paste0(PPS, multiplicity)) # Create new PPS id with multiplicity - named Protein.IDs so that the code doesn't break downstream
    
  } else {
    
    df <- raw %>%
      dplyr::select(-Q.value) %>%
      dplyr::filter(!Reverse %in% "+",
                    !Only.identified.by.site %in% "+",
                    !Potential.contaminant %in% "+") %>% # Remove contaminants
      dplyr::select(-Reverse,
                    -Only.identified.by.site,
                    -Potential.contaminant) %>%
      dplyr::mutate(Protein.IDs = sub(";.*", "", Protein.IDs))
  }
  
  # General handling
  df <- df %>%
    dplyr::mutate(
      dplyr::across(tidyselect::all_of({{ reporter_names }}), ~ dplyr::na_if(.x, 0)),
      dplyr::across(tidyselect::all_of({{ reporter_names }}), ~ dplyr::na_if(.x, Inf)),
      dplyr::across(tidyselect::all_of({{ reporter_names }}), ~ dplyr::na_if(.x, -Inf))
    ) # Transform to NAs
  
  return(df)
  
}


process_mq <- function(raw, metadata, phospho = "no"){
  
  # Main MQ processing
  
  # Clean reporter names
  reporter_names_clean <- make.names(metadata$key)
  
  # Select variables
  df <- select_mq_vars(raw, metadata, reporter_names = reporter_names_clean, phospho)
  
  # Store those proteins with unreliable info, to check later if they were removed
  # names_q <- df$Protein.IDs[which(df$Q.value > 0.01)]
  
  # Clean variables
  df <- clean_mq_vars(df, reporter_names_clean, phospho)
  
  
  # Checking Q values
  # names_q_after <- df$Protein.IDs[df$Protein.IDs %in% names_q]
  #
  # if (length(names_q_after) != 0){
  #
  #   df <- df[-which(df$Protein.IDs %in% names_q_after),]
  # }
  
  # Check / remove duplicates
  # if (any(duplicated(na.omit(df$Gene.names)))) {
  #
  #   df <- df %>%
  #     dplyr::distinct(Gene.names, .keep_all = TRUE)
  #
  # }
  
  
  if (any(duplicated(na.omit(df$Protein.IDs)))) {
    
    df <- df %>%
      dplyr::distinct(Protein.IDs, .keep_all = TRUE)
    
  }
  
  return(df)
  
}



process_fp <- function(raw, metadata){
  
  # Main FragPipe processing
  
  # Clean raw data colnames
  raw <- raw %>%
    dplyr::rename_with(.,
                       ~ make.names(.),
                       dplyr::everything())
  
  # Clean reporter names
  reporter_names_clean <- make.names(metadata$key)
  
  # Make sure that the sample ids match with those in the metadata
  if (all(reporter_names_clean %in% colnames(raw))){
    
    if ("Index" %in% colnames(raw)){
      
      df <- raw %>%
        dplyr::select(Index,
                      {{ reporter_names_clean }}) %>%
        dplyr::mutate(
          dplyr::across(tidyselect::where(is.numeric), ~ dplyr::na_if(.x, 0)),
          dplyr::across(tidyselect::where(is.numeric), ~ dplyr::na_if(.x, Inf)),
          dplyr::across(tidyselect::where(is.numeric), ~ dplyr::na_if(.x, -Inf))
        ) %>%  # Transform to NAs
        dplyr::rename(Protein.IDs = Index)
      
    } else {
      
      df <- raw %>%
        dplyr::mutate(
          dplyr::across(tidyselect::where(is.numeric), ~ dplyr::na_if(.x, 0)),
          dplyr::across(tidyselect::where(is.numeric), ~ dplyr::na_if(.x, Inf)),
          dplyr::across(tidyselect::where(is.numeric), ~ dplyr::na_if(.x, -Inf))
        ) %>%  # Transform to NAs
        dplyr::rename(Protein.IDs = Protein.ID) %>%
        dplyr::select(Protein.IDs,
                      {{ reporter_names_clean }})
    }
    
    if (all(sapply(df[, reporter_names_clean], function(x) class(x) %in% c("integer","numeric", "double")))){
      
      df <- df
      
    } else { # In case there is any problem with the locale and type.convert() does not work
      
      if (any(grepl("\\.", df[, reporter_names_clean]))){
        
        df[,-1] <- apply(df[,-1], 2, function(x){
          
          temp <- gsub("\\.", "", x)
          
          as.numeric(temp)
          
        })
        
      } else { # If it doesn't coerce to numeric, then I have no clue
        
        stop("The intensity values cannot be read as numeric values. Please check format.")
        
      }
      
    }
    
    return(df)
    
  } else {
    
    stop("Names in metadata 'key' do not match names in data matrix")
    
  }
}

process_other <- function(raw, metadata){
  
  # "Other" processing
  
  # Clean raw data colnames
  raw <- raw %>%
    dplyr::rename_with(.,
                       ~ make.names(.),
                       dplyr::everything())
  
  # Clean reporter names
  reporter_names_clean <- make.names(metadata$key)
  
  # Make sure that the sample ids match with those in the metadata
  if (all(reporter_names_clean %in% colnames(raw))){
    
    if (any(c("Index", "Accession", "Protein.IDs") %in% colnames(raw))){
      
      df <- raw %>%
        dplyr::select(dplyr::starts_with("Index"),
                      dplyr::starts_with("Accession"),
                      dplyr::starts_with("Protein.IDs"),
                      {{ reporter_names_clean }}) %>%
        dplyr::mutate(
          dplyr::across(tidyselect::where(is.numeric), ~ dplyr::na_if(.x, 0)),
          dplyr::across(tidyselect::where(is.numeric), ~ dplyr::na_if(.x, Inf)),
          dplyr::across(tidyselect::where(is.numeric), ~ dplyr::na_if(.x, -Inf))
        )  # Transform to NAs
      
      colnames(df) <- gsub("Index|Accession", "Protein.IDs", colnames(df))
      
      return(df)
      
    } else {
      
      stop("There is no 'Index' or 'Accession' variable to mark Protein IDs")
    }
    
  } else {
    
    stop("Names in metadata 'key' do not match names in data matrix")
    
  }
}


#' Clean proteomic data matrix
#'
#' @description `process_raw()` prepares the proteomic data matrix, depending on
#'   its origin (FragPipe, MaxQuant, or others), so that it can be work with
#'   safely. It cleans variable names, transforms 0 into NA, and removes
#'   contaminants.
#'
#' @param software String of length 1 indicating the software of origin, either
#'   FragPipe ("fp"), MaxQuant ("mq") or others ("other").
#' @param raw Dataframe with the proteomic data. If "other" was selected in
#'   `software`, the matrix should contain only the numeric matrix where the
#'   column names correspond to what's indicated in the metadata; and a column
#'   with the protein identifications (called either "Index", "Accession" or
#'   "Protein.IDs").
#' @param metadata Dataframe with the corresponding metadata.
#' @param phospho TRUE or FALSE (default), whether `raw` is a phosphoproteomic
#'   data matrix.
#'
#' @returns A dataframe with the proteomic data ready to work with.
#' @export process_raw
process_raw <- function(software, raw, metadata, phospho = FALSE){
  
  # Depending on the software of origin, use the corresponding function for processing
  
  switch(software,
         "mq" = process_mq(raw, metadata, phospho), # MaxQuant
         "fp" = process_fp(raw, metadata), # FragPipe
         "other" = process_other(raw, metadata)
  )
  
}