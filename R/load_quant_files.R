#' Analyze quant fluorescence data
#'
#' @description
#' Prepare the raw data files for fluorescence analysis. This function reads the raw data and plate map, validates the file format and data structure, merges the files, and then prepares a list that separates the samples from the standards.
#'
#' This is Step 1 of 4 in a series of functions for analyzing fluorescent data.
#'
#' @import tidyverse
#' @import janitor
#'
#' @param raw_data
#' path to the raw plate reader data saved as a comma-separated value (csv) file. This is the file exported from the SpectraMax software. It contains the well and fluorescence information. This file should not be manipulated once it has been exported.
#'
#' @param plate_map
#' path to the assay plate map. This contains the information on how the samples were laid out in the 364-well plate and its associated metadata. This file follows a specific data structure. The following columns are expected:
#'
#' - plate_id: character string. Name of the sample plate/strip assayed. All standards used in the assay should be named "standard".
#'
#'
#' - sample_id: character string. Name of the samples in the plate. This usually corresponds to well position. Standards used in this assay should be named based on the concentration (ng/ul) used.
#'
#'
#' - replicate: numeric. The replicate used in the assay. A value of 2 is expected for "double-quants"
#'
#'
#' - quant_row: character. Typically a letter indicating the row of the sample in the quant plate.
#'
#'
#' - quant_column: numeric. Typically a number indicating the column of the sample in the quant plate.
#'
#'
#' - sample_volume: numeric. Indicates the volume (uL) of samples and standards used in the assay.
#'
#'
#' @returns
#' Returns a list containing two data frames: `sample` and `standard`
#'
#' @examples
#' require (tamuccGCL)
#'
#' # Import data files
#' raw_data <- system.file("extdata", "raw_data.csv", package = "tamuccGCL")
#' plate_map <- system.file("extdata", "plate_map.csv", package = "tamuccGCL")
#'
#' quant_data <- load_quant_files(raw_data, plate_map)
#'
#' @export

load_quant_files <- function (raw_data, plate_map) {

  # Step 1. Validate input files
  if (!file.exists(raw_data)) {
    stop ("The specified raw data file does not exist.")
  }

  if (!file.exists(plate_map)) {
    stop ("The specified plate map file does not exist.")
  }

  # Step 2: Read input files

  ## Check if the provided file are csv files.
  if (!grepl("\\.csv$", raw_data, ignore.case = TRUE)) {
    stop("The `raw_data` file must be a CSV file.")
  }
  if (!grepl("\\.csv$", plate_map, ignore.case = TRUE)) {
    stop("The `plate_map` file must be a CSV file.")
  }

  raw_data_df <- read.csv (raw_data) %>% clean_names() %>% select (wells, value)

  plate_map_df <- read.csv (plate_map) %>% clean_names()

  # Step 3: Data validation
  required_raw_cols <- c("wells", "value")
  required_plate_cols <- c("plate_id", "sample_id", "replicate", "quant_row", "quant_column", "sample_volume")

  ## Validate and modify the raw data file

  missing_raw_cols <- setdiff(required_raw_cols, names(raw_data_df))
  if (length(missing_raw_cols) > 0) {
    stop("The raw data file is missing the following required columns: ", paste(missing_raw_cols, collapse = ", "))
  } else {
    raw_data_df <- raw_data_df %>% select(wells, value)
  }

  ## Validate and modify the plate map file
  missing_plate_cols <- setdiff(required_plate_cols, names(plate_map_df))
  if (length(missing_plate_cols) > 0) {
    stop("The plate map file is missing the following required columns: ", paste(missing_plate_cols, collapse = ", "))
  } else {
    # Check if "standards" are present in the plate_id column
    if (!"standard" %in% plate_map_df$plate_id) {
      stop ("No standards detected.\n  Standards used in the assay should be named as 'standard' under the plate_id column in the plate map. Morevoer, the standard concentration (ng/ul) should be provided in the sample_id column.")
    } else {
      plate_map_df <- plate_map_df %>%
        mutate(wells = paste0(quant_row, quant_column))
    }
  }

  # Step 4: Report the number of standards used and their concentrations
  standard_count <- plate_map_df %>% filter (plate_id == "standard") %>% nrow()
  message (paste(standard_count, "standards detected"))

  standard_samples <- plate_map_df %>% filter (plate_id == "standard") %>% pull(sample_id)

  ## Report the standards used
  message(paste("The standards used (in ng/ul) are:", paste(standard_samples, collapse = ", ")))


  # Step 5: Merge data frames and rename "value" to "rfu"
  quant_data <-
    left_join(plate_map_df, raw_data_df, by = "wells") %>%
    rename (rfu = value) %>%
    mutate (category = case_when (plate_id == "standard" ~ "standard",
                                  T ~ "sample")) %>%
    split (.$category) %>%
    map (., ~select (.x, -category))


}
