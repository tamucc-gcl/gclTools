load_files_for_quant <- function (raw_data, plate_map) {

  # Step 1. Validate input files
  if (!file.exists(raw_data)) {
    stop ("The specified raw data file does not exist.")
  }

  if (!file.exists(plate_map)) {
    stop ("The specified plate map file does not exist.")
  }

  # Step 2 Load required libraries
  required_packages <- c("tidyverse", "janitor")

  ## Helper funtion to ensure required packages are installed
  check_and_load_packages <- function(packages) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }

  check_and_load_packages(required_packages)

  # Step 3: Read input files

  ## Check if the provided file are csv files.
  if (!grepl("\\.csv$", raw_data, ignore.case = TRUE)) {
    stop("The `raw_data` file must be a CSV file.")
  }
  if (!grepl("\\.csv$", plate_map, ignore.case = TRUE)) {
    stop("The `plate_map` file must be a CSV file.")
  }

  raw_data_df <- read.csv (raw_data) %>% clean_names()

  plate_map_df <- read.csv (plate_map) %>% clean_names()

  # Step 4: Data validation
  required_raw_cols <- c("sample", "wells", "value", "r")
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

  # Step 5: Report the number of standards used and their concentrations
  standard_count <- plate_map_df %>% filter (plate_id == "standard") %>% nrow()
  message (paste(standard_count, "standards detected"))

  standard_samples <- plate_map_df %>% filter (plate_id == "standard") %>% pull(sample_id)

  ## Report the standards used
  message(paste("The standards used (in ng/ul) are:", paste(standard_samples, collapse = ", ")))


  # Step 6: Merge data frames and rename "value" to "rfu"
  quant_data <-
    left_join(plate_map_df, raw_data_df, by = "wells") %>%
    rename (rfu = value) %>%
    mutate (category = case_when (plate_id == "standard" ~ "standard",
                                  T ~ "sample")) %>%
    split (.$category) %>%
    map (., ~select (.x, -category))


}
