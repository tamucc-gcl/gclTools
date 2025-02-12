#### format_plate_map.R ####
## Kevin Labrador
## 2025-02-04

# Set-up the working directory in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Clear global environment.
rm(list = ls())

#### USER DEFINED VARIABLES ####
## Input file path
infile <- "./rbd_edna-plate05_AccuClear_2024-08-27.xlsx"

# List sheets in the excel file
sheets <- excel_sheets(infile)

# Use regex to match sheet names
sheet_rawdata <- grep("(?i)^rawdata$", sheets, value = TRUE)
sheet_datacalcs <- grep("(?i)^datacalcs$", sheets, value = TRUE)

## Output file path
outfile <- 
  basename(infile) %>% 
  gsub (".xlsx", "_plate-map_formatted.csv", .)

#### PROGRAM DEFINED VARIABLES

#### Load Libraries ####
pacman::p_load(
  janitor,
  readxl,
  tidyverse
)

#### Prepare plate map ####
plate_map <- 
  read_xlsx(
    infile,
    sheet = sheet_rawdata
  ) %>% 
  clean_names() %>% 
  
  # Change "Standards" to "standard"
  mutate (
    plate_id =
      case_when (
        grepl ("[Ss]tandards?", sample_plate) ~ "standard",
        T ~ sample_plate
      )
  ) %>% 
  
  # Add leading zero to plate number (if it's a single digit)
  mutate(
    plate_id = sub("(rbd-eDNA-Plate)(\\d)$", "\\10\\2", plate_id)
  ) %>%
  
  # Add sample_id 
  mutate (
    sample_id = paste0(
      sample_row, 
      sample_column)
  ) %>% 
  
  # Select only the necessary files
  select (
    plate_id,
    sample_id,
    sample_row,
    sample_column,
    quant_row,
    quant_column
  ) 


#### Extract replicate information ####
data_calcs <- 
  read_xlsx(
    infile,
    sheet = sheet_datacalcs
  ) %>% 
  clean_names() 

replicate_information <-  
  data_calcs %>%  
  # Select only the necessary rows with replicate information
  select (
    sample_plate,
    sample_column, 
    sample_row,
    quant_row_1,
    quant_row_2,
    quant_col_1,
    quant_col_2) %>% 
  
  # Pivot longer
  pivot_longer(
    cols = starts_with ("quant_"),
    names_to = c(".value", "replicate"),
    names_pattern = "quant_(row|col)_(\\d+)") %>% 
  
  # Rename columns
  rename (
    quant_row = row,
    quant_column = col,
    plate_id = sample_plate
  ) %>% 
  
  # Coerce replicate into numeric
  mutate (
    replicate = as.numeric(replicate)
  )%>% 
  
  # Add leading zero to plate number (if it's a single digit)
  mutate(
    plate_id = sub("(rbd-eDNA-Plate)(\\d)$", "\\10\\2", plate_id)
  )

#### Check if there are samples that need to be dropped ####
samples_to_remove <- 
  data_calcs %>%   
  mutate(
    across(
      matches("_(1|2)$"), 
      as.character)
  ) %>%  # Convert numeric columns to character
  mutate(
    sample_id = 
      paste0(sample_row, sample_column)
  ) %>%  
  pivot_longer(
    cols = matches("_(1|2)$"),
    names_to = c(".value", "replicate"), 
    names_pattern = "(.*)_(\\d)$"
  ) %>% 
  mutate (raw_rfu = as.numeric (raw_rfu)) %>% 
  filter (is.na(raw_rfu))

  
  #### Extract sample volume information ####
sample_volume <- 
  data_calcs %>% 
  select (c(
    "ng_well_1",
    "ng_ul_1")
  ) %>% 
  mutate (
    sample_volume = ng_well_1 / ng_ul_1
  ) %>% 
  pull (sample_volume) %>% 
  unique()

#### Extract sample information ####
plate_map_sample <- 
  plate_map %>% 
  filter (plate_id != "standard") %>% 
  mutate (
    sample_column = as.numeric(sample_column),
    sample_volume = sample_volume) %>% 
  filter(
    ! paste0(quant_row, quant_column) %in% paste0(samples_to_remove$quant_row, samples_to_remove$quant_col)
    ) %>% 
  left_join (replicate_information) 

#### Extract standard information ####
plate_map_standard <- 
  plate_map %>% 
  filter (plate_id == "standard") %>% 
  mutate (sample_id = as.character(sample_column),
          sample_volume = mean( (as.numeric(sample_row) / sample_column), na.rm = T ),
          replicate = 1
  ) %>% 
  select (-c(
    sample_column, 
    sample_row)
  )

#### Merge sample and standard plate maps ####
plate_map_formatted <- 
  full_join(
    plate_map_sample,
    plate_map_standard) 


# Export file
write_csv(
  plate_map_formatted,
  file = outfile)



