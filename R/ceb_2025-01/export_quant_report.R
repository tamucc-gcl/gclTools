#' Analyze quant fluorescence data
#'
#' @description
#' Wrangles and exports the quant data in a format used by the Genomics Core Laboratory (GCL). This takes the object generated via [quant_dna()] as input.
#'
#' This is step 4 of 4 in a series of functions for analyzing fluorescent data.
#'
#' @import tidyverse
#' @import openxlsx
#'
#' @param quant_dna list. Output of [quant_dna()]
#'
#' @returns
#' Generates an excel file of the quant report.
#'
#' @seealso [load_quant_files()], [train_quant_standards()], [quant_dna()]
#'
#' @examples
#' require (tamuccGCL)
#'
#' # Import data files
#' raw_data <- system.file("extdata", "raw_data.csv", package = "tamuccGCL")
#' plate_map <- system.file("extdata", "plate_map.csv", package = "tamuccGCL")
#'
#' quant_data <- load_quant_files(raw_data, plate_map)
#' trained_model <- train_quant_standards (quant_data)
#' quant_report <- quant_dna(quant_data = quant_data, trained_model = trained_model)
#'
#' \dontrun{
#' export_quant_report (quant_report)
#'}
#'
#' @export
export_quant_report <- function (quant_dna) {

  # Step 1: Extract the quant output
  quant_output <- quant_dna$quant_output
  quant_standard <- quant_dna$standard
  quant_plot <- quant_dna$plot

  # Step 2: Wrangle the output_quant
  quant_report <- quant_output %>%

    # Split sample_id into sample column (numeric) and sample_row (character)
    mutate (
      sample_column = as.numeric(str_extract(sample_id, "\\d+")),
      sample_row = str_extract (sample_id, "[A-Za-z]+")
    ) %>%

    # Rename rfu to raw_rfu
    rename (raw_rfu = rfu,
            ng_well = dna_per_well,
            ng_ul = dna_per_ul) %>%

    # Remove unnecessary columns
    select (-c(
      wells,
      sample_volume
      )) %>%

    # Pivot wider based on replicate
    pivot_wider (names_from = replicate,
                 values_from = c(quant_row,
                                 quant_column,
                                 raw_rfu,
                                 zeroed_rfu,
                                 ng_well,
                                 ng_ul)
                 ) %>%

    # Calculate difference, average, and determine pass/redo
    mutate (
      difference = abs(ng_ul_1 - ng_ul_2),
      average = (ng_ul_1 + ng_ul_2) / 2,
      pass_redo = ifelse (difference < (average / 2), "pass", "redo")
    ) %>%

    # reorder columns
    select (
      plate_id,
      sample_row,
      sample_column,
      quant_row_1,
      quant_column_1,
      quant_row_2,
      quant_column_2,
      raw_rfu_1,
      zeroed_rfu_1,
      ng_well_1,
      ng_ul_1,
      raw_rfu_2,
      zeroed_rfu_2,
      ng_well_2,
      ng_ul_2,
      difference,
      average,
      pass_redo)


  # Step 3: Export in workbook
  wb <- createWorkbook()

  addWorksheet(wb, "quant_data")
  writeData(wb, "quant_data", quant_output)

  addWorksheet(wb, "quant_standard")
  writeData(wb, "quant_standard", quant_standard)

  print (quant_plot)
  insertPlot(wb,
             sheet = "quant_standard",
             startRow = 1,
             startCol = 10,
             width = 6,
             height = 5,
             units = "in",
             dpi = 300,
             fileType = "png")

  addWorksheet(wb, "quant_report")
  writeData(wb, "quant_report", quant_report)

  # Save the workbook
  saveWorkbook (wb, "quant_report.xlsx", overwrite = T)

  message("Quant report exported successfully!")

  return(quant_report)



}
