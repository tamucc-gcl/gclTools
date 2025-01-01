export_quant_report <- function (quant_dna) {

  # Step 1: Load required libraries
  required_packages <- c("openxlsx")

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

  # Step 2: Extract the quant output
  quant_output <- quant_dna$quant_output

  # Step 3: Wrangle the output_quant
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


  # Step 4: Export in workbook
  wb <- createWorkbook()

  addWorksheet(wb, "quant_data")
  writeData(wb, "quant_data", quant_output)

  addWorksheet(wb, "quant_report")
  writeData(wb, "quant_report", quant_report)

  # Save the workbook
  saveWorkbook (wb, "~/quant_report.xlsx", overwrite = T)

  message("Quant report exported successfully!")

  return(quant_report)



}
