# Read the raw quant text file exported from SpectraMax


require(tidyverse)

filepath <- "inst/extdata/Tepolt_mesopelagic_fish_plate1_accublue_2025-01-10.txt"

filepath <- "inst/extdata/raw_data.txt"

# Function to find the section starting with the header
extract_section <- function(file_path, header) {
  # Open the file and read it line by line
  con <- file(file_path, open = "r")
  lines <- readLines(con)
  close(con)

  # Find the line number where the header is located
  header_line <- grep(header, lines)

  # If header is found, extract lines after it
  if (length(header_line) > 0) {
    # Extract the header and lines below it
    start_line <- header_line[1]  # Take the first match of the header
    section_lines <- lines[start_line:length(lines)]

    # Remove any non-data lines that might follow the header (e.g., blank lines or metadata)
    section_lines <- section_lines[section_lines != ""]

    return(section_lines)
  } else {
    return(NULL)  # If header is not found
  }
}

# Example usage
file_path <- filepath
header <- "Sample\tWells\tValue\tR\tResult\tMeanResult\tSD\tCV"

section_data <- extract_section(file_path, header)

# Print the extracted section data
if (!is.null(section_data)) {
  print(section_data)
} else {
  print("Header not found.")
}
