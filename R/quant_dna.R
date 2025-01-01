#' Analyze quant fluorescence data
#'
#' @description
#' Use a trained regression model to predict DNA concentrations for sample data in a fluorescent quant assay. This function processes the sample data from the quant assay, calculates zeroed RFUs, predicts DNA concentration (ng/well) based on the regression model, and normalizes it based on the sample volume (ng/uL). This takes the object generated via `train_quant_standards` as input.
#'
#' This is Step 3 of 4 in a series of functions for analyzing fluorescent data.
#'
#' @import tidyverse
#'
#' @param quant_data
#' list. Output of `load_quant_files`, containing both the raw data and plate map information.
#'
#' @param trained_model
#' list. Output of `train_quant_standards`, containing the trained regression model, background RFU, and a base ggplot object.
#'
#' @details
#' The function uses the background RFU from the trained model to calculate zeroed RFUs for the sample data. Predictions are based on the regression model from `trained_model`. For power models, predictions are transformed back from log space.
#'
#' @returns
#' Returns a list containing:
#' - `quant_output`: data frame. Contains the sample data with predicted DNA concentrations (ng/well and ng/uL).
#' - `plot`: ggplot object. Updated visualization of the regression model with sample predictions overlaid as jittered points.

#'
#' @examples
#' # Import data files
#' raw_data <- system.file("extdata", "raw_data.csv", package = "tamuccGCL")
#' plate_map <- system.file("extdata", "plate_map.csv", package = "tamuccGCL")
#'
#' quant_data <- load_quant_files(raw_data, plate_map)
#' trained_model <- train_quant_standards (quant_data)
#' quant_report <- quant_dna(quant_data = quant_data, trained_model = trained_model)
#
quant_dna <- function (quant_data, trained_model) {

  # Step 1: Extract the necessary components from trained_data
  background_rfu <- trained_model$background_rfu
  model_fit <- trained_model$model_fit
  plot <- trained_model$plot


  # Step 2: Extract the sample data and calculate zeroed rfu
  sample_data <- quant_data$sample %>%
    mutate (zeroed_rfu = rfu - background_rfu)

  # Step 3: Predict the DNA concentration per well
  ## If the model is a power model, we must transform the prediction back from log space
  if (inherits(model_fit, "lm") && any(grepl("log", formula(model_fit)))) {
    sample_data <- sample_data %>%
      mutate (dna_per_well = exp(predict(model_fit, newdata = .)))
  } else {
    sample_data <- sample_data %>%
      mutate (dna_per_well = predict (model_fit, newdata = .))
  }

  # Step 4: Calculate dna_per_ul
    sample_data <- sample_data %>%
      mutate (dna_per_ul = dna_per_well / sample_volume)

  # Step 5: Add predicted values to plot
  gg <- plot +
    geom_jitter (
      data = sample_data,
      aes(x = zeroed_rfu, y = dna_per_well),
      fill = "gray50",
      col = "black",
      pch = 23,
      size = 3,
      alpha = 0.75 )

  # Extract function output
  return (list(
    quant_output = sample_data,
    plot = gg)
  )
}


