#' Analyze quant fluorescence data
#'
#' @description
#' Train a regression model to predict DNA concentration (ng/well) using relative fluorescence unit (RFU) measured from the standards in the quant assay. This takes the object generated via [load_quant_files()] as input.
#'
#' This is Step 2 of 4 in a series of functions for analyzing fluorescent data.
#'
#' @import tidyverse
#'
#' @param quant_data
#' list. Output of [load_quant_files()].
#'
#' @param model
#' character string. Select the regression model to use, where \eqn{y} = DNA concentration (ng/well) and \eqn{x} = Zeroed_RFU = RFU - background RFU; background_RFU, obtained by subtracting the RFU of the zero standard against all the other standards, is done internally.
#'
#' Valid options are:
#' - "kit_reference" (\eqn{y = 0 + bx}; a linear model with zero intercept as per the instructions in the user manual. This is the default option)
#' - "linear" (\eqn{y = a + bx}; a linear model with non-zero intercept)
#' - "power" (\eqn{y = exp(a)*x^b}; a power model with the formula definition \eqn{log(y) = log(x)}; used when there is a need to normalize the standard concentrations)
#'
#' @param remove_standard
#' numeric vector. Select the standards that need to be removed to improve model fit. This is done when standards become outliers thus lowering the model's R2 value. Default is set to `FALSE` indicating that none of the standards will be removed.
#'
#' @returns
#' Returns a list containing the following:
#' - `background_rfu`: numeric. The RFU of the zero-standard.
#' - `standard`: data frame. Contains the standard data used in the regression.
#' - `model_fit`: the summary of the regression model, as though running `summary(model)`.
#' - `r_squared`: numeric. The R squared value of the regression model.
#' - `plot`: ggplot visualization of the regression showing the standards (point) and the model fit (regression line).
#'
#' @details
#' Note that when running the power model, the standards with zeroed RFUs \eqn{<= 0} are removed to make the model work.
#'
#' @seealso [load_quant_files()]
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
#' # Train using the zero-intercept model (default)
#' trained_model <- train_quant_standards (quant_data)
#'
#' # Train using the power model
#' trained_model <- train_quant_standards (quant_data, model = "power")
#'
#' # Train using the linear model with the 10 ng/uL standard removed
#' trained_model <- train_quant_standards (quant_data, model = "linear", remove_standard = 10)
#'
#' @export
train_quant_standards <- function (quant_data, model = c("kit_reference", "linear", "power"), remove_standard = FALSE) {

  # Step 1: Extract the standard data frame and modify columns
  std <- quant_data$standard %>%
    mutate (sample_id = as.numeric(sample_id),
            dna_per_ul = as.numeric (sample_id),
            dna_per_well = dna_per_ul * sample_volume) %>%

    select (plate_id, sample_id, replicate, rfu, dna_per_ul, dna_per_well)

  # Step 2: Calculate the mean of replicates if necessary.
  replicate_max <- max(std$replicate, na.rm = T)

  if (replicate_max == 1) {
    message ("No replicate standards detected. Original RFU values are kept.")
    std <- std %>% select (-replicate)

  } else {
    message ("Replicate standards detected. All values shown hereafter are the mean of the replicates.")

    std <- std %>%
      group_by (plate_id, sample_id) %>%
      summarize (rfu = mean(rfu, na.rm = T),
                 dna_per_ul = mean(dna_per_ul, na.rm = T),
                 dna_per_well = mean(dna_per_well, na.rm = T)
      ) %>%
      ungroup()
  }

  # Step 3: Extract and subtract the background RFU from the rest of the data
  background_rfu <- std %>%
    filter (dna_per_ul == 0) %>%
    pull (rfu)

  std <- std %>%
    mutate (zeroed_rfu = rfu - background_rfu)

  # Step 4 [OPTIONAL]: Remove standard to improve model fit
  if (!isFALSE (remove_standard)) {
    # Check if the input is a single value or a vector
    standards_to_remove <- as.character (remove_standard)

    # Validate that all specified standards exist in the dataset
    missing_standards <- setdiff(standards_to_remove, std$sample_id)
    if (length(missing_standards) > 0) {
      stop ("The following standard/s to remove do not exist in the dataset: ",
            paste(missing_standards, collapse = ", "))
    }

    # Filter out the specified standards
    message("Removing the following standards: ", paste(standards_to_remove, collapse = ", "))
    std <- std %>% filter(!sample_id %in% standards_to_remove)
  }

  # Step 5: Fit the model

  model <- match.arg (model)  ## Ensure the model option is valid; this also takes the first argument as the default

  fit <- switch(
    model,

    # zero-intercept linear model
    kit_reference = lm (dna_per_well ~ 0 + zeroed_rfu, data = std),

    # non-zero intercept linear model
    linear = lm (dna_per_well ~ zeroed_rfu, data = std),

    # power model; requires that only zeroed_rfu > 0 are retained.
    power = {

      # Count rows with zeroed_rfu <=0
      invalid_count <- sum (std$zeroed_rfu <= 0, na.rm = T)
      if (invalid_count > 0) {
        warning(invalid_count, " standards with zeroed_rfu <=0 were removed for the power model. See plot for more details.")
      }

      # Filter rows with zeroed_rfu > 0
      valid_data <- filter (std, zeroed_rfu > 0)
      if (nrow(valid_data) == 0) {
        stop("No valid rows with zeroed_rfu > 0 for the power model.")
      }
      lm (log(dna_per_well) ~ log(zeroed_rfu), data = valid_data)
    }
  )

  fitted_values <- if (model == "power") {
    # For power model, expand fitted values back to the full dataset
    fitted <- exp(fit$fitted.values)
    fitted_full <- rep(NA, nrow(std))
    fitted_full[std$zeroed_rfu > 0] <- fitted
    fitted_full
  } else {
    fit$fitted.values
  }

  ## Add fitted values to training data
  std <- std %>%
    mutate (fit = ifelse (row_number () <= length(fitted_values), fitted_values, NA))

  ## Extract coefficients and R2
  coefficients <- coef (fit)
  r_squared <- summary(fit)$r.squared


  ## Provide feedback
  message ("Fitted a ", model, " model.")
  message("Coefficients: ", paste(names(coefficients), round(coefficients, 4), sep = " = ", collapse = ", "))
  message("R-squared: ", round(r_squared, 4))

  # Step 6. Plot the model

  ## Construct regression equation string for plot subtitle
  subtitle <- if (model == "power") {
    # Convert to y = a * x^b form for power model
    a <- exp(coefficients[1])  # a = exp(log(a))
    b <- coefficients[2]       # b remains the same
    paste0(
      "Model: y = ", round(a, 4), " * x^", round(b, 4), "; ",
      "R2 = ", round(r_squared, 4)
    )
  } else if (model == "linear") {
    paste0(
      "Model: y = ", round(coefficients[1], 4),
      " + ", round(coefficients[2], 4), " * x; ",
      "R2 = ", round(r_squared, 4)
    )
  } else {
    paste0(
      "Model: y = ", round(coefficients[1], 4),
      " * x; ",
      "R2 = ", round(r_squared, 4)
    )
  }

  ## Generate plot for visualization
  gg <- ggplot (std,
                aes (x = zeroed_rfu, y = dna_per_well)) +
    geom_point(color = "blue", size = 3, alpha = 0.7) +
    geom_line (aes(y = fit), color = "red", size = 1) +
    geom_text(aes(label = sample_id), vjust = -1, hjust = 0.5, size = 4, color = "black") +
    labs (title = paste("Fitted", model, "model"),
          subtitle = paste0(subtitle, "\nblue points = standards; grey polygons = samples"),
          x = "Zeroed RFU",
          y = "DNA concentration (ng/well)") +
    theme_bw()
  print (gg)

  ## Extract function output

  return (list (
    background_rfu = background_rfu,
    standard = std,
    model_fit = fit,
    r_squared = r_squared,
    plot = gg)
  )

}



