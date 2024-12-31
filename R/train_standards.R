train_standards <- function (quant_data, model = c("kit_reference", "linear", "power")) {

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

  ## Step 4: Fit the model

  model <- match.arg (model)  ## Ensure the model option is valid; this also takes the first argument as the default

  fit <- switch(
    model,
    kit_reference = lm (dna_per_well ~ 0 + zeroed_rfu, data = std),
    linear = lm (dna_per_well ~ zeroed_rfu, data = std),
    power = lm (log(dna_per_well) ~ log(zeroed_rfu), data = filter (std, dna_per_well > 0))
  )

  fitted_values <- if (model == "power") {
    c(NA, exp(fit$fitted.values))
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

  ## Construct regression equation string for plot subtitle
  subtitle <- if (model == "power") {
    # Convert to y = a * x^b form for power model
    a <- exp(coefficients[1])  # a = exp(log(a))
    b <- coefficients[2]       # b remains the same
    paste0(
      "Model: y = ", round(a, 4), " * x^", round(b, 4), "; ",
      "R² = ", round(r_squared, 4)
    )
  } else if (model == "linear") {
    paste0(
      "Model: y = ", round(coefficients[1], 4),
      " + ", round(coefficients[2], 4), " * x; ",
      "R² = ", round(r_squared, 4)
    )
  } else {
    paste0(
      "Model: y = ", round(coefficients[1], 4),
      " * x; ",
      "R² = ", round(r_squared, 4)
    )
  }

  ## Generate plot for visualization
  gg <- ggplot (std,
                aes (x = zeroed_rfu, y = dna_per_well)) +
    geom_point(color = "blue", size = 3, alpha = 0.7) +
    geom_line (aes(y = fit), color = "red", size = 1) +
    geom_text(aes(label = sample_id), vjust = -1, hjust = 0.5, size = 4, color = "black") +
    labs (title = paste("Fitted", model, "model"),
          subtitle = subtitle,
          x = "Zeroed RFU",
          y = "DNA concentration (ng) per well") +
    theme_bw()
  print (gg)


  return (list (
    standard = std,
    model_fit = fit,
    r_squared = r_squared)
  )

}



