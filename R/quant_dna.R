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


