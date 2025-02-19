##### Remove replicates with high CVs ####

remove_high_cv_replicates <- function(
    data,
    group_var = "sample_id",
    response_var = "rfu",
    cv_threshold = 15
) {
  # Ensure required packages are loaded
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required")
  }

  data %>%
    group_by(across(all_of(group_var))) %>%
    group_modify(~ {
      group_data <- .x
      n_reps <- nrow(group_data)

      # Only process groups with >2 replicates
      if (n_reps <= 2) return(group_data)

      # Calculate original CV
      original_cv <- (sd(group_data[[response_var]], na.rm = TRUE) /
                        mean(group_data[[response_var]], na.rm = TRUE)) * 100

      # Skip if original CV is already acceptable
      if (original_cv <= cv_threshold) return(group_data)

      # Calculate CV when removing each replicate
      cv_values <- sapply(1:n_reps, function(i) {
        remaining <- group_data[[response_var]][-i]
        (sd(remaining)/mean(remaining)) * 100
      })

      # Find which removal gives best CV improvement
      best_removal <- which.min(cv_values)

      # Only remove if improvement meets threshold
      if (cv_values[best_removal] <= cv_threshold) {
        return(group_data[-best_removal, ])
      } else {
        return(group_data)  # Keep original if no single removal helps enough
      }
    }) %>%
    ungroup()
}


#### Report Standards ####
report_standards <-
  function(
    quant_plate_map,
    quant_kit
  ) {
    # Count the number of standards
    standard_count <-
      quant_plate_map %>%
      filter(grepl("standard", plate_id, ignore.case = TRUE)) %>%
      nrow()
    message(paste(standard_count, "standards detected"))

    # Extract the standard sample IDs
    standard_samples <-
      quant_plate_map %>%
      filter(grepl("standard", plate_id, ignore.case = TRUE)) %>%
      pull(sample_id)

    # Determine the standard unit based on the quant_kit
    standard_unit <- if (quant_kit == "accublue-nextgen") "pg/ul" else "ng/ul"
    message(paste("The standards used (in", standard_unit, ") are:", paste(standard_samples, collapse = ", ")))

    # Return a list with the results
    list(
      standard_count = standard_count,
      standard_samples = standard_samples,
      standard_unit = standard_unit
    )
  }

##### Regression model without samples ####
# Fit regression model and plot
fit_and_plot <-
  function(
    x_var,
    y_var,
    data,
    model_type,
    degree = NULL
  ) {

    model_type <-
      match.arg(
        model_type,
        c("Linear", "Power", "Polynomial"))

    # Define model formula
    model_formula <-
      switch(
        model_type,
        "Linear" = as.formula(paste(y_var, "~", x_var)),
        "Power" = as.formula(paste("log10(", y_var, ") ~ log10(", x_var, ")")),
        "Polynomial" = as.formula(paste("log10(", y_var, ") ~ poly(log10(", x_var, "), degree =", degree, ")"))
      )

    # Fit model
    fit <- lm(model_formula, data = data %>% filter(include_in_model))

    # Summary statistics
    model_summary <- summary(fit)
    aic <- round(AIC(fit), 4)
    bic <- round(BIC(fit), 4)
    rse <- round(model_summary$sigma, 4)
    subtitle <- as.character(formula(fit))

    # Filter data to include only standards marked for modeling
    filtered_data <- data %>% filter(include_in_model)

    # Generate prediction data
    pred_data <-
      data.frame(
        x_var = seq(min(filtered_data[[x_var]], na.rm = TRUE),
                    max(filtered_data[[x_var]], na.rm = TRUE),
                    length.out = 100)
      ) %>%
      setNames(x_var)

    # Apply transformations as needed
    pred_data <-
      pred_data %>%
      mutate(
        log_rfu = log10(.data[[x_var]]),
        log_dna_pred = predict(fit, newdata = .),
        dna_pred = if (model_type == "Linear") log_dna_pred else 10^log_dna_pred
      )

    # Generate plot
    plot_model <-
      ggplot(
        data,
        aes(x = log10(.data[[x_var]]),
            y = log10(.data[[y_var]]))
      ) +

      geom_errorbarh(
        aes(xmin = log10(.data[[x_var]] - rfu_sd),
            xmax = log10(.data[[x_var]] + rfu_sd),
            color = include_in_model),
        height = 0.15) +

      geom_point(
        aes(shape = "Standard",
            color = include_in_model),
        size = 2) +

      geom_text(
        aes(label = .data[[y_var]]),
        size = 3.5,
        nudge_x = 0.01,
        nudge_y = 0.1) +

      geom_line(
        data = pred_data,
        aes(x = log_rfu,
            y = log10(dna_pred)),
        color = "red",
        lwd = 1) +

      labs(
        title = paste(model_type, "Model"),
        subtitle =
          paste("Model:",
                subtitle[2],
                subtitle[1],
                subtitle[3],
                "\nRSE =", rse,
                "; AIC =", aic,
                "; BIC =", bic),
        x = expression(log[10] * "(RFU)"),
        y = if (quant_kit == "accublue-nextgen") {
          expression(log[10] * "(DNA concentration; pg/well)")
        } else {
          expression(log[10] * "(DNA concentration; ng/well)")
        }
      ) +
      theme_bw()

    list(fit = fit, plot = plot_model)
  }


##### Add samples to regression plot ####
plot_regression_with_samples <-
  function(
    base_plot,
    model_type,
    min_rfu,
    max_rfu,
    predicted_sample_concentration
  ) {
    standard_unit <- if (quant_kit == "accublue-nextgen") "pg/well" else "ng/well"

    base_plot +

      geom_point(
        data = predicted_sample_concentration %>% filter(model == model_type),
        aes(x = log10(rfu),
            y = log10(dna_per_well),
            shape = "Sample"),
        alpha = 0.50,
        size = 2
      ) +

      labs(color = "Included in model") +

      scale_shape_manual(
        name = "Classification",
        values = c("Standard" = 15, "Sample" = 16),
        labels = c("Standard" = paste0("Standard (", standard_unit, ")"), "Sample" = "Sample")
      ) +

      geom_vline(xintercept = log10(min_rfu),
                 linetype = 2) +

      # annotate("text",
      #          x = log10(min_rfu) - 0.025,
      #          y = log10(max_rfu) + 0.5,
      #          label = "below lod",
      #          angle = 90,
      #          size = 3) +

      geom_vline(xintercept = log10(max_rfu),
                 linetype = 2)

      # annotate("text",
      #          x = log10(max_rfu) + 0.025,
      #          y = 1,
      #          label = "above lod",
      #          angle = 90,
      #          size = 3)
  }


##### Compare, rank, and select best model ####
rank_models <- function(...) {
  models <- list(...)
  model_names <- names(models)

  if (is.null(model_names)) {
    stop("Please provide named arguments for each model (e.g., linear = fit_linear).")
  }

  metrics_df <- tibble(
    model = model_names,
    AIC = map_dbl(models, AIC),
    BIC = map_dbl(models, BIC),
    RSE = map_dbl(models, ~ summary(.x)$sigma)
  ) %>%
    mutate(
      AIC_rank = dense_rank(AIC),
      BIC_rank = dense_rank(BIC),
      RSE_rank = dense_rank(RSE),
      total_rank = AIC_rank + BIC_rank + RSE_rank
    )

  model_priority <- c("power", "polynomial", "linear")
  metrics_df %>%
    mutate(
      model_factor = factor(model,
                            levels = model_priority,
                            ordered = T)
      ) %>%
    arrange(
      total_rank,
      BIC,
      AIC,
      RSE,
      model_factor
      )
}



# Predict sample concentration
##### Predict sample concentration ####
predict_sample_concentration <-
  function(
    data,
    model,
    model_type
  ) {
    data %>%
      rename(rfu_mean = rfu) %>%
      mutate(
        dna_per_well = if (model_type == "linear") predict(model, newdata = .) else 10^predict(model, newdata = .),
        pg_per_ul = if (quant_kit == "accublue-nextgen") dna_per_well / sample_volume else NA,
        ng_per_ul = if (quant_kit == "accublue-nextgen") pg_per_ul / 1000 else dna_per_well / sample_volume,
        model = model_type
      )
  }


#### Identify Outlier Standards #### - JDS
identify_outliers <- function(x, test = dixon.test, alpha = 0.05, ...){
  # Returns a logical vector indicating which values are outliers found by the specified test (either dixon.test or grubbs.test) at the specified alpha.
  # Specifically testing for outliers above the average

  p <- 1/Inf
  outlier_standards <- logical(length(x))

  while(p < alpha){
    test_out <- test(x[!outlier_standards], two.sided = FALSE, )

    p <- test_out$p.value
    if(p < alpha){
      outlier_standards[!outlier_standards][which.max(x[!outlier_standards])] <- TRUE
    }

  }

  outlier_standards
}
