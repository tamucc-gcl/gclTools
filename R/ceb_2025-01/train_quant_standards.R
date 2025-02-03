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
train_quant_standards <-
  function (
    quant_data,
    model =
      c(
        "kit_reference",
        "linear",
        "linear_nozerostd_zeroedrfu",
        "power",
        "power_rawrfu"
      ),
    remove_standard = FALSE
  ) {

    # Step 1: Extract the standard data frame and modify columns
    std <-
      quant_data$standard %>%
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

      std <-
        std %>%
        group_by (plate_id, sample_id) %>%
        summarize (rfu = mean(rfu, na.rm = T),
                   dna_per_ul = mean(dna_per_ul, na.rm = T),
                   dna_per_well = mean(dna_per_well, na.rm = T)
        ) %>%
        ungroup()
    }

    # Step 3: Extract and subtract the background RFU from the rest of the data
    background_rfu <-
      std %>%
      filter (dna_per_ul == 0) %>%
      pull (rfu)

    std <-
      std %>%
      mutate (
        zeroed_rfu = rfu - background_rfu
      )

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

    # filter std data for desired model
    model <- match.arg (model)  ## Ensure the model option is valid; this also takes the first argument as the default

    if (model == "linear_nozerostd_zeroedrfu" | str_detect(model,"power")){
      std <-
        std %>%
        filter(sample_id != 0)
    }

    print(std)
    # Step 5: Fit the model


    fit <-
      switch(
        model,

        # zero-intercept kit_reference model
        kit_reference = lm (dna_per_well ~ 0 + zeroed_rfu, data = std),

        # linear model (kit_reference without forced 0 RFU at zero dna conc)
        linear =
          lm (
            dna_per_well ~ zeroed_rfu,
            data = std
          ),

        # linear model nozerostd (linear model without the zero dna conc std )
        linear_nozerostd_zeroedrfu =
          lm (
            dna_per_well ~ zeroed_rfu,
            data = std
          ),

        # power model; requires that only zeroed_rfu > 0 are retained.
        power = {

          # # Count rows with zeroed_rfu <=0
          # invalid_count <- sum (std$zeroed_rfu <= 0, na.rm = T)
          # if (invalid_count > 0) {
          #   warning(invalid_count, " standards with zeroed_rfu <=0 were removed for the power model. See plot for more details.")
          # }
          #
          # # Filter rows with zeroed_rfu > 0
          # valid_data <- filter (std, zeroed_rfu > 0)
          # if (nrow(valid_data) == 0) {
          #   stop("No valid rows with zeroed_rfu > 0 for the power model.")
          # }
          # lm (log(dna_per_well) ~ log(zeroed_rfu), data = valid_data)


          lm (log(dna_per_well) ~ log(zeroed_rfu), data = std)

        },
        power_rawrfu = {

          # # Count rows with rfu <=0
          # invalid_count <- sum (std$rfu <= 0, na.rm = T)
          # # invalid_count <- sum(std$rfu <= 0 | std$dna_per_well <= 0, na.rm = TRUE)
          # print(invalid_count)
          #
          # if (invalid_count > 0) {
          #   warning(invalid_count, " standards with rfu <=0 were removed for the power_rawrfu model. See plot for more details.")
          # }
          #
          # # Filter rows with rfu > 0
          # valid_data <- filter (std, rfu > 0)
          # # valid_data <- filter(std, rfu > 0, dna_per_well > 0)
          # print(valid_data)
          #
          # if (nrow(valid_data) == 0) {
          #   stop("No valid rows with rfu > 0 for the power_rawrfu model.")
          # }
          # lm (log(dna_per_well) ~ log(rfu), data = valid_data)

          lm (log(dna_per_well) ~ log(rfu), data = std)
        }
      )
    print(fit)

    # fitted_values <-
    if (str_detect(model, "power")) {
      # For power model, expand fitted values back to the full dataset
      fitted_values <<- exp(fit$fitted.values)
      # fitted_full <- rep(NA, nrow(std))
      # fitted_full[std$zeroed_rfu > 0] <- fitted
      # fitted_full
      #
      #       } else if (model == "power_rawrfu") {
      #         # For power model, expand fitted values back to the full dataset
      #         fitted <- exp(fit$fitted.values)
      #         fitted_full <- rep(NA, nrow(std))
      #         fitted_full[std$rfu > 0] <- fitted
      #         fitted_full
    } else {
      # fit$fitted.values
      fitted_values <<- fit$fitted.values
    }
    print(fitted_values)

    ## Add fitted values to data
    std_fits <-
      std %>%
      # mutate (fit = ifelse (row_number () <= length(fitted_values), fitted_values, NA))
      mutate (fit = fitted_values)

    print(std_fits)
    ## Extract coefficients and R2
    coefficients <- coef (fit)
    r_squared <- summary(fit)$r.squared


    ## Provide feedback
    message ("Fitted a ", model, " model.")
    message("Coefficients: ", paste(names(coefficients), round(coefficients, 4), sep = " = ", collapse = ", "))
    message("R-squared: ", round(r_squared, 4))

    # Step 6. Plot the model

    ## Construct regression equation string for plot subtitle
    subtitle <-
      if (model == "power") {
        # Convert to y = a * x^b form for power model
        a <- exp(coefficients[1])  # a = exp(log(a))
        b <- coefficients[2]       # b remains the same
        paste0(
          "Model: y = ", round(a, 4), " * x^", round(b, 4), "; ",
          "R2 = ", round(r_squared, 4)
        )
      } else if (model == "power_rawrfu") {
        # Convert to y = a * x^b form for power model
        a <- exp(coefficients[1])  # a = exp(log(a))
        b <- coefficients[2]       # b remains the same
        paste0(
          "Model: y = ", round(a, 4), " * x^", round(b, 4), "; ",
          "R2 = ", round(r_squared, 4)
        )
      }
    else if (model == "linear") {
      paste0(
        "Model: y = ", round(coefficients[1], 4),
        " + ", round(coefficients[2], 4), " * x; ",
        "R2 = ", round(r_squared, 4)
      )
    } else if (model == "linear_nozerostd_zeroedrfu") {
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


    #### plotting function CEB ####

    make_regression_plot <-
      function(
    data,
    x_var,       # character string for the column name to use on x-axis
    y_var,       # character string for the column name to use on y-axis
    fit_var = "fit",     # character string for the column containing fitted values
    model_fit = fit, # pass an 'lm' or 'nls' model if you want a smooth line
    point_color = "blue",
    line_color  = "red",
    point_size  = 3,
    line_size   = 1,
    alpha       = 0.7,
    plot_title  = "Fitted Model",
    plot_subtitle = NULL,
    x_label     = NULL,
    y_label     = NULL
      ) {

        # Convert the column names from character strings into something ggplot can use.
        # aes_string() or .data pronoun can be used to refer to columns by string.
        p <-
          ggplot(
            data,
            aes_string(
              x = x_var,
              y = y_var)
          ) +
          # Scatter points
          geom_point(color = point_color, size = point_size, alpha = alpha) +
          # Fitted line
          # geom_line(
          #   aes_string(y = fit_var),
          #   color = line_color,
          #   size  = line_size
          # ) +
          labs(
            title    = plot_title,
            subtitle = plot_subtitle,
            x        = x_label %||% x_var,
            y        = y_label %||% y_var
          ) +
          theme_bw()


        # 1. Extract x values (optionally removing NAs or negative values if needed)
        xvals <- data[[x_var]]
        xvals <- xvals[!is.na(xvals)]         # remove NA
        # xvals <- xvals[xvals > 0]           # optionally remove x <= 0, if desired

        # 2. Sort the x values
        xvals_sorted <- sort(xvals)

        # 3. Check that we actually have at least 2 unique values
        if (length(xvals_sorted) < 2) {
          stop("Cannot find a second‐smallest x-value. Need at least 2 distinct x-values.")
        }

        # 4. Define x_min (the absolute smallest) and second_smallest (the upper bound)
        x_min          <- xvals_sorted[1]
        x_second_smallest <- xvals_sorted[2]
        x_max <- max(xvals_sorted)

        # x_min <- min(data[[x_var]], na.rm = TRUE)
        # x_max <- max(data[[x_var]], na.rm = TRUE)
        print(x_min)
        # print(x_max)
        print(x_second_smallest)

        # 3b. Generate a new data frame with n_points equally spaced x-values
        # generate several values, but remove zeros and negs
        new_x    <-
          c(
            seq(x_min, x_max, length.out = 100)[seq(x_min, x_max, length.out = 100) > 0],
            seq(x_min, x_second_smallest, length.out = 1000)[seq(x_min, x_second_smallest, length.out = 1000) > 0]
          ) %>%
          sort() %>%
          unique()

        print(new_x)
        new_data <- data.frame(x_val = new_x)
        # rename 'x_val' to the actual column name so predict() can see it
        names(new_data) <- x_var
        print(new_data)

        # 3c. Predict on this new_data
        preds <- predict(model_fit, newdata = new_data)

        # If it's a log-based model (like lm(log(y) ~ log(x))),
        # then 'preds' is on the log scale.
        # So we exponentiate:
        # preds <- exp(preds)

        # 3d. Attach preds as a column in 'new_data'
        new_data$predicted <- preds

        # 3e. Add another line to the plot
        p <-
          p + geom_line(
            data = new_data,
            aes_string(x = x_var, y = "predicted"),
            color = "darkgreen",
            size  = line_size,
            linetype = "dashed"
          )
      }

    #### PLOT STANDARD CURVE CEB ####

    plot_title <-
      paste(
        "Fitted",
        model,
        "model"
      )

    plot_subtitle <-
      paste0(
        subtitle,          # your existing model equation string
        "\nblue points = standards; grey polygons = samples"
      )

    # Decide on x_var based on the model chosen
    if (model == "power_rawrfu") {
      x_for_plot <- "rfu"
      x_label    <- "Raw RFU"
      y_for_plot <- "dna_per_well"
      y_label    <- "DNA concentration (ng/well)"
    }
    else if (model == "power") {
      x_for_plot <- "zeroed_rfu"
      x_label    <- "Zeroed RFU"
      y_for_plot <- "dna_per_well"
      y_label    <- "DNA concentration (ng/well)"
    }
    else {
      x_for_plot <- "zeroed_rfu"
      x_label    <- "Zeroed RFU"
      y_for_plot <- "dna_per_well"
      y_label    <- "DNA concentration (ng/well)"
    }

    gg <-
      make_regression_plot(
        data         = std_fits,
        x_var        = x_for_plot,
        y_var        = y_for_plot,
        fit_var      = "fit",
        model_fit    = fit,
        point_color  = "blue",
        line_color   = "red",
        plot_title   = plot_title,
        plot_subtitle= plot_subtitle,
        x_label      = x_label,
        y_label      = y_label
      )
    print (gg)

    ## Generate plot for visualization
    # gg <- ggplot (std,
    #               aes (x = zeroed_rfu, y = dna_per_well)) +
    #   geom_point(color = "blue", size = 3, alpha = 0.7) +
    #   geom_line (aes(y = fit), color = "red", size = 1) +
    #   geom_text(aes(label = sample_id), vjust = -1, hjust = 0.5, size = 4, color = "black") +
    #   labs (title = paste("Fitted", model, "model"),
    #         subtitle = paste0(subtitle, "\nblue points = standards; grey polygons = samples"),
    #         x = "Zeroed RFU",
    #         y = "DNA concentration (ng/well)") +
    #   theme_bw()
    # print (gg)



    ## Extract function output

    return (list (
      background_rfu = background_rfu,
      standard = std_fits,
      model_fit = fit,
      r_squared = r_squared,
      plot = gg)
    )

  }



