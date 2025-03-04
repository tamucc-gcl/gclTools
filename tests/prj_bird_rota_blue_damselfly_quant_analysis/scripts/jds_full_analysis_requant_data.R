library(tidyverse)
library(readxl)
library(janitor)
library(outliers)
library(glmmTMB)

source ("helper_functions.R")

#Tweaked/new functions to make work
read_quant_data <- function(path_raw_data, path_plate_map){
  quant_raw_data <-
    read_csv(path_raw_data, show_col_types = FALSE) %>%
    clean_names() %>%
    select (wells, value) %>%
    mutate(across(where(is.character), str_to_lower))

  # Quant plate map
  quant_plate_map <-
    read_csv(path_plate_map, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate (
      wells = paste0(quant_row, quant_column),
      sample_id = trimws(sample_id),
      across(where(is.character), str_to_lower),
      # If sample_id starts with digits followed by letters (e.g. "12b"), swap them to become "b12"
      sample_id = str_replace(sample_id, "^([0-9]+)([a-z]+)$", "\\2\\1")
    )


  ## Merge data frames
  quant_data <-
    left_join(
      quant_plate_map,
      quant_raw_data,
      by = "wells"
    ) %>%
    rename (rfu = value) %>%
    mutate (category = case_when (
      plate_id == "standard" ~ "standard",
      T ~ "sample")) %>%
    split (.$category) %>%
    map (., ~select (.x, -category))

  quant_data
}

fit_and_plot <- function(
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

    # geom_errorbarh(
    #   aes(xmin = log10(.data[[x_var]] - rfu_sd),
    #       xmax = log10(.data[[x_var]] + rfu_sd),
    #       color = include_in_model),
    #   height = 0.15) +

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


fit_models <- function(x_var, y_var, data){
  ## Fit big models
  fit_linear <-
    fit_and_plot(
      x_var = x_var,
      y_var = y_var,
      data = data,
      model_type = "Linear",
      degree = NULL
    )

  fit_power <-
    fit_and_plot(
      x_var = x_var,
      y_var = y_var,
      data = data,
      model_type = "Power",
      degree = NULL
    )

  fit_polynomial <-
    fit_and_plot(
      x_var = x_var,
      y_var = y_var,
      data = data,
      model_type = "Polynomial",
      degree = 2
    )

  fit_polynomial3 <-
    fit_and_plot(
      x_var = x_var,
      y_var = y_var,
      data = data,
      model_type = "Polynomial",
      degree = 3
    )

  ## Jacknife to check bad fits
  jacknife_standards <- mutate(data,
                               sample_id = row_number()) %>%
    filter(include_in_model) %>%

    expand_grid(tibble(model = c('Linear', 'Power', 'Polynomial', 'Polynomial'),
                       degree = c(NA_integer_, NA_integer_, 2, 3),
                       full_model = list(fit_linear$fit, fit_power$fit, fit_polynomial$fit, fit_polynomial3$fit))) %>%
    nest(removed_standard = -c(sample_id, model, degree, full_model)) %>%
    rowwise %>%

    #Make dataset with one standard removed
    mutate(other_standards = list(anti_join(data,
                                            removed_standard,
                                            by = colnames(removed_standard)))) %>%

    #Fit model without standard
    mutate(jacknife_model = list(fit_and_plot(x_var = x_var,
                                              y_var = y_var,
                                              data = other_standards,
                                              model_type = model,
                                              degree = degree) %>%
                                   pluck('fit'))) %>%


    #Calculate predicted values with/without point
    mutate(observed_standard = removed_standard[["dna_per_well"]],

           across(c(full_model, jacknife_model),
                  list(prediction = ~predict(., newdata = removed_standard)),
                  .names = "{stringr::str_replace(.col, 'model', .fn)}")) %>%
    ungroup %>%

    #Calculate relative difference between the standard's known value and the predicted DNA
    mutate(across(ends_with('prediction'),
                  ~sqrt((. - observed_standard)^2) / observed_standard,
                  .names = "{stringr::str_replace(.col, 'prediction', 'relDiff')}"),
           .keep = 'unused') %>%

    #Smaller is better - larger differences mean relative to the y-value the predicted value without that point was much different than it was with that point
    mutate(rel_improvement = sqrt((full_relDiff - jacknife_relDiff)^2),
           .keep = 'unused') %>%
    select(sample_id, model, degree, rel_improvement)

  ## Refit without outliers
  standards_without_outliers <- jacknife_standards %>%
    mutate(is_outlier = identify_outliers(rel_improvement),
           .by = c('model', 'degree'), .keep = 'unused') %>%
    filter(all(is_outlier),
           .by = 'sample_id') %>%
    select(-is_outlier) %>%
    summarise(outlier_standards = list(c(sample_id) %>% unique)) %>%
    rowwise %>%
    mutate(new_standards = list(mutate(data,
                                       sample_id = row_number(),
                                       include_in_model = !sample_id %in% outlier_standards & include_in_model)),
           .keep = 'unused') %>%
    ungroup %>%
    unnest(new_standards)

  fit_linear <-
    fit_and_plot(
      x_var = x_var,
      y_var = y_var,
      data = standards_without_outliers,
      model_type = "Linear",
      degree = NULL
    )

  fit_power <-
    fit_and_plot(
      x_var = x_var,
      y_var = y_var,
      data = standards_without_outliers,
      model_type = "Power",
      degree = NULL
    )

  fit_polynomial <-
    fit_and_plot(
      x_var = x_var,
      y_var = y_var,
      data = standards_without_outliers,
      model_type = "Polynomial",
      degree = 2
    )

  fit_polynomial3 <-
    fit_and_plot(
      x_var = x_var,
      y_var = y_var,
      data = standards_without_outliers,
      model_type = "Polynomial",
      degree = 3
    )



  tibble(model = c('linear', 'power', 'poly2', 'poly3'),
         fit = list(fit_linear$fit,
                    fit_power$fit,
                    fit_polynomial$fit,
                    fit_polynomial3$fit),
         plot = list(fit_linear$plot,
                     fit_power$plot,
                     fit_polynomial$plot,
                     fit_polynomial3$plot))
}

rank_models <- function(models) {
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

# Settings
path_dna_plate_map <-
  paste0(
    "../data_raw/",
    "rbd_edna-extraction_plate-map.xlsx"
  )
sheet_metadata <- "Plate-map-tidy"
quant_kit <- "accuclear" # Choices: accublue-nextgen or accuclear
remove_these_standards <- c(0)

path_quant_report_summarized <-
  paste0(
    "../data_processed/",
    "rbd_edna_requant_report_merged_summarized.csv"
  )


path_quant_plot <-
  paste0(
    "../results/",
    "rbd_requant_sample_assessment.png"
  )

minimum_pipettable_volume <- 0.75 #uL
maximum_low_volume <- 4
target_dna <- 2 #multiplier_on mean

#### Read in Data ####
dna_plate_map <-
  read_xlsx(
    path_dna_plate_map,
    sheet = sheet_metadata
  ) %>%
  clean_names() %>%
  mutate (wells = paste0(plate_row, plate_col)) %>%
  mutate(across(where(is.character), str_to_lower))

all_quant_data <- tibble(data_plate = list.files('../data_raw/', pattern = 'requants.*\\.csv$', full.names = TRUE),
       plate_map = list.files('../data_processed/', pattern = 'requants.*\\.csv$', full.names = TRUE)) %>%
  mutate(name = str_extract(data_plate, 'plate[0-9]+(_Column-9)?|requants')) %>%
  rowwise(name) %>%
  reframe(quant_data = read_quant_data(data_plate, plate_map)) %>%
  mutate(type = rep(c('sample', 'standard'), 1)) %>%
  pivot_wider(names_from = type, values_from = quant_data) %>%
  rowwise %>%
  mutate(standard = mutate(standard,
                           sample_id = as.numeric(sample_id),
                           dna_per_ul = as.numeric(sample_id),
                           dna_per_well = dna_per_ul * sample_volume,
                           include_in_model = !dna_per_well %in% remove_these_standards) %>%
           filter (sample_id != 0) %>%
           list) %>%
  ungroup


## Model standards & predict DNA concentrations
all_quant_models <- all_quant_data %>%
  rowwise %>%
  mutate(tmp = fit_models(data = standard,
                    x_var = "rfu",
                    y_var = "dna_per_well") %>%
           list) %>%
  unnest(tmp) %>%
  rowwise %>%
  mutate(min_rfu = min(standard$rfu[standard$include_in_model]),
         max_rfu = max(standard$rfu[standard$include_in_model]),
         sample_predictions = predict_sample_concentration(sample,
                                                           fit,
                                                           model,
                                                           x_var = 'rfu') %>%
           mutate(model = model) %>%
           list,

         plot = plot_regression_with_samples(
           base_plot = plot,
           predicted_sample_concentration = sample_predictions,
           model_type = model,
           min_rfu = min_rfu,
           max_rfu = max_rfu
         ) %>%
           list())

## ID Best model for each quant plate
investigate_standard_dropping <- all_quant_models %>%
  ungroup %>%
  group_by(name) %>%
  summarise(top_model = set_names(fit, model) %>%
                                rank_models() %>%
              slice(1) %>%
              pull(model),
            sample_predictions = bind_rows(sample_predictions) %>%
              filter(model == top_model) %>%
              list,
            plot = list(patchwork::wrap_plots(plot) +
                          patchwork::plot_annotation(title = str_c('QuantPlate: ', name, '\nTop Model: ', top_model))))


with(investigate_standard_dropping, walk2(plot, name, ~ggsave(plot = .x,
                                                              file = str_c('../results/', .y, '_standardFitting.png'),
                                                              width = 7.5,
                                                              height = 10,
                                                              units = "in",
                                                              dpi = 330)))



## Get DNA Quant Table
dna_quants <- investigate_standard_dropping %>%
  select(name, sample_predictions) %>%
  unnest(sample_predictions) %>%

  select(-wells)  %>%
  left_join (
    dna_plate_map,
    by = c(
      "sample_row" = "zymo_plate_row",
      "sample_column" = "zymo_plate_col",
      "plate_id" = "zymo_plate_id"
    )
  ) %>%
  select (
    plate_id,
    replicate,
    plate_row,
    plate_column = plate_col,
    quant_row,
    quant_column,
    sample_volume,
    model,
    rfu,
    dna_per_well,
    ng_per_ul,
    dna_extract_tube_id,
    sample_type
  ) %>%
  mutate (plate_row = toupper(plate_row),
          quant_row = toupper(quant_row)
  )

#### Model individual DNA amounts per sample #### - summarize_quant_results.R
dna_amounts <- mutate(dna_quants,
                      ID = str_c(plate_id,
                                 plate_column,
                                 plate_row,
                                 sample_type,
                                 dna_extract_tube_id,
                                 sep = ';;'),
                      is_control = str_detect(sample_type, 'control'),
                      .keep = 'unused') %>%
  select(ID, is_control, ng_per_ul) %>%
  na.omit()

dna_model_freq <- glmmTMB(ng_per_ul ~ is_control + (1 | ID),
                     dispformula = ~is_control,
                     family = Gamma(link = 'log'),
                     dna_amounts)

library(brms)
library(cmdstanr)
library(tidybayes)

get_prior(bf(ng_per_ul ~ is_control + (1 | ID),
             shape ~ is_control + (1 | ID)),
          family = Gamma(link = 'log'),
          data = dna_amounts)

dna_model_bayes <- brm(bf(ng_per_ul ~ is_control + (1 | ID),
                          shape ~ is_control + (1 | ID)),
                       prior = prior(student_t(3, -1.1, 2.5),
                                     class = 'Intercept') +
                         prior(student_t(3, 0, 2.5),
                                     class = 'Intercept',
                                     dpar = 'shape') +

                         prior(student_t(3, 0, 2.5),
                               class = 'b') +
                         prior(student_t(3, 0, 2.5),
                               class = 'sd'),

                       family = Gamma(link = 'log'),
                       data = dna_amounts,
                       backend = 'cmdstanr',
                       chains = 4,
                       iter = 5000,
                       warmup = 2000,
                       thin = 10,
                       cores = parallelly::availableCores() - 1)

plot(dna_model_bayes)
summary(dna_model_bayes)


#### Predict DNA amounts ####
quant_files_summarized <- distinct(dna_amounts,
         is_control,
         ID) %>%
  add_epred_draws(dna_model_bayes) %>%
  point_interval(.width = 0.95) %>%
  select(ID, is_control,
         ng_per_ul_mean = .epred,
         ng_per_ul_lwr95 = .lower,
         ng_per_ul_upr95 = .upper) %>%

  separate(ID, sep = ';;',
           into = c('plate_id',
                    'plate_column',
                    'plate_row',
                    'sample_type',
                    'dna_extract_tube_id'),
           convert = TRUE) %>%
  mutate(ng_per_ul_normspread = (ng_per_ul_upr95 - ng_per_ul_lwr95) / ng_per_ul_mean,
         sample_well_id = as.factor(sprintf("%s%02d", plate_row, plate_column))) %>%
  relocate(sample_well_id,
           .before = ng_per_ul_mean) %>%
  arrange(plate_id, plate_column, plate_row)

# quant_files_summarized <- predict(dna_model_freq,
#                                  newdata = distinct(dna_amounts,
#                                                     is_control,
#                                                     ID),
#                                  re.form = NULL,
#                                  type = 'link',
#                                  se.fit = TRUE) %>%
#   bind_cols(distinct(dna_amounts,
#                      is_control,
#                      ID),
#             .) %>%
#   mutate(lwr = fit - 1.96 * se.fit,
#          upr = fit + 1.96 * se.fit,
#          across(c(fit, lwr, upr), exp)) %>%
#   select(ID, is_control,
#          ng_per_ul_mean = fit,
#          ng_per_ul_lwr95 = lwr,
#          ng_per_ul_upr95 = upr) %>%
#   separate(ID, sep = ';;',
#            into = c('plate_id',
#                     'plate_column',
#                     'plate_row',
#                     'sample_type',
#                     'dna_extract_tube_id'),
#            convert = TRUE) %>%
#   mutate(ng_per_ul_normspread = (ng_per_ul_upr95 - ng_per_ul_lwr95) / ng_per_ul_mean,
#          sample_well_id = as.factor(sprintf("%s%02d", plate_row, plate_column))) %>%
#   relocate(sample_well_id,
#            .before = ng_per_ul_mean) %>%
#   arrange(plate_id, plate_column, plate_row)


mean_ng_per_ul <- add_epred_draws(newdata = data.frame(is_control = FALSE),
                dna_model_bayes,
                re.form = NA) %>%
  point_interval() %>%
  pull(.epred)


# predict(dna_model_freq, re.form = NA,
#         newdata = data.frame(is_control = FALSE),
#         se.fit = FALSE) %>%
#   exp

#### Flag Problem Samples ####
dna_quantity_interval <- lm(log(ng_per_ul_mean) ~ is_control, quant_files_summarized) %>%
  predict(newdata = data.frame(is_control = c(TRUE, FALSE)),
          interval = 'prediction', level = 0.95) %>%
  as_tibble() %>%
  mutate(is_control = c(TRUE, FALSE),
         across(where(is.numeric), exp)) %>%
  select(is_control, mean_dna = fit, lwr_limit = lwr, upr_limit = upr)

dna_variability_interval <- lm(log(ng_per_ul_normspread) ~ is_control,
                               data = quant_files_summarized) %>%
  predict(newdata = data.frame(is_control = c(TRUE, FALSE)),
          interval = 'prediction', level = 0.95) %>%
  as_tibble() %>%
  mutate(is_control = c(TRUE, FALSE),
         across(where(is.numeric), exp)) %>%
  select(is_control, lwr_limit = lwr, upr_limit = upr)


#### Calculate Dilutions to get to target DNA amounts ####
quant_files_flagged <- quant_files_summarized %>%

  mutate(ul_per_rxn = (target_dna * mean_ng_per_ul) / ng_per_ul_mean,
         ul_per_rxn = case_when(ul_per_rxn > maximum_low_volume & is_control ~ maximum_low_volume,
                                ul_per_rxn > maximum_low_volume & !is_control ~ maximum_low_volume,
                                ul_per_rxn < minimum_pipettable_volume ~ NA_real_,
                                TRUE ~ ul_per_rxn),
         ul_per_rxn = round(ul_per_rxn),

         dilution_factor = case_when(is.na(ul_per_rxn) ~ ng_per_ul_mean / mean_ng_per_ul,
                                     TRUE ~ 0),
         dilution_factor = round(dilution_factor),

         postDilution_ng_per_ul = if_else(dilution_factor > 0,
                                          ng_per_ul_mean / dilution_factor,
                                          NA_real_),
         postDilution_ul_per_rxn = (target_dna * mean_ng_per_ul) / postDilution_ng_per_ul,
         postDilution_ul_per_rxn = round(postDilution_ul_per_rxn),

         rxn_ng = if_else(is.na(ul_per_rxn),
                          postDilution_ng_per_ul * postDilution_ul_per_rxn,
                          ng_per_ul_mean * ul_per_rxn)) %>%

  left_join(select(dna_variability_interval, is_control, upr_limit),
            by = 'is_control') %>%

  mutate(flags = case_when(is_control & ng_per_ul_upr95 > mean_ng_per_ul ~ "Contaminated Control",
                           !is.na(postDilution_ng_per_ul) & ng_per_ul_normspread > upr_limit ~ 'Excess & Variable DNA',
                           !is.na(postDilution_ng_per_ul) ~ 'Excess DNA',
                           ng_per_ul_normspread > upr_limit ~ 'Variable DNA',
                           TRUE ~ 'Good Sample'),
         .after = sample_type) %>%

  select(plate_id, plate_column, plate_row,
         dna_extract_tube_id, sample_well_id,
         sample_type,
         ng_per_ul_mean, ng_per_ul_lwr95, ng_per_ul_upr95,
         flags,
         ul_per_rxn, dilution_factor,
         starts_with('postDilution'),
         rxn_ng)


#### Visualize quant results ####
(quant_plot <-
   quant_files_flagged  %>%
   ggplot (
     aes (
       y = fct_reorder (
         sample_well_id,
         desc(sample_well_id)
       ),
       x = ng_per_ul_mean,
       shape = sample_type,
       col = flags
     )
   )  +
   geom_vline(data = dna_quantity_interval,
              aes(xintercept = mean_dna,
                  linetype = is_control)) +
   geom_point(size = 1.5) +
   scale_colour_manual(
     values = c(
       'Good Sample' = 'Black',
       'Excess DNA' = '#F8766D',
       'Variable DNA' = '#00BFC4',
       "Contaminated Control" = 'purple',
       'Excess & Variable DNA' = 'orange'
     )
   ) +
   geom_errorbar(
     aes(xmin = ng_per_ul_lwr95,
         xmax = ng_per_ul_upr95
     )
   ) +
   scale_shape_manual(
     values = c(
       15,
       17,
       18,
       19
     )
   ) +
   scale_linetype_manual(values = c('TRUE' = 'dashed',
                                    'FALSE' = 'solid'),
                         labels = c('TRUE' = 'Control',
                                    'FALSE' = 'Sample')) +
   facet_wrap (~plate_id,
               ncol = 5,
               scales =  "free_y") +
   theme_bw() +
   labs (
     x = "DNA Concentration (ng/uL) \u00b1 95% CI",
     y = "Sample Well ID",
     color = "Sample Flag",
     shape = "Sample Type",
     linetype = 'Sample Type'
   ) +
   theme (axis.text.x = element_text (size = 7))
)



#### OUTPUT RESULTS ####
ggsave(
  quant_plot,
  file = path_quant_plot,
  width = 10,
  height = 5,
  units = "in",
  dpi = 330)

ggsave(
  quant_plot +
    scale_x_continuous(transform = scales::log10_trans(),
                       labels = scales::comma_format()),
  file = str_replace(path_quant_plot, 'assessment.png', 'assessment_log.png'),
  width = 10,
  height = 5,
  units = "in",
  dpi = 330)


write_csv(mutate(quant_files_flagged,
                 flags = na_if(flags, 'Good Sample')),
          file = path_quant_report_summarized,
          na = ''
)
