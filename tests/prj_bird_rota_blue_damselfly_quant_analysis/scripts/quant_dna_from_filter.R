#### quant_dna_from_filter.R ####
## Kevin Labrador
## 2025-01-27

####-----INTRODUCTION-----####


####-----INITIALIZE-----####

#### HOUSEKEEPING ####
# Clear global environment.
rm(list = ls())

# Set-up the working directory in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Load Libraries ####
pacman::p_load(
  janitor,
  readxl,
  cowplot,
  ggpubr,
  outliers,
  tidyverse
)

#### Source helper functions ####
source ("helper_functions.R")

#### USER DEFINED VARIABLES ####
# Assign input file paths
path_quant_plate_map <-
  paste0(
    "../data_processed/",
    "rbd_edna-plate01_AccuClear_2024-05-14_plate-map_formatted.csv"
  )

path_quant_raw_data <-
  paste0(
    "../data_raw/",
    "rbd_edna-plate01_AccuClear_2024-05-14_raw-data.csv"
  )

path_dna_plate_map <-
  paste0(
    "../data_raw/",
    "rbd_edna-extraction_plate-map.xlsx"
  )

sheet_metadata <- "Plate-map-tidy"

# Assign output file paths
basename <- path_quant_raw_data %>%
  basename() %>%
  gsub (".csv", "", .)

path_quant_report <-
  paste0(
    "../data_processed/",
    basename,
    "_quant_report.csv"
  )

path_quant_model_comparison <-
  paste0(
    "../results/",
    basename,
    "_quant_model_comparison.jpg"
  )


# Set standards that needs to be removed
remove_these_standards <- c(0)

# Assign the quant kit used
quant_kit <- "accuclear" # Choices: accublue-nextgen or accuclear


#### Load Datasets ####
# Quant raw data
quant_raw_data <-
  read.csv (path_quant_raw_data) %>%
  clean_names() %>%
  select (wells, value) %>%
  mutate(across(where(is.character), str_to_lower))

# Quant plate map
quant_plate_map <-
  read.csv (path_quant_plate_map) %>%
  clean_names() %>%
  mutate (
    wells = paste0(quant_row, quant_column),
    sample_id = trimws(sample_id),
    across(where(is.character), str_to_lower),
    # If sample_id starts with digits followed by letters (e.g. "12b"), swap them to become "b12"
    sample_id = str_replace(sample_id, "^([0-9]+)([a-z]+)$", "\\2\\1")
  )

# DNA plate map
dna_plate_map <-
  read_xlsx(
    path_dna_plate_map,
    sheet = sheet_metadata
  ) %>%
  clean_names() %>%
  mutate (wells = paste0(plate_row, plate_col)) %>%
  mutate(across(where(is.character), str_to_lower))


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



#### Assess Standards ####
standard_report <-
  report_standards(
    quant_plate_map = quant_plate_map,
    quant_kit = quant_kit
  )

#### Standard Quant Regression ####
standard <-
  quant_data$standard %>%
  mutate (
    sample_id = as.numeric(sample_id),
    dna_per_ul = as.numeric(sample_id),
    dna_per_well = dna_per_ul * sample_volume
  ) %>%

  # Remove 0 standard
  filter (sample_id != 0)

standard_stats <-
  standard %>%

  # calculate the mean, sd, and coef.var of the standards
  group_by (sample_id) %>%
  summarize (
    n_reps = n(),

    # Means
    rfu_mean = mean(rfu),
    dna_per_well = mean(dna_per_well),
    dna_per_ul = mean(dna_per_ul),

    # SD
    rfu_sd = sd(rfu),

    # %CV
    rfu_cv_pct = rfu_sd / rfu_mean*100,

    # CI
    rfu_ci_95 =
      ifelse(n_reps > 1,
             qt(0.975, df = n_reps - 1) * (rfu_sd / sqrt(n_reps)),
             NA_real_),  # Avoids errors when n_reps == 1
    .groups = "drop"
  ) %>%
  mutate(
    include_in_model =
      case_when (
        dna_per_well %in% remove_these_standards ~ F,
        T ~ T)
  )


#### Assess Model Fit ####
fit_linear <-
  fit_and_plot(
    x_var = "rfu_mean",
    y_var = "dna_per_well",
    data = standard_stats,
    model_type = "Linear",
    degree = NULL
  )

fit_power <-
  fit_and_plot(
    x_var = "rfu_mean",
    y_var = "dna_per_well",
    data = standard_stats,
    model_type = "Power",
    degree = NULL
  )

fit_polynomial <-
  fit_and_plot(
    x_var = "rfu_mean",
    y_var = "dna_per_well",
    data = standard_stats,
    model_type = "Polynomial",
    degree = 2
  )

#### Jacknife Fit & Identify Standards to Remove #### - JDS
jacknife_standards <- mutate(standard_stats,
              sample_id = row_number()) %>%
  filter(include_in_model) %>%
  
  expand_grid(tibble(model = c('Linear', 'Power', 'Polynomial'),
                     full_model = list(fit_linear$fit, fit_power$fit, fit_polynomial$fit))) %>%
  nest(removed_standard = -c(sample_id, model, full_model)) %>%
  rowwise %>%
  
  #Make dataset with one standard removed
  mutate(other_standards = list(anti_join(standard_stats,
                                          removed_standard,
                                          by = colnames(removed_standard)))) %>%

  #Fit model without standard
  mutate(jacknife_model = list(fit_and_plot(x_var = "rfu_mean",
                                            y_var = "dna_per_well",
                                            data = other_standards,
                                            model_type = model,
                                            degree = 2) %>%
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
  select(sample_id, model, rel_improvement)


#### Identify Outlier Standards #### - JDS
# x <- filter(jacknife_standards, model == 'Power') %>% pull(rel_improvement) %>% unname
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

# jacknife_standards %>%
#   mutate(is_outlier = identify_outliers(rel_improvement),
#          .by = 'model') %>%
#   ggplot(aes(x = sample_id, y = rel_improvement, colour = is_outlier)) +
#   geom_point() +
#   facet_wrap(~model, scales = 'free_y')


#### Refit models excluding poorly fit Samples #### - JDS
standards_without_outliers <- jacknife_standards %>%
  mutate(is_outlier = identify_outliers(rel_improvement),
         .by = 'model', .keep = 'unused') %>%
  filter(is_outlier) %>%
  select(-is_outlier) %>%
  summarise(outlier_standards = list(c(sample_id)),
            .by = model) %>%
  rowwise %>%
  mutate(new_standards = list(mutate(standard_stats, 
                                     sample_id = row_number(),
                                     include_in_model = !sample_id %in% outlier_standards & include_in_model)),
         .keep = 'unused') %>%
  ungroup

fit_linear <-
  fit_and_plot(
    x_var = "rfu_mean",
    y_var = "dna_per_well",
    data = standards_without_outliers$new_standards[[which(standards_without_outliers$model == 'Linear')]],
    model_type = "Linear",
    degree = NULL
  )

fit_power <-
  fit_and_plot(
    x_var = "rfu_mean",
    y_var = "dna_per_well",
    data = standards_without_outliers$new_standards[[which(standards_without_outliers$model == 'Power')]],
    model_type = "Power",
    degree = NULL
  )

fit_polynomial <-
  fit_and_plot(
    x_var = "rfu_mean",
    y_var = "dna_per_well",
    data = standards_without_outliers$new_standards[[which(standards_without_outliers$model == 'Polynomial')]],
    model_type = "Polynomial",
    degree = 2
  )



## Output plots
(plot_model_comparison <-
    plot_grid(
      fit_linear$plot + guides(col = "none"),
      fit_power$plot + guides(col = "none"),
      fit_polynomial$plot + guides(col = "none"),
      nrow = 3,
      labels = "AUTO"
    ))

#### Sample Quant Prediction ####

# Get the min and max values of the standards that were included in the model
min_rfu <-
  standard_stats %>%
  filter (include_in_model == T) %>%
  pull(rfu_mean) %>%
  min()

max_rfu <-
  standard_stats %>%
  filter (include_in_model == T) %>%
  pull(rfu_mean) %>%
  max()

# Designate samples that are beyond the quant's limit of detection (LOD)
quant_data_for_predictions <-
  quant_data$sample %>%
  mutate (
    quant_status =
      as.factor (
        case_when (
          rfu < min_rfu ~ "below_lod",
          rfu > max_rfu ~ "above_lod",
          T ~ "within_lod"
        )
      )
  )

# Predict sample dataset
predicted_sample_concentration_linear <-
  predict_sample_concentration(
    data = quant_data_for_predictions,
    model = fit_linear$fit,
    model_type =  "linear"
  )


predicted_sample_concentration_power <-
  predict_sample_concentration(
    data = quant_data_for_predictions,
    model = fit_power$fit,
    model_type =  "power"
  )

predicted_sample_concentration_polynomial <-
  predict_sample_concentration(
    data = quant_data_for_predictions,
    model = fit_polynomial$fit,
    model_type =  "polynomial"
  )



predicted_sample_concentration <-
  rbind (
    predicted_sample_concentration_linear,
    predicted_sample_concentration_power,
    predicted_sample_concentration_polynomial
  ) %>%
  select (-wells) %>%
  rename (
    # dna_extract_tube_id = sample_id,
    rfu = rfu_mean
  ) %>%
  left_join (
    dna_plate_map,
    by = c(
      "sample_id" = "wells",
      "plate_id" = "plate_id"
    )
  ) %>%
  select (
    plate_id,
    # dna_extract_tube_id,
    # sample_id,
    replicate,
    # preservative,
    plate_row,
    plate_col,
    quant_row,
    quant_column,
    sample_volume,
    model,
    rfu,
    dna_per_well,
    ng_per_ul,
    quant_status,
    dna_extract_tube_id,
    sample_type
  ) %>%
  rename (plate_column = plate_col) %>%
  mutate (plate_row = toupper(plate_row),
          quant_row = toupper(quant_row)
  )

#### Visualize dataset ####

## Linear
(plot_linear_with_samples <-
   plot_regression_with_samples(
     base_plot = fit_linear$plot,
     predicted_sample_concentration = predicted_sample_concentration,
     model_type = "linear",
     min_rfu = min_rfu,
     max_rfu = max_rfu
   )
)

## Power
(plot_power_with_samples <-
    plot_regression_with_samples(
      base_plot = fit_power$plot,
      predicted_sample_concentration = predicted_sample_concentration,
      model_type = "power",
      min_rfu = min_rfu,
      max_rfu = max_rfu
    )
)

## Polynomial
(plot_polynomial_with_samples <-
    plot_regression_with_samples(
      base_plot = fit_polynomial$plot,
      predicted_sample_concentration = predicted_sample_concentration,
      model_type = "polynomial",
      min_rfu = min_rfu,
      max_rfu = max_rfu
    )
)


#### Merge Plots ####
(plot_with_samples <-
   ggarrange(
     plot_linear_with_samples,
     plot_power_with_samples,
     plot_polynomial_with_samples,
     nrow = 3,
     labels = "AUTO",
     common.legend = T,
     legend = "right"
   ))


#### Select best model ####
model_rank_df <-
  rank_models(
    linear = fit_linear$fit,
    power = fit_power$fit,
    polynomial = fit_polynomial$fit
  )
print(model_rank_df)

# The first row is our best model
best_model <- model_rank_df$model[1]
message("The ", toupper(best_model), " model was selected based on following prioritization total_rank > BIC > AIC > RSE > Power > Polynomial > Logistic > Linear")


#### Filter the predictions to keep only the best model’s results ####
final_predicted_samples <-
  predicted_sample_concentration %>%
  dplyr::filter(model == best_model)



#### Export output tables and figures ####
ggsave(
  plot_with_samples,
  file = path_quant_model_comparison,
  width = 7.5,
  height = 10,
  units = "in",
  dpi = 330)


write_csv(
  final_predicted_samples,
  file = path_quant_report
)

