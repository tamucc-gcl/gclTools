#### quant_dna_from_filter.R ####
## Kevin Labrador
## 2025-01-27

# Set-up the working directory in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#### USER DEFINED ####
# Clear global environment.
rm(list = ls())

# Load Datasets
path_quant_plate_map <- "./rbd_edna-plate05_AccuClear_2024-08-27_plate-map_formatted.csv"

path_quant_raw_data <- "./rbd_edna-plate05_AccuClear_2024-08-27_raw-data.csv"

path_dna_plate_map <- "./rbd_edna-extraction_plate-map.xlsx"

sheet_metadata <- "Plate-map-tidy"

outdir_processed_data <- "./"
outdir_plots <- "./"

remove_these_standards <- c(0.015, 0.05, 0.15)

# scripts <- "../scripts/"



#### PROGRAM DEFINED VARIABLES

#### Load Libraries ####
pacman::p_load(
  janitor,
  readxl,
  cowplot,
  ggpubr,
  tidyverse
)

#### Source helper functions ####
source ("helper_functions.R")

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


#### Assess Standards ####

# Report the number of standards used and their concentrations
standard_count <-
  quant_plate_map %>%
  filter (grepl("standard", plate_id)) %>%
  nrow()
message (paste(standard_count, "standards detected"))

## Report the standards used
standard_samples <- quant_plate_map %>%
  filter (grepl("standard", plate_id)) %>%
  pull(sample_id)
message(paste("The standards used (in ng/ul) are:", paste(standard_samples, collapse = ", ")))


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

standard_filtered <-
  standard %>%
  remove_high_cv_replicates(
    cv_threshold = 15,
    group_var = "sample_id",
    response_var = "rfu")

standard_stats <-
  standard_filtered %>%
  
  # calculate the mean, sd, and coef.var of the standards
  group_by (sample_id) %>%
  summarize (
    n_replicates = n(),
    rfu_mean = mean(rfu),
    rfu_sd = sd(rfu),
    rfu_cv_pct = rfu_sd / rfu_mean*100,
    dna_per_well = mean(dna_per_well),
    dna_per_ul = mean(dna_per_ul),
  ) %>%
  ungroup() %>%
  mutate(
    include_in_model = 
      case_when (
        dna_per_well %in% remove_these_standards ~ F,
        T ~ T)
  )


data_to_model <-
  standard_stats %>%
  filter (include_in_model)


#### Assess Linear Model Fit ####
fit <- lm (
  dna_per_well ~ rfu_mean,
  data = data_to_model)

summary(fit)

pred_data <- data.frame(
  rfu_mean = seq(
    min(data_to_model$rfu_mean),
    max(data_to_model$rfu_mean),
    length.out = 100)
) %>%
  mutate (
    log_rfu = log10(rfu_mean),
    dna_pred = predict(fit, newdata = list(rfu_mean = rfu_mean))  # No need to exponentiate
  )

## Generate plot for visualization
subtitle <- formula(fit) %>% as.character
aic <- AIC(fit) %>% round(4)
bic <- BIC(fit) %>% round(4)
rse <- summary(fit)$sigma %>% round(4)

(plot_linear <-
    
    ggplot (standard_stats,
            aes (x = log10(rfu_mean),
                 y = log10(dna_per_well))) +
    
    geom_errorbarh(
      aes (xmin = log10(rfu_mean - rfu_sd),
           xmax = log10(rfu_mean + rfu_sd),
           color = include_in_model),
      height = 0.25) +
    
    geom_point(
      aes(pch = "Standard",
          color = include_in_model),
      size = 2) +
    
    geom_text(
      aes(label = dna_per_well),
      size = 3.5,
      nudge_x = 0.01,
      nudge_y = 0.1) +
    
    
    geom_line (
      data = pred_data,
      aes(x = log_rfu,
          y = log10(dna_pred)),
      color = "red",
      lwd = 1) +
    
    labs (
      title = "Linear Model",
      subtitle = paste(
        "Model:",
        subtitle [2],
        subtitle [1],
        subtitle [3],
        "\nRSE =",
        rse,
        "; AIC =",
        aic,
        "; BIC =",
        bic
      ),
      x = expression(log[10]*"(RFU)"),
      y = expression(log[10]*"(DNA concentration; ng/well)")) +
    theme_bw()
)

# Save model parameters
fit_linear <- fit

#### Assess Power Model Fit ####
fit <- lm (
  log10(dna_per_well) ~ log10(rfu_mean),
  data = data_to_model)

summary(fit)

pred_data <- data.frame(
  rfu_mean = seq(
    min(data_to_model$rfu_mean),
    max(data_to_model$rfu_mean),
    length.out = 100)
) %>%
  mutate (
    log_rfu = log10(rfu_mean),
    log_dna_pred = predict (fit,
                            newdata = list(rfu_mean = rfu_mean)),
    dna_pred = 10^log_dna_pred
  )

## Generate plot for visualization
subtitle <- formula(fit) %>% as.character
aic <- AIC(fit) %>% round(4)
bic <- BIC(fit) %>% round(4)
rse <- summary(fit)$sigma %>% round(4)

(plot_power <-
    
    ggplot (standard_stats,
            aes (x = log10(rfu_mean),
                 y = log10(dna_per_well))) +
    
    geom_errorbarh(
      aes (xmin = log10(rfu_mean - rfu_sd),
           xmax = log10(rfu_mean + rfu_sd),
           color = include_in_model),
      height = 0.25) +
    
    geom_point(
      aes(pch = "Standard",
          color = include_in_model),
      size = 2) +
    
    geom_text(
      aes(label = dna_per_well),
      size = 3.5,
      nudge_x = 0.01,
      nudge_y = 0.1) +
    
    
    geom_line (
      data = pred_data,
      aes(x = log_rfu,
          y = log_dna_pred),
      color = "red",
      lwd = 1) +
    
    labs (
      title = "Power Model",
      subtitle = paste(
        "Model:",
        subtitle [2],
        subtitle [1],
        subtitle [3],
        "\nRSE =",
        rse,
        "; AIC =",
        aic,
        "; BIC =",
        bic
      ),
      x = expression(log[10]*"(RFU)"),
      y = expression(log[10]*"(DNA concentration; ng/well)")) +
    theme_bw()
)

# Save model parameters
fit_power <- fit

#### Assess Polynomial Model Fit ####
fit <- lm (
  log10(dna_per_well) ~ poly (log10(rfu_mean), degree = 2),
  data = data_to_model)

summary(fit)

pred_data <- data.frame(
  rfu_mean = seq(
    min(data_to_model$rfu_mean),
    max(data_to_model$rfu_mean),
    length.out = 100)
) %>%
  mutate (
    log_rfu = log10(rfu_mean),
    log_dna_pred = predict (fit,
                            newdata = list(rfu_mean = rfu_mean)),
    dna_pred = 10^log_dna_pred
  )

## Generate plot for visualization
subtitle <- formula(fit) %>% as.character
aic <- AIC(fit) %>% round(4)
bic <- BIC(fit) %>% round(4)
rse <- summary(fit)$sigma %>% round(4)

(plot_polynomial <-
    
    ggplot (standard_stats,
            aes (x = log10(rfu_mean),
                 y = log10(dna_per_well))) +
    
    geom_errorbarh(
      aes (xmin = log10(rfu_mean - rfu_sd),
           xmax = log10(rfu_mean + rfu_sd),
           color = include_in_model),
      height = 0.25) +
    
    geom_point(
      aes(pch = "Standard",
          color = include_in_model),
      size = 2) +
    
    geom_text(
      aes(label = dna_per_well),
      size = 3.5,
      nudge_x = 0.01,
      nudge_y = 0.1) +
    
    
    geom_line (
      data = pred_data,
      aes(x = log_rfu,
          y = log_dna_pred),
      color = "red",
      lwd = 1) +
    
    labs (
      title = "Polynomial Model",
      subtitle = paste(
        "Model:",
        subtitle [2],
        subtitle [1],
        subtitle [3],
        "\nRSE =",
        rse,
        "; AIC =",
        aic,
        "; BIC =",
        bic
      ),
      x = expression(log[10]*"(RFU)"),
      y = expression(log[10]*"(DNA concentration; ng/well)")) +
    theme_bw()
)

# Save model parameters
fit_polynomial <- fit

#### Assess Logistic Model Fit ####

# Assign nls starting parameters

## asymptotic length
L <- max(log10(standard_stats$dna_per_well))

## growth rate coefficient
k <- 1.3

## midpoint or inflection point
x0 <- mean(log10(standard_stats$dna_per_well))


fit <- nls(
  log10(dna_per_well) ~ L / (1 + exp(-k * (log10(rfu_mean) - x0))),
  data = data_to_model,
  algorithm = "port",
  start = list (
    L = L, # Asymptote
    k = k, # Growth rate coefficient
    x0 = x0 # Inflection point
  )
)

summary(fit)

pred_data <- data.frame(
  rfu_mean = seq(
    min(data_to_model$rfu_mean),
    max(data_to_model$rfu_mean),
    length.out = 100)
) %>%
  mutate (
    log_rfu = log10(rfu_mean),
    log_dna_pred = predict (fit,
                            newdata = list(rfu_mean = rfu_mean)),
    dna_pred = 10^log_dna_pred
  )

## Generate plot for visualization
subtitle <- formula(fit) %>% as.character
aic <- AIC(fit) %>% round(4)
bic <- BIC(fit) %>% round(4)
rse <- summary(fit)$sigma %>% round(4)

(plot_logistic <-
    
    ggplot (standard_stats,
            aes (x = log10(rfu_mean),
                 y = log10(dna_per_well))) +
    
    geom_errorbarh(
      aes (xmin = log10(rfu_mean - rfu_sd),
           xmax = log10(rfu_mean + rfu_sd),
           color = include_in_model),
      height = 0.25) +
    
    geom_point(
      aes(pch = "Standard",
          color = include_in_model),
      size = 2) +
    
    geom_text(
      aes(label = dna_per_well),
      size = 3.5,
      nudge_x = 0.01,
      nudge_y = 0.1) +
    
    
    geom_line (
      data = pred_data,
      aes(x = log_rfu,
          y = log_dna_pred),
      color = "red",
      lwd = 1) +
    
    labs (
      title = "Logistic Model",
      subtitle = paste(
        "Model:",
        subtitle [2],
        subtitle [1],
        subtitle [3],
        "\nRSE =",
        rse,
        "; AIC =",
        aic,
        "; BIC =",
        bic
      ),
      x = expression(log[10]*"(RFU)"),
      y = expression(log[10]*"(DNA concentration; ng/well)")) +
    theme_bw()
)

# Save model parameters
fit_logistic <- fit


## Output plots
(plot_model_comparison <-
    plot_grid(
      plot_linear + guides(col = "none"),
      plot_power + guides(col = "none"),
      plot_polynomial + guides(col = "none"),
      plot_logistic + guides(col = "none"),
      nrow = 2,
      labels = "AUTO"
    ))



#### Sample Quant Prediction  ####
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
  quant_data_for_predictions %>%
  rename (rfu_mean = rfu) %>%  # Rename rfu to rfu_mean to coincide with model's formula
  mutate (
    dna_per_well =
      predict(fit_linear,
              newdata = .),
    ng_per_ul = dna_per_well/sample_volume,
    model = "linear"
  )

predicted_sample_concentration_power <-
  quant_data_for_predictions %>%
  rename (rfu_mean = rfu) %>%  # Rename rfu to rfu_mean to coincide with model's formula
  mutate (
    dna_per_well =
      10^predict(fit_power,
                 newdata = .),
    ng_per_ul = dna_per_well/sample_volume,
    model = "power"
  )

predicted_sample_concentration_polynomial <-
  quant_data_for_predictions %>%
  rename (rfu_mean = rfu) %>%  # Rename rfu to rfu_mean to coincide with model's formula
  mutate (
    dna_per_well =
      10^predict(fit_polynomial,
                 newdata = .),
    ng_per_ul = dna_per_well/sample_volume,
    model = "polynomial"
  )

predicted_sample_concentration_logistic <-
  quant_data_for_predictions %>%
  rename (rfu_mean = rfu) %>%  # Rename rfu to rfu_mean to coincide with model's formula
  mutate (
    dna_per_well =
      10^predict(fit_logistic,
                 newdata = .),
    ng_per_ul = dna_per_well/sample_volume,
    model = "logistic"
  )

predicted_sample_concentration <-
  rbind (
    predicted_sample_concentration_linear,
    predicted_sample_concentration_power,
    predicted_sample_concentration_polynomial,
    predicted_sample_concentration_logistic
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
plot_linear_with_samples <-
  plot_linear +
  geom_point(
    data = predicted_sample_concentration %>% filter (model == "linear"),
    aes (x = log10(rfu),
         pch = "Sample"),
    alpha = 0.25,
    size = 2.5) +
  labs (col = "Included in model") +
  scale_shape_manual(
    name = "Classification",
    values = c(
      "Standard" = 15,
      "Sample"= 16),
    labels = c(
      "Standard" = "Standard (ng/well)",
      "Sample" = "Sample")
  ) +
  geom_vline(
    xintercept = log10(min_rfu), 
    lty = 2
  ) + 
  annotate (
    "text",
    x = log10(min_rfu)-0.05,
    y = log10(max_rfu)/2.25,
    label = "below lod",
    angle = 90,
    size = 3) + 
  geom_vline(
    xintercept = log10(max_rfu), 
    lty = 2
  ) + 
  annotate (
    "text",
    x = log10(max_rfu)+0.05,
    y = 0,
    label = "above lod",
    angle = 90,
    size = 3) 

plot_linear_with_samples

plot_power_with_samples <-
  plot_power +
  geom_point(
    data = predicted_sample_concentration %>% filter (model == "power"),
    aes (x = log10(rfu),
         pch = "Sample"),
    alpha = 0.25,
    size = 2.5) +
  labs (col = "Included in model") +
  scale_shape_manual(
    name = "Classification",
    values = c(
      "Standard" = 15,
      "Sample"= 16),
    labels = c(
      "Standard" = "Standard (ng/well)",
      "Sample" = "Sample")
  ) +
  geom_vline(
    xintercept = log10(min_rfu), 
    lty = 2
  ) + 
  annotate (
    "text",
    x = log10(min_rfu)-0.05,
    y = log10(max_rfu)/2.25,
    label = "below lod",
    angle = 90,
    size = 3) + 
  geom_vline(
    xintercept = log10(max_rfu), 
    lty = 2
  ) + 
  annotate (
    "text",
    x = log10(max_rfu)+0.05,
    y = 0,
    label = "above lod",
    angle = 90,
    size = 3) 
plot_power_with_samples

plot_polynomial_with_samples <-
  plot_polynomial +
  geom_point(
    data = predicted_sample_concentration %>% filter (model == "polynomial"),
    aes (x = log10(rfu),
         pch = "Sample"),
    alpha = 0.25,
    size = 2.5) +
  labs (col = "Included in model") +
  scale_shape_manual(
    name = "Classification",
    values = c(
      "Standard" = 15,
      "Sample"= 16),
    labels = c(
      "Standard" = "Standard (ng/well)",
      "Sample" = "Sample")
  ) +
  geom_vline(
    xintercept = log10(min_rfu), 
    lty = 2
  ) + 
  annotate (
    "text",
    x = log10(min_rfu)-0.05,
    y = log10(max_rfu)/2.25,
    label = "below lod",
    angle = 90,
    size = 3) + 
  geom_vline(
    xintercept = log10(max_rfu), 
    lty = 2
  ) + 
  annotate (
    "text",
    x = log10(max_rfu)+0.05,
    y = 0,
    label = "above lod",
    angle = 90,
    size = 3) 
plot_polynomial_with_samples

plot_logistic_with_samples <-
  plot_logistic +
  geom_point(
    data = predicted_sample_concentration %>% filter (model == "logistic"),
    aes (x = log10(rfu),
         pch = "Sample"),
    alpha = 0.25,
    size = 2.5) +
  labs (col = "Included in model") +
  scale_shape_manual(
    name = "Classification",
    values = c(
      "Standard" = 15,
      "Sample"= 16),
    labels = c(
      "Standard" = "Standard (ng/well)",
      "Sample" = "Sample")
  ) + 
  geom_vline(
    xintercept = log10(min_rfu), 
    lty = 2
  ) + 
  annotate (
    "text",
    x = log10(min_rfu)-0.05,
    y = log10(max_rfu)/2.25,
    label = "below lod",
    angle = 90,
    size = 3) + 
  geom_vline(
    xintercept = log10(max_rfu), 
    lty = 2
  ) + 
  annotate (
    "text",
    x = log10(max_rfu)+0.05,
    y = 0,
    label = "above lod",
    angle = 90,
    size = 3) 
plot_logistic_with_samples

#### Output results ####
(plot_with_samples <-
   ggarrange(
     plot_linear_with_samples,
     plot_power_with_samples,
     plot_polynomial_with_samples,
     plot_logistic_with_samples,
     ncol = 2,
     nrow = 2,
     labels = "AUTO",
     common.legend = T,
     legend = "right"
   ))

basename <- path_quant_raw_data %>%
  basename() %>%
  gsub ("_raw-data.csv", "", .)

ggsave(
  plot_with_samples,
  file = paste0(outdir_plots,
                basename,
                "_quant_model_comparison_replicates.jpg"),
  width = 12,
  height = 10,
  units = "in",
  dpi = 330)


#### Select Best model and summarize dataset ####

#### 1. Gather the model metrics
aic_linear   <- AIC(fit_linear)
bic_linear   <- BIC(fit_linear)
rse_linear   <- summary(fit_linear)$sigma

aic_power <- AIC(fit_power)
bic_power <- BIC(fit_power)
rse_power <- summary(fit_power)$sigma  # residual standard error

aic_poly  <- AIC(fit_polynomial)
bic_poly  <- BIC(fit_polynomial)
rse_poly  <- summary(fit_polynomial)$sigma

aic_log   <- AIC(fit_logistic)
bic_log   <- BIC(fit_logistic)
rse_log   <- summary(fit_logistic)$sigma

#### 2. Create a summary data frame for each model
metrics_df <- tibble::tibble(
  model = c("linear", "power", "polynomial", "logistic"),
  AIC   = c(aic_linear, aic_power, aic_poly, aic_log),
  BIC   = c(bic_linear, bic_power, bic_poly, bic_log),
  RSE   = c(rse_linear, rse_power, rse_poly, rse_log)
)

#### 3. Rank each metric
# Lower = better, so we rank with ties.method = "first" (or "min") so the best (lowest) gets rank 1
metrics_df <- metrics_df %>%
  dplyr::mutate(
    AIC_rank = dplyr::dense_rank(AIC),
    BIC_rank = dplyr::dense_rank(BIC),
    RSE_rank = dplyr::dense_rank(RSE),
    total_rank = AIC_rank + BIC_rank + RSE_rank
  )

#### 4. Choose the model with the best (lowest) total_rank

# Make a factor to reflect final preference if total_rank, BIC, AIC, and RSE are all tied
model_priority <- c("power", "polynomial", "logistic", "linear")

metrics_df <-
  metrics_df %>%
  mutate(
    model_factor = factor(model, levels = model_priority)
  ) %>%
  arrange(
    total_rank,  # 1) primary
    BIC,         # 2) tiebreak
    AIC,         # 3) tiebreak
    RSE,         # 4) tiebreak
    model_factor # 5) final fallback
  )

# Examine the summary table (optional):
metrics_df


# The first row is our best model
best_model <- metrics_df$model[1]

message("The ", toupper(best_model), " model was selected based on following prioritization total_rank > BIC > AIC > RSE > Power > Polynomial > Logistic > Linear")


#### 5. Filter the predictions to keep only the best model’s results (Section moved to `summarize_quant_results.R`)
final_predicted_samples <-
  predicted_sample_concentration %>%
  dplyr::filter(model == best_model)


# Export data frame
write_csv(
  final_predicted_samples,
  file = paste0(outdir_processed_data,
                basename,
                "_quant_report_replicates.csv")
)


#### 5. Summarize quant results (Section moved to `summarize_quant_results.R`)
# final_predicted_samples_summarized <-
#   predicted_sample_concentration %>%
#   dplyr::filter(model == best_model) %>%
#   group_by(
#     plate_id,
#     plate_column,
#     plate_row,
#     sample_type,
#     dna_extract_tube_id,
#     model,
#   ) %>%
#   summarize(
#     # Number of replicates
#     n_reps = n(),
#     
#     # Means (na.rm = TRUE ensures that NA values won't break the calculation)
#     sample_volume_mean = mean(sample_volume, na.rm = TRUE),
#     rfu_mean           = mean(rfu, na.rm = TRUE),
#     dna_per_well_mean  = mean(dna_per_well, na.rm = TRUE),
#     ng_per_ul_mean     = mean(ng_per_ul, na.rm = TRUE),
#     
#     # Standard deviation of ng_per_ul
#     ng_per_ul_sd = sd(ng_per_ul, na.rm = TRUE),
#     
#     # CV of ng_per_ul (SD / Mean * 100)
#     ng_per_ul_cv = (sd(ng_per_ul, na.rm = TRUE) / mean(ng_per_ul, na.rm = TRUE)) * 100,
#     
#     # 95% CI for ng_per_ul (using a t-based approach, suitable for small sample sizes)
#     ng_per_ul_ci_95 = qt(0.975, df = n() - 1) * (sd(ng_per_ul, na.rm = TRUE) / sqrt(n())),
#     
#     # % of replicates that are "within_lod"
#     pct_reps_within_lod = (sum(quant_status == "within_lod", na.rm = TRUE) / n()) * 100
#   ) %>%
#   ungroup()
# 
# 
# #### OUTPUT RESULTS ####
# 
# # Export data frame
# final_predicted_samples_summarized %>%
#   select(
#     -contains("_cv"),
#     -contains("_rfu"),
#     -contains("per_well"),
#     -contains("_sd")
#   ) %>%
#   arrange (dna_extract_tube_id) %>% 
#   write_csv(
#     file = paste0(outdir_processed_data,
#                   basename,
#                   "_quant_report_summarized.csv")
#   )
# 

#### Visualize Summary ####

# Suppose you already created final_summarized via something like:
# final_summarized <- final_predicted_samples %>%
#   group_by(dna_extract_tube_id) %>%
#   summarize(
#     n_reps              = n(),
#     sample_volume_mean  = mean(sample_volume, na.rm = TRUE),
#     rfu_mean            = mean(rfu, na.rm = TRUE),
#     dna_per_well_mean   = mean(dna_per_well, na.rm = TRUE),
#     ng_per_ul_mean      = mean(ng_per_ul, na.rm = TRUE),
#     ng_per_ul_sd        = sd(ng_per_ul, na.rm = TRUE),
#     ng_per_ul_cv        = (sd(ng_per_ul, na.rm = TRUE) /
#                           mean(ng_per_ul, na.rm = TRUE)) * 100,
#     ng_per_ul_ci_95     = qt(0.975, df = n() - 1) *
#                          (sd(ng_per_ul, na.rm = TRUE) / sqrt(n())),
#     pct_reps_within_lod = (sum(quant_status == "within_lod", na.rm = TRUE) / n()) * 100
#   ) %>%
#   ungroup()
# 
# # 1) Identify "suspicious" or "high-variation" samples to highlight.
# #    For example, if CV > 20% or not all replicates are within LOD:
# final_summarized_flagged <- final_summarized %>%
#   mutate(
#     flag = case_when(
#       ng_per_ul_cv > 20                 ~ "check_high_cv",         # High variation
#       pct_reps_within_lod < 100         ~ "check_partial_outside", # Some replicates outside LOD
#       TRUE                              ~ "ok"
#     )
#   )
# 
# # 2) Plot mean ng/ul with 95% CI error bars, sorted by mean
# #    We'll color or shape points by the 'flag' status.
# ggplot(final_summarized_flagged,
#        aes(x = reorder(dna_extract_tube_id, ng_per_ul_mean),
#            y = ng_per_ul_mean,
#            color = flag)) +
# 
#   # Point + error bar (95% CI)
#   geom_pointrange(
#     aes(ymin = ng_per_ul_mean - ng_per_ul_ci_95,
#         ymax = ng_per_ul_mean + ng_per_ul_ci_95),
#     size = 0.5
#   ) +
# 
#   # Flip coordinates so sample labels go down the y-axis
#   coord_flip() +
# 
#   labs(
#     x = "DNA Extract Tube ID",
#     y = "Mean DNA Concentration (ng/µl)",
#     color = "Sample Status",
#     title = "Mean DNA Conc. ± 95% CI",
#     subtitle = "Flagged if CV > 20% or partial replicates outside LOD"
#   ) +
#   theme_bw(base_size = 12)
