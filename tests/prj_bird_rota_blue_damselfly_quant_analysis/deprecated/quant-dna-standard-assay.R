#### quant-dna ####
## Kevin Labrador
## 2024-02-22

#### INTRODUCTION ####
# This script uses R to analyze the flourescence data exported from the SpectraMax software.
# This is a version of the quant-dna script where we tried to perform an assay exclusively for the standards.

## 1. raw plate reader data: this is the file exported from the SpectraMax software. It contains the well and fluorescence information.

## 2. plate map: this contains the information on how the samples were laid out in the 364-well plate.

# These files should be located in the '../data/quants' directory.

#### INITIALIZE ####

#### Housekeeping ####
# Clear global environment.
rm(list = ls())

# Set-up the working directory in the source file location: 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

wd <- getwd() # Do this after the working directory was set to source file location

#### Load Libraries ####
require (pacman)
pacman::p_load(janitor,
               ggpubr,
               readxl,
               lme4,
               mixtox,
               REAT,               
               cowplot,
               rstatix,
               ggpmisc,
               ggbeeswarm,
               nlme,
               tidyverse

)

#### USER-DEFINED VARIABLES ####
# Assign directories to objects in the global environment
indir <- "../data/quants/"
outdir <- "../output/quants/"
scripts <- "../scripts/"
# All user-defined variables should be here.

#### Load Dataset ####
# Check files to load
list.files(indir)

# Load raw quant data

raw_quant_data_5min_normal <- 
  paste0(indir, "rme_std-assay_AccuBlueNextGen_5min-normal-plate-orientation_2024-02-22_raw-data.xlsx") %>% 
  read_excel() %>% 
  clean_names() %>% 
  # Remove unwanted columns
  select (-c(sample, r)) %>% 
  mutate (setup = "normal-5min",
          plate_orientation = "normal",
          incubation_time = 5)


raw_quant_data_5min_reverse <- 
  paste0(indir, "rme_std-assay_AccuBlueNextGen_5min-reverse-plate-orientation_2024-02-22_raw-data.xlsx") %>% 
  read_excel() %>% 
  clean_names() %>% 
  # Remove unwanted columns
  select (-c(sample, r)) %>% 
  mutate (setup = "reverse-5min",
          plate_orientation = "reverse",
          incubation_time = 5)

raw_quant_data_20min_normal <- 
  paste0(indir, "rme_std-assay_AccuBlueNextGen_20min-normal-plate-orientation_2024-02-22_raw-data.xlsx") %>% 
  read_excel() %>% 
  clean_names() %>% 
  # Remove unwanted columns
  select (-c(sample, r)) %>% 
  mutate (setup = "normal-20min",
          plate_orientation = "normal",
          incubation_time = 20)


raw_quant_data_20min_reverse <- 
  paste0(indir, "rme_std-assay_AccuBlueNextGen_20min-reverse-plate-orientation_2024-02-22_raw-data.xlsx") %>% 
  read_excel() %>% 
  clean_names() %>% 
  # Remove unwanted columns
  select (-c(sample, r)) %>% 
  mutate (setup = "reverse-20min",
          plate_orientation = "reverse",
          incubation_time = 20)

raw_quant_data_overnight_normal <- 
  paste0(indir, "rme_std-assay_AccuBlueNextGen_overnight-normal-plate-orientation_2024-02-23_raw-data.xlsx") %>% 
  read_excel() %>% 
  clean_names() %>% 
  # Remove unwanted columns
  select (-c(sample, r)) %>% 
  mutate (setup = "normal-overnight",
          plate_orientation = "normal",
          incubation_time = 1380)


raw_quant_data_overnight_reverse <- 
  paste0(indir, "rme_std-assay_AccuBlueNextGen_overnight-reverse-plate-orientation_2024-02-23_raw-data.xlsx") %>% 
  read_excel() %>% 
  clean_names() %>% 
  # Remove unwanted columns
  select (-c(sample, r)) %>% 
  mutate (setup = "reverse-overnight",
          plate_orientation = "reverse",
          incubation_time = 1380)

raw_quant_data <- rbind(
  raw_quant_data_5min_normal,
  raw_quant_data_5min_reverse,
  raw_quant_data_20min_normal,
  raw_quant_data_20min_reverse,
  raw_quant_data_overnight_normal,
  raw_quant_data_overnight_reverse)
  
  
# Load quant plate map
quant_plate_map_normal <- 
  paste0(indir, "rme_std-assay_AccuBlueNextGen_2024-02-22_plate-map.xlsx") %>% 
  read_excel(sheet = "plate-map") %>% 
  clean_names() %>% 
  mutate (wells = paste0(quant_row, quant_col),
          plate_orientation = "normal")

quant_plate_map_reverse <- 
  paste0(indir, "rme_std-assay_AccuBlueNextGen_2024-02-22_plate-map.xlsx") %>% 
  read_excel(sheet = "plate-map_reverse") %>% 
  clean_names() %>% 
  mutate (wells = paste0(quant_row, quant_col),
          plate_orientation = "reverse")

quant_plate_map <- 
  rbind (quant_plate_map_normal,
         quant_plate_map_reverse)

# Merge quant_data and plate_map
quant_data <- 
  full_join(quant_plate_map, raw_quant_data) %>% 
  rename (raw_rfu = value) %>% 
  na.omit() %>% 
  group_by (sample_id, plate_orientation, incubation_time) %>% 
  mutate (dna_per_ul = as.numeric(str_remove(sample_id, "^std\\-")),
          dna_per_well = dna_per_ul * sample_volume,
          strip_id = as.factor(strip_id),
          replicate = as.factor (replicate),
          plate_orientation = as.factor (plate_orientation)
          ) %>% 
  ungroup()

#### SUMMARIZE DATA ####   

# Summarize data
summary_stats <- 
  quant_data %>% 
  group_by ( strip_id, sample_id, setup, plate_orientation, incubation_time) %>% 
  summarize (n = n(),
             dna_per_ul = mean(dna_per_ul),
             sample_volume = mean(sample_volume),
             dna_per_well = mean(dna_per_well),
             mean_raw_rfu = mean(raw_rfu),
             log10_mean_raw_rfu = log10(mean_raw_rfu),
             log10_dna_per_well = log10(dna_per_well+1)
  ) %>% 
  arrange(dna_per_well) %>% 
  ungroup()


#### CURVE FITTING ####
train_data <- quant_data %>% 
  # Remove all values <= 0
  filter (!dna_per_ul <= 0)  %>%
  droplevels()

# Untransformed variables
t <- train_data %>% 
  mutate (log_dna_per_ul = log10(dna_per_ul+1),
          log_raw_rfu = log10(raw_rfu+1)) %>% 
  split (.$setup)


map2 (t, names(t),
     ~ curvefit (.x$dna_per_ul,
                 .x$raw_rfu,
                 plot.title = .y)
)

map2 (t, names(t),
      ~ curvefit (.x$log_dna_per_ul,
                  .x$raw_rfu,
                  plot.title = .y)
)

map2 (t, names(t),
      ~ curvefit (.x$log_dna_per_ul,
                  .x$log_raw_rfu,
                  plot.title = .y)
)

#### PREPARE TRAIN (Set A) AND TEST DATASET (Sets B and C) FOR REGRESSION ####

# Use only Strip A samples as the train dataset
train_data <- quant_data %>% 
  filter (plate_orientation == "normal") %>% 
  filter (strip_id == "std-setA") %>% 
  filter (incubation_time != 1380) %>% 
  droplevels() %>% 
  mutate (incubation_time = as.factor(incubation_time))

# Prepare log-transformed train data 
train_data_transformed <- 
  train_data %>%
  filter(dna_per_ul > 0) %>%
  mutate(dna_per_ul = log10(dna_per_ul),
         dna_per_well = log10(dna_per_well),
         trans_rfu = log10(raw_rfu))

# Use Strips B and C to test model
test_data <- 
  quant_data %>% 
  filter (plate_orientation == "normal") %>% 
  filter(!strip_id == "std-setA") %>% 
  filter (incubation_time != 1380) %>% 
  droplevels() %>% 
  mutate (incubation_time = as.factor(incubation_time))

# Prepare log-transformed train data 
test_data_transformed <- 
  test_data %>%
  filter(dna_per_ul > 0) %>%
  mutate(dna_per_ul = log10(dna_per_ul),
         dna_per_well = log10(dna_per_well),
         trans_rfu = log10(raw_rfu))

#### USE TRAIN DATASET TO CREATE MODEL, THEN PREDICT TEST DATASET ####

#### Fit a linear mixed effects model ####
model_linear <- 
  lme(dna_per_ul ~ raw_rfu, 
      random = list(~ 1|incubation_time),
      data = train_data)

# Displaying the model summary
summary(model_linear)

# Plot the fitted line to the model
plot_linear <- 
  ggplot (
    data = train_data,
    aes (x = raw_rfu, 
         y = dna_per_ul,
         col = incubation_time)
    ) + 
  geom_point() +
  geom_line(aes(y = predict(model_linear,
                            level = 1))
  ) +
  labs(x = "RFU (Relative Fluorescent Units)", 
       y = "DNA Concentration (pg/ul)", 
       title = "Observed vs. Predicted DNA concentration (pg/ul)",
       subtitle = "Strip A") +   
  theme_bw()
plot_linear

# Assess model residuals 
trained_model <- model_linear
plot(trained_model)
residuals(trained_model) %>% hist


#### Predict test-set ####
test_data <- 
  test_data %>% 
  mutate (predicted_dna_per_ul = predict(trained_model,
                                         newdata = .,
                                         level = 0),
          predicted_dna_per_ul_with_random_effects = predict (trained_model,
                                                              newdata = .,
                                                              level = 1)
  )


# Perform a regression between actual and predicted DNA concentration
plot_prediction_linear <- 
  ggplot(test_data, 
         aes(x = dna_per_ul, 
             y = predicted_dna_per_ul_with_random_effects,  
             #col=incubation_time,
             #group = incubation_time,
             #fill = incubation_time
         )
  ) +
  geom_smooth (method = "lm",
               se = T,
               alpha = 0.10) +
  geom_point(aes(col = incubation_time),
             alpha = 0.5,
             size = 2) + 
  theme_bw() +
  stat_poly_eq() + 
  facet_grid(~strip_id) +
  geom_abline (intercept = c(0,0), lty = 2) + 
  labs(x = "Actual Concentration (pg/ul)", y = "Predicted Concentration (pg/ul)", 
       title = "Actual vs. Predicted pg/ul using Linear Mixed Model",
       subtitle = "Model: dna_per_ul ~ raw_rfu; ~ 1 | incubation_time")
plot_prediction_linear




#### Fit a polynomial mixed effects model (2nd order) ####
model_poly2 <- 
  lme(dna_per_ul ~ trans_rfu + I(trans_rfu^2), 
      random = list(~ 1|incubation_time),
      data = train_data_transformed)

# Displaying the model summary
summary(model_poly2)

# Plot the fitted line to the model
plot_model_poly2 <- 
  ggplot (
    data = train_data_transformed,
    aes (x = trans_rfu, 
         y = dna_per_ul,
         col = incubation_time)
  ) + 
  geom_point() +
  geom_line(aes(y = predict(model_poly2,
                            level = 1))
  ) +
  labs(x = "RFU (Relative Fluorescent Units)", 
       y = "DNA Concentration (pg/ul)", 
       title = "Observed vs. Predicted DNA concentration (pg/ul)",
       subtitle = "Strip A") +   
  theme_bw()
plot_model_poly2

# Assess model residuals 
trained_model <- model_poly2
plot(trained_model)
residuals(trained_model) %>% hist


#### Predict test-set ####
test_data <- 
  test_data_transformed %>% 
  mutate (predicted_dna_per_ul = predict(trained_model,
                                         newdata = .,
                                         level = 0),
          predicted_dna_per_ul_with_random_effects = predict (trained_model,
                                                              newdata = .,
                                                              level = 1)
  )


# Perform a regression between actual and predicted DNA concentration
plot_prediction_poly2 <- 
  ggplot(test_data, 
         aes(x = 10^dna_per_ul, 
             y = 10^predicted_dna_per_ul_with_random_effects,  
             #col=incubation_time,
             #fill = incubation_time,
             #group = incubation_time
         )
  ) +
  geom_smooth (method = "lm",
               se = T,
               alpha = 0.10) +
  geom_point(aes(col = incubation_time),
             alpha = 0.5,
             size = 2) + 
  theme_bw() +
  stat_poly_eq() + 
  facet_grid(~strip_id) +
  geom_abline (intercept = c(0,0), lty = 2) + 
  labs(x = "Actual Concentration (pg/ul)", y = "Predicted Concentration (pg/ul)", 
       title = "Actual vs. Predicted pg/ul using Polynomial (2nd Order) Mixed Model",
       subtitle = "Model: log10(dna_per_ul) ~ log10(raw_rfu) + I(log10(raw_rfu^2)); ~ 1 | incubation_time")
plot_prediction_poly2


#### Fit a polynomial mixed effects model (3rd order) ####
model_poly3 <- 
  lme(dna_per_ul ~ trans_rfu + I(trans_rfu^2) + I(trans_rfu^3), 
      random = list(~ 1|incubation_time),
      data = train_data_transformed)

# Displaying the model summary
summary(model_poly3)

# Plot the fitted line to the model
plot_model_poly3 <- 
  ggplot (
    data = train_data_transformed,
    aes (x = trans_rfu, 
         y = dna_per_ul,
         col = incubation_time)
  ) + 
  geom_point() +
  geom_line(aes(y = predict(model_poly3,
                            level = 1))
  ) +
  labs(x = "RFU (Relative Fluorescent Units)", 
       y = "DNA Concentration (pg/ul)", 
       title = "Observed vs. Predicted DNA concentration (pg/ul)",
       subtitle = "Strip A") +   
  theme_bw()
plot_model_poly3

# Assess model residuals 
trained_model <- model_poly3
plot(trained_model)
residuals(trained_model) %>% hist


#### Predict test-set ####
test_data <- 
  test_data_transformed %>% 
  mutate (predicted_dna_per_ul = predict(trained_model,
                                         newdata = .,
                                         level = 0),
          predicted_dna_per_ul_with_random_effects = predict (trained_model,
                                                              newdata = .,
                                                              level = 1)
  )


# Perform a regression between actual and predicted DNA concentration
plot_prediction_poly3 <- 
  ggplot(test_data, 
         aes(x = 10^dna_per_ul, 
             y = 10^predicted_dna_per_ul_with_random_effects,  
             #col=incubation_time,
             #fill = incubation_time,
             #group = incubation_time
         )
  ) +
  geom_smooth (method = "lm",
               se = T,
               alpha = 0.10) +
  geom_point(aes(col = incubation_time),
             alpha = 0.5,
             size = 2) + 
  theme_bw() +
  stat_poly_eq() + 
  facet_grid(~strip_id) +
  geom_abline (intercept = c(0,0), lty = 2) + 
  labs(x = "Actual Concentration (pg/ul)", y = "Predicted Concentration (pg/ul)", 
       title = "Actual vs. Predicted pg/ul using Polynomial (3rd Order) Mixed Model",
       subtitle = "Model: log10(dna_per_ul) ~ log10(raw_rfu) + I(log10(raw_rfu^2)) + I(log10(raw_rfu^3)); ~ 1 | incubation_time")
plot_prediction_poly3


#### Combine Plots ####
plot_predictions <- 
  ggarrange(plot_prediction_linear,
            plot_prediction_poly2,
            plot_prediction_poly3,
            common.legend = T,
            labels = "AUTO",
            legend = "bottom"
  )
plot_predictions


ggsave (plot = plot_predictions,
        filename = paste0(outdir, "std_quant_predictions_A_", Sys.Date(), ".jpg"),
        height = 10,
        width = 13,
        dpi = 330)



#### PREPARE TRAIN (Set B) AND TEST DATASET (Sets A and C) FOR REGRESSION ####

# Use only Strip B samples as the train dataset
train_data <- quant_data %>% 
  filter (plate_orientation == "normal") %>% 
  filter (strip_id == "std-setB") %>% 
  filter (incubation_time != 1380) %>% 
  droplevels() %>% 
  mutate (incubation_time = as.factor(incubation_time))

# Prepare log-transformed train data 
train_data_transformed <- 
  train_data %>%
  filter(dna_per_ul > 0) %>%
  mutate(dna_per_ul = log10(dna_per_ul),
         dna_per_well = log10(dna_per_well),
         trans_rfu = log10(raw_rfu))

# Use Strips A and C to test model
test_data <- 
  quant_data %>% 
  filter (plate_orientation == "normal") %>% 
  filter(!strip_id == "std-setB") %>% 
  filter (incubation_time != 1380) %>% 
  droplevels() %>% 
  mutate (incubation_time = as.factor(incubation_time))

# Prepare log-transformed train data 
test_data_transformed <- 
  test_data %>%
  filter(dna_per_ul > 0) %>%
  mutate(dna_per_ul = log10(dna_per_ul),
         dna_per_well = log10(dna_per_well),
         trans_rfu = log10(raw_rfu))

#### USE TRAIN DATASET TO CREATE MODEL, THEN PREDICT TEST DATASET ####

#### Fit a linear mixed effects model ####
model_linear <- 
  lme(dna_per_ul ~ raw_rfu, 
      random = list(~ 1|incubation_time),
      data = train_data)

# Displaying the model summary
summary(model_linear)

# Plot the fitted line to the model
plot_linear <- 
  ggplot (
    data = train_data,
    aes (x = raw_rfu, 
         y = dna_per_ul,
         col = incubation_time)
  ) + 
  geom_point() +
  geom_line(aes(y = predict(model_linear,
                            level = 1))
  ) +
  labs(x = "RFU (Relative Fluorescent Units)", 
       y = "DNA Concentration (pg/ul)", 
       title = "Observed vs. Predicted DNA concentration (pg/ul)",
       subtitle = "Strip B") +   
  theme_bw()
plot_linear

# Assess model residuals 
trained_model <- model_linear
plot(trained_model)
residuals(trained_model) %>% hist


#### Predict test-set ####
test_data <- 
  test_data %>% 
  mutate (predicted_dna_per_ul = predict(trained_model,
                                         newdata = .,
                                         level = 0),
          predicted_dna_per_ul_with_random_effects = predict (trained_model,
                                                              newdata = .,
                                                              level = 1)
  )


# Perform a regression between actual and predicted DNA concentration
plot_prediction_linear <- 
  ggplot(test_data, 
         aes(x = dna_per_ul, 
             y = predicted_dna_per_ul_with_random_effects,  
             #col=incubation_time,
             #group = incubation_time,
             #fill = incubation_time
         )
  ) +
  geom_smooth (method = "lm",
               se = T,
               alpha = 0.10) +
  geom_point(aes(col = incubation_time),
             alpha = 0.5,
             size = 2) + 
  theme_bw() +
  stat_poly_eq() + 
  facet_grid(~strip_id) +
  geom_abline (intercept = c(0,0), lty = 2) + 
  labs(x = "Actual Concentration (pg/ul)", y = "Predicted Concentration (pg/ul)", 
       title = "Actual vs. Predicted pg/ul using Linear Mixed Model",
       subtitle = "Model: dna_per_ul ~ raw_rfu; ~ 1 | incubation_time")
plot_prediction_linear




#### Fit a polynomial mixed effects model (2nd order) ####
model_poly2 <- 
  lme(dna_per_ul ~ trans_rfu + I(trans_rfu^2), 
      random = list(~ 1|incubation_time),
      data = train_data_transformed)

# Displaying the model summary
summary(model_poly2)

# Plot the fitted line to the model
plot_model_poly2 <- 
  ggplot (
    data = train_data_transformed,
    aes (x = trans_rfu, 
         y = dna_per_ul,
         col = incubation_time)
  ) + 
  geom_point() +
  geom_line(aes(y = predict(model_poly2,
                            level = 1))
  ) +
  labs(x = "RFU (Relative Fluorescent Units)", 
       y = "DNA Concentration (pg/ul)", 
       title = "Observed vs. Predicted DNA concentration (pg/ul)",
       subtitle = "Strip B") +   
  theme_bw()
plot_model_poly2

# Assess model residuals 
trained_model <- model_poly2
plot(trained_model)
residuals(trained_model) %>% hist


#### Predict test-set ####
test_data <- 
  test_data_transformed %>% 
  mutate (predicted_dna_per_ul = predict(trained_model,
                                         newdata = .,
                                         level = 0),
          predicted_dna_per_ul_with_random_effects = predict (trained_model,
                                                              newdata = .,
                                                              level = 1)
  )


# Perform a regression between actual and predicted DNA concentration
plot_prediction_poly2 <- 
  ggplot(test_data, 
         aes(x = 10^dna_per_ul, 
             y = 10^predicted_dna_per_ul_with_random_effects,  
             #col=incubation_time,
             #fill = incubation_time,
             #group = incubation_time
         )
  ) +
  geom_smooth (method = "lm",
               se = T,
               alpha = 0.10) +
  geom_point(aes(col = incubation_time),
             alpha = 0.5,
             size = 2) + 
  theme_bw() +
  stat_poly_eq() + 
  facet_grid(~strip_id) +
  geom_abline (intercept = c(0,0), lty = 2) + 
  labs(x = "Actual Concentration (pg/ul)", y = "Predicted Concentration (pg/ul)", 
       title = "Actual vs. Predicted pg/ul using Polynomial (2nd Order) Mixed Model",
       subtitle = "Model: log10(dna_per_ul) ~ log10(raw_rfu) + I(log10(raw_rfu^2)); ~ 1 | incubation_time")
plot_prediction_poly2


#### Fit a polynomial mixed effects model (3rd order) ####
model_poly3 <- 
  lme(dna_per_ul ~ trans_rfu + I(trans_rfu^2) + I(trans_rfu^3), 
      random = list(~ 1|incubation_time),
      data = train_data_transformed)

# Displaying the model summary
summary(model_poly3)

# Plot the fitted line to the model
plot_model_poly3 <- 
  ggplot (
    data = train_data_transformed,
    aes (x = trans_rfu, 
         y = dna_per_ul,
         col = incubation_time)
  ) + 
  geom_point() +
  geom_line(aes(y = predict(model_poly3,
                            level = 1))
  ) +
  labs(x = "RFU (Relative Fluorescent Units)", 
       y = "DNA Concentration (pg/ul)", 
       title = "Observed vs. Predicted DNA concentration (pg/ul)",
       subtitle = "Strip B") +   
  theme_bw()
plot_model_poly3

# Assess model residuals 
trained_model <- model_poly3
plot(trained_model)
residuals(trained_model) %>% hist


#### Predict test-set ####
test_data <- 
  test_data_transformed %>% 
  mutate (predicted_dna_per_ul = predict(trained_model,
                                         newdata = .,
                                         level = 0),
          predicted_dna_per_ul_with_random_effects = predict (trained_model,
                                                              newdata = .,
                                                              level = 1)
  )


# Perform a regression between actual and predicted DNA concentration
plot_prediction_poly3 <- 
  ggplot(test_data, 
         aes(x = 10^dna_per_ul, 
             y = 10^predicted_dna_per_ul_with_random_effects,  
             #col=incubation_time,
             #fill = incubation_time,
             #group = incubation_time
         )
  ) +
  geom_smooth (method = "lm",
               se = T,
               alpha = 0.10) +
  geom_point(aes(col = incubation_time),
             alpha = 0.5,
             size = 2) + 
  theme_bw() +
  stat_poly_eq() + 
  facet_grid(~strip_id) +
  geom_abline (intercept = c(0,0), lty = 2) + 
  labs(x = "Actual Concentration (pg/ul)", y = "Predicted Concentration (pg/ul)", 
       title = "Actual vs. Predicted pg/ul using Polynomial (3rd Order) Mixed Model",
       subtitle = "Model: log10(dna_per_ul) ~ log10(raw_rfu) + I(log10(raw_rfu^2)) + I(log10(raw_rfu^3)); ~ 1 | incubation_time")
plot_prediction_poly3


#### Combine Plots ####
plot_predictions <- 
  ggarrange(plot_prediction_linear,
            plot_prediction_poly2,
            plot_prediction_poly3,
            common.legend = T,
            labels = "AUTO",
            legend = "bottom"
  )
plot_predictions


ggsave (plot = plot_predictions,
        filename = paste0(outdir, "std_quant_predictions_B_", Sys.Date(), ".jpg"),
        height = 10,
        width = 13,
        dpi = 330)


#### PREPARE TRAIN (Set C) AND TEST DATASET (Sets A and B) FOR REGRESSION ####

# Use only Strip B samples as the train dataset
train_data <- quant_data %>% 
  filter (plate_orientation == "normal") %>% 
  filter (strip_id == "std-setC") %>% 
  filter (incubation_time != 1380) %>% 
  droplevels() %>% 
  mutate (incubation_time = as.factor(incubation_time))

# Prepare log-transformed train data 
train_data_transformed <- 
  train_data %>%
  filter(dna_per_ul > 0) %>%
  mutate(dna_per_ul = log10(dna_per_ul),
         dna_per_well = log10(dna_per_well),
         trans_rfu = log10(raw_rfu))

# Use Strips A and C to test model
test_data <- 
  quant_data %>% 
  filter (plate_orientation == "normal") %>% 
  filter(!strip_id == "std-setC") %>% 
  filter (incubation_time != 1380) %>% 
  droplevels() %>% 
  mutate (incubation_time = as.factor(incubation_time))

# Prepare log-transformed train data 
test_data_transformed <- 
  test_data %>%
  filter(dna_per_ul > 0) %>%
  mutate(dna_per_ul = log10(dna_per_ul),
         dna_per_well = log10(dna_per_well),
         trans_rfu = log10(raw_rfu))

#### USE TRAIN DATASET TO CREATE MODEL, THEN PREDICT TEST DATASET ####

#### Fit a linear mixed effects model ####
model_linear <- 
  lme(dna_per_ul ~ raw_rfu, 
      random = list(~ 1|incubation_time),
      data = train_data)

# Displaying the model summary
summary(model_linear)

# Plot the fitted line to the model
plot_linear <- 
  ggplot (
    data = train_data,
    aes (x = raw_rfu, 
         y = dna_per_ul,
         col = incubation_time)
  ) + 
  geom_point() +
  geom_line(aes(y = predict(model_linear,
                            level = 1))
  ) +
  labs(x = "RFU (Relative Fluorescent Units)", 
       y = "DNA Concentration (pg/ul)", 
       title = "Observed vs. Predicted DNA concentration (pg/ul)",
       subtitle = "Strip C") +   
  theme_bw()
plot_linear

# Assess model residuals 
trained_model <- model_linear
plot(trained_model)
residuals(trained_model) %>% hist


#### Predict test-set ####
test_data <- 
  test_data %>% 
  mutate (predicted_dna_per_ul = predict(trained_model,
                                         newdata = .,
                                         level = 0),
          predicted_dna_per_ul_with_random_effects = predict (trained_model,
                                                              newdata = .,
                                                              level = 1)
  )


# Perform a regression between actual and predicted DNA concentration
plot_prediction_linear <- 
  ggplot(test_data, 
         aes(x = dna_per_ul, 
             y = predicted_dna_per_ul_with_random_effects,  
             #col=incubation_time,
             #group = incubation_time,
             #fill = incubation_time
         )
  ) +
  geom_smooth (method = "lm",
               se = T,
               alpha = 0.10) +
  geom_point(aes(col = incubation_time),
             alpha = 0.5,
             size = 2) + 
  theme_bw() +
  stat_poly_eq() + 
  facet_grid(~strip_id) +
  geom_abline (intercept = c(0,0), lty = 2) + 
  labs(x = "Actual Concentration (pg/ul)", y = "Predicted Concentration (pg/ul)", 
       title = "Actual vs. Predicted pg/ul using Linear Mixed Model",
       subtitle = "Model: dna_per_ul ~ raw_rfu; ~ 1 | incubation_time")
plot_prediction_linear




#### Fit a polynomial mixed effects model (2nd order) ####
model_poly2 <- 
  lme(dna_per_ul ~ trans_rfu + I(trans_rfu^2), 
      random = list(~ 1|incubation_time),
      data = train_data_transformed)

# Displaying the model summary
summary(model_poly2)

# Plot the fitted line to the model
plot_model_poly2 <- 
  ggplot (
    data = train_data_transformed,
    aes (x = trans_rfu, 
         y = dna_per_ul,
         col = incubation_time)
  ) + 
  geom_point() +
  geom_line(aes(y = predict(model_poly2,
                            level = 1))
  ) +
  labs(x = "RFU (Relative Fluorescent Units)", 
       y = "DNA Concentration (pg/ul)", 
       title = "Observed vs. Predicted DNA concentration (pg/ul)",
       subtitle = "Strip C") +   
  theme_bw()
plot_model_poly2

# Assess model residuals 
trained_model <- model_poly2
plot(trained_model)
residuals(trained_model) %>% hist


#### Predict test-set ####
test_data <- 
  test_data_transformed %>% 
  mutate (predicted_dna_per_ul = predict(trained_model,
                                         newdata = .,
                                         level = 0),
          predicted_dna_per_ul_with_random_effects = predict (trained_model,
                                                              newdata = .,
                                                              level = 1)
  )


# Perform a regression between actual and predicted DNA concentration
plot_prediction_poly2 <- 
  ggplot(test_data, 
         aes(x = 10^dna_per_ul, 
             y = 10^predicted_dna_per_ul_with_random_effects,  
             #col=incubation_time,
             #fill = incubation_time,
             #group = incubation_time
         )
  ) +
  geom_smooth (method = "lm",
               se = T,
               alpha = 0.10) +
  geom_point(aes(col = incubation_time),
             alpha = 0.5,
             size = 2) + 
  theme_bw() +
  stat_poly_eq() + 
  facet_grid(~strip_id) +
  geom_abline (intercept = c(0,0), lty = 2) + 
  labs(x = "Actual Concentration (pg/ul)", y = "Predicted Concentration (pg/ul)", 
       title = "Actual vs. Predicted pg/ul using Polynomial (2nd Order) Mixed Model",
       subtitle = "Model: log10(dna_per_ul) ~ log10(raw_rfu) + I(log10(raw_rfu^2)); ~ 1 | incubation_time")
plot_prediction_poly2


#### Fit a polynomial mixed effects model (3rd order) ####
model_poly3 <- 
  lme(dna_per_ul ~ trans_rfu + I(trans_rfu^2) + I(trans_rfu^3), 
      random = list(~ 1|incubation_time),
      data = train_data_transformed)

# Displaying the model summary
summary(model_poly3)

# Plot the fitted line to the model
plot_model_poly3 <- 
  ggplot (
    data = train_data_transformed,
    aes (x = trans_rfu, 
         y = dna_per_ul,
         col = incubation_time)
  ) + 
  geom_point() +
  geom_line(aes(y = predict(model_poly3,
                            level = 1))
  ) +
  labs(x = "RFU (Relative Fluorescent Units)", 
       y = "DNA Concentration (pg/ul)", 
       title = "Observed vs. Predicted DNA concentration (pg/ul)",
       subtitle = "Strip C") +   
  theme_bw()
plot_model_poly3

# Assess model residuals 
trained_model <- model_poly3
plot(trained_model)
residuals(trained_model) %>% hist


#### Predict test-set ####
test_data <- 
  test_data_transformed %>% 
  mutate (predicted_dna_per_ul = predict(trained_model,
                                         newdata = .,
                                         level = 0),
          predicted_dna_per_ul_with_random_effects = predict (trained_model,
                                                              newdata = .,
                                                              level = 1)
  )


# Perform a regression between actual and predicted DNA concentration
plot_prediction_poly3 <- 
  ggplot(test_data, 
         aes(x = 10^dna_per_ul, 
             y = 10^predicted_dna_per_ul_with_random_effects,  
             #col=incubation_time,
             #fill = incubation_time,
             #group = incubation_time
         )
  ) +
  geom_smooth (method = "lm",
               se = T,
               alpha = 0.10) +
  geom_point(aes(col = incubation_time),
             alpha = 0.5,
             size = 2) + 
  theme_bw() +
  stat_poly_eq() + 
  facet_grid(~strip_id) +
  geom_abline (intercept = c(0,0), lty = 2) + 
  labs(x = "Actual Concentration (pg/ul)", y = "Predicted Concentration (pg/ul)", 
       title = "Actual vs. Predicted pg/ul using Polynomial (3rd Order) Mixed Model",
       subtitle = "Model: log10(dna_per_ul) ~ log10(raw_rfu) + I(log10(raw_rfu^2)) + I(log10(raw_rfu^3)); ~ 1 | incubation_time")
plot_prediction_poly3


#### Combine Plots ####
plot_predictions <- 
  ggarrange(plot_prediction_linear,
            plot_prediction_poly2,
            plot_prediction_poly3,
            common.legend = T,
            labels = "AUTO",
            legend = "bottom"
  )
plot_predictions


ggsave (plot = plot_predictions,
        filename = paste0(outdir, "std_quant_predictions_C_", Sys.Date(), ".jpg"),
        height = 10,
        width = 13,
        dpi = 330)

