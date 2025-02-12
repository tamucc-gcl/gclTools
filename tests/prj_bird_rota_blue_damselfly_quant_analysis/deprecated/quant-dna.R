#### quant-dna ####
## Kevin Labrador
## 2024-02-07

#### INTRODUCTION ####
# This script uses R to analyze the flourescence data exported from the SpectraMax software.
# There are two files needed for this pipeline:

## 1. raw plate reader data: this is the file exported from the SpectraMax software. It contains the well and fluourescence information.

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
pacman::p_load(tidyverse,
               janitor,
               ggpubr,
               readxl,
               scales,
               magrittr,
               broom,
               cowplot,
               rstatix,
               ggpmisc
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
raw_quant_data <- 
  paste0(indir, "rme_preservative_dna-extracts_accublue-nextgen_raw-data_2024-11-21.xlsx") %>% 
  read_excel() %>% 
  clean_names() %>% 
  
  # Remove unwanted columns
  select (-c(sample, r))
  
# Load quant plate map
quant_plate_map <- 
  paste0(indir, "rme_preservative_dna-extracts_accublue-nextgen_plate-map_2024-11-21.xlsx") %>% 
  read_excel(sheet = "plate-map") %>% 
  clean_names() %>% 
  mutate (wells = paste0(quant_row, quant_col))

# Load dna plate map 
dna_plate_map <- 
  paste0("../data/rme_preservative_dna-extract-plate-map.xlsx") %>%  
  read_excel(sheet = "plate-map_tidy") %>% 
  clean_names() %>% 
  mutate (sample_id = paste0(plate_row, plate_col)) %>% 
  select (c(plate_id, dna_extract_tube_id, preservative, sample_id))

# Merge quant_data and plate_map
quant_data <- 
  full_join(quant_plate_map, raw_quant_data, by = "wells") %>% 
  rename (rfu = value) %>% 
  na.omit()

#### REGRESSION MODEL FROM STANDARD ####

# Pull out the standards from the dataset
standard_df <- 
  quant_data %>% 
  filter (grepl("std", sample_id),
          grepl ("setA", plate_id)
          ) %>%
  mutate (dna_per_ul = as.numeric(str_remove(sample_id, "std-")),
          dna_per_well = dna_per_ul*sample_volume)


#### Following the kit instructions ####
## Average the triplicate values for each standard sample and subtract the average zero DNA value from each data point. 
standard_mean <- 
  standard_df %>% 
  group_by (sample_id) %>% 
  summarize (n = n(),
             dna_per_ul = mean(dna_per_ul),
             sample_volume = mean(sample_volume),
             dna_per_well = mean(dna_per_well),
             mean_rfu = mean(rfu),
             )%>% 
  arrange(dna_per_well)

# Extract the background rfu 
background_rfu <- 
  standard_mean %>% 
  filter (sample_id == "std-0") %>% 
  pull (mean_rfu)

# Prepare training dataset by subtracting the background rfu to the rest of the standards
train_data <- 
  standard_mean %>% 
  mutate (zeroed_rfu = mean_rfu - background_rfu) 

# Prepare test dataset by subtracting the background rfu to the samples
test_data <- 
  quant_data %>% 
  filter (!grepl("std", sample_id)) %>% 
  mutate (zeroed_rfu = rfu - background_rfu)

# Prepare a linear regression model
train_data_linear <- train_data 

model_linear <- lm (dna_per_well ~  0 + zeroed_rfu,
                    data = train_data_linear)

# Add fitted model to train data
train_data_linear$fitted_values <- model_linear$fitted.values

# Predict test data using model
test_data_linear <- test_data %>%
  mutate (dna_per_well = predict (model_linear, 
                                    newdata = .), 
          dna_per_ul = dna_per_well/sample_volume
  ) %>% 
  group_by (sample_id) %>% 
  summarize (n = n(),
             mean_zeroed_rfu = mean(zeroed_rfu),
             sd_zeroed_rfu = sd(zeroed_rfu),
             mean_dna_per_well = mean(dna_per_well),
             sd_dna_per_well = sd(dna_per_well),
             mean_dna_per_ul = mean(dna_per_ul),
             sd_dna_per_ul = sd(dna_per_ul)) %>% 
  full_join(dna_plate_map)



# Plot the standard curve with raw data
ggplot (data = train_data,
        aes (x = zeroed_rfu,
             y = dna_per_well)) +
  geom_text (aes(label = sample_id),
             nudge_x = 2,
             nudge_y = -2,
             size = 3) + 
  # geom_errorbarh (aes(xmin = zeroed_rfu - sd_zeroed_rfu,
  #                     xmax = zeroed_rfu + sd_zeroed_rfu,
  #                     height = 20)
  #                 )+
  geom_smooth (method = "lm",
               formula = y ~ 0 + x,
               alpha = 0.25) +
  stat_poly_eq(use_label(c("eq", "R2", "n"))) +
  geom_point(col = "red", size = 3) +
  geom_point(data = standard_df,
             aes (x = rfu - background_rfu, 
                  y = dna_per_well),
             col = "gray50",
             alpha = 1,
             size = 2) +
  theme_bw() +
  labs (title = "Linear Regression of Standards",
        x = "Zeroed RFU",
        y = "DNA Concentration (pg/well)") 


# Plot the standard curve with samples
ggplot (data = train_data,
        aes (x = zeroed_rfu,
             y = dna_per_well)) +
  geom_smooth (method = "lm",
               alpha = 0.25) +
  stat_poly_eq(use_label(c("eq", "R2", "n"))) +
  geom_point(col = "red", size = 3) +
  geom_point(data = test_data_linear,
             aes(x = mean_zeroed_rfu, 
                 y = mean_dna_per_well),
             col = "gray50") +
  geom_errorbar(data = test_data_linear,
                aes(x = mean_zeroed_rfu, 
                    y = mean_dna_per_well,
                    ymin = mean_dna_per_well - sd_dna_per_well,
                    ymax = mean_dna_per_well + sd_dna_per_well),
                col = "gray50") +  
  geom_errorbarh(data = test_data_linear,
                 aes(x = mean_zeroed_rfu, 
                     y = mean_dna_per_well,
                     xmin = mean_zeroed_rfu - sd_zeroed_rfu,
                     xmax = mean_zeroed_rfu + sd_zeroed_rfu),
                 col = "gray50") +
  geom_text(data = test_data_linear,
            aes(x = mean_zeroed_rfu, 
                y = mean_dna_per_well,
                label = dna_extract_tube_id),
            nudge_x = 5,
            col = "black") +  
  theme_bw() +
  labs (title = "Linear Regression of Standards",
        x = "Zeroed RFU",
        y = "DNA Concentration (pg/well)")



#### Getting the log of observations ####
standard_mean <- 
  standard_df %>% 
  group_by (sample_id) %>% 
  summarize (n = n(),
             dna_per_ul = mean(dna_per_ul),
             sample_volume = mean(sample_volume),
             dna_per_well = mean(dna_per_well),
             mean_rfu = mean(rfu),
  )%>% 
  arrange(dna_per_well)

# Extract the background rfu 
background_rfu <- 
  standard_mean %>% 
  filter (sample_id == "std-0") %>% 
  pull (mean_rfu)

# Prepare training dataset by subtracting the background rfu to the rest of the standards
train_data <- 
  standard_df %>% 
  mutate (zeroed_rfu = rfu - background_rfu) %>% 
  group_by (sample_id) %>% 
  summarize (n = n(),
             dna_per_ul = mean(dna_per_ul),
             sample_volume = mean(sample_volume),
             dna_per_well = mean(dna_per_well),
             mean_rfu = mean(rfu),
             sd_rfu = sd(rfu),
             mean_zeroed_rfu = mean(zeroed_rfu),
             sd_zeroed_rfu = sd(zeroed_rfu)
  )%>% 
  rename (zeroed_rfu = mean_zeroed_rfu) %>% 
  arrange(dna_per_well)

# Prepare test dataset by subtracting the background rfu to the samples
test_data <- 
  quant_data %>% 
  filter (plate_id != "standard") %>% 
  mutate (zeroed_rfu = rfu - background_rfu)


#### REGRESSION MODEL FROM STANDARDS ####

# We use the DNA concentration of the standards to formulate a regression from which to calculate the DNA concentration of the samples.
# We will fit both a linear model (y = mx+b; kit-based method) and a power model (y = ax^b; Chris' method).

#### LINEAR REGRESSION MODEL ####
train_data_linear <- train_data

model_linear <- lm (dna_per_well ~  zeroed_rfu,
                    data = train_data_linear)

# Add fitted model to train data
train_data_linear %<>% 
  mutate (fitted_values = model_linear$fitted.values) #%>% 
  #filter (!dna_per_well <= 0)


# Predict test data using model
test_data_linear <- test_data %>%
  mutate (dna_per_well_linear = predict (model_linear, 
                             newdata = .), 
            dna_per__ul_linear = dna_per_well_linear/sample_volume
  )

# Plot the dataset
plot_linear <- 
ggplot () + 
  geom_point (data = train_data_linear, 
              aes (x = log10(zeroed_rfu+11),
                   y = log10(dna_per_well))) +
  geom_line (data = train_data_linear, 
              aes (x = log10(zeroed_rfu),
                   y = log10(fitted_values)),
             col = "red") + 
  geom_point (data = test_data_linear,
              aes (x = log10(zeroed_rfu),
                   y = log10(dna_per_well_linear)),
              col = "gray"
  ) +

  # Annotate
  annotate ("text",
            x = -Inf,
            y = Inf,
            hjust = -0.5,
            vjust = 2,
            label = sprintf("y = %.4fx\nR² = %.3f", 
                            coef(model_linear)[1],
                            summary(model_linear)$r.squared),
            size = 3.5
  ) +
  
  theme_bw() +
  labs (x = "Zeroed RFU",
        y = "DNA concentration (pg/well)"
  )
plot_linear


ggplot (data = train_data_linear, 
        aes (x = log10(zeroed_rfu+ 11),
             y = log10(dna_per_well))) + 
  geom_point () + 
  geom_smooth(data = train_data_linear %>% 
                filter (!zeroed_rfu <= 0), 
              method = "lm")
##  
ggplot (data = train_data_linear, 
        aes (x = log10(zeroed_rfu+ 11),
             y = log10(dna_per_well+0.01))) + 
  geom_point () +
  geom_text (aes(label = sample_id),
             nudge_x = 0.20,
             nudge_y = 0.20) +
  geom_smooth(data = train_data_linear %>% 
                filter (!zeroed_rfu == 0), 
              aes (x = log10(zeroed_rfu+11),
                   y = log10(dna_per_well+0.01)),
              method = "loess",
              inherit.aes =  F)


#### POWER REGRESSION MODEL ####

train_data_power <- train_data %>% 
  
  # Remove all values <= 0
  filter (!zeroed_rfu <= 0) #%>% 
  
  # Do only if the 125 ng standard did not turn out great.
  #filter (!zeroed_rfu >= max(.$zeroed_rfu))

model_power <- lm (log(dna_per_well) ~ log(zeroed_rfu),
                    data = train_data_power)

# Add fitted model to train data
train_data_power %<>% 
  mutate (fitted_values = exp(model_power$fitted.values))


# Predict test data using model
test_data_power<- test_data %>%
  mutate (dna_per_well_power = exp(predict (model_power, 
                                    newdata = .)), 
          dna_per_ug_power = dna_per_well_power/sample_volume
  )

# Plot the dataset
plot_power <-  

ggplot () + 
  geom_point (data = train_data_power, 
              aes (x = log(zeroed_rfu),
                   y = log(ng_well))) +
  geom_line (data = train_data_power, 
             aes (x = log(zeroed_rfu),
                  y = log(fitted_values)),
             col = "red") + 
  geom_point (data = test_data_power,
              aes (x = log(zeroed_rfu),
                   y = log(ng_well_power)),
              col = "gray" 
  ) +
  geom_text (data = test_data_power,
             aes(label = sample_id,
             x = log(zeroed_rfu),
             y = log(ng_well_power))) + 
  
  # Annotate
  annotate ("text",
            x = -Inf,
            y = Inf,
            hjust = -0.5,
            vjust = 2,
            label = sprintf("y = %.4fx^%.4f\n   R² = %.3f", 
                            exp(coef(model_power)[1]),
                            coef(model_power)[2],
                            summary(model_power)$r.squared),
            parse = F,
            size = 3.5
  ) +
  
  
  theme_bw() +
  labs (x = expression(ln~"[Zeroed RFU]"),
        y = expression(ln~"[DNA concentration (ng/well)]")
  )

plot_power

#### COMPARE CONCENTRATION BETWEEN REGRESSION MODELS ####

quant_data <- full_join (test_data_linear, test_data_power) 
quant_data_pivot <- quant_data %>% 
  pivot_longer(cols = c(ng_ul_linear, ng_ul_power),
               names_to = "model",
               names_prefix = "ng_ul_",
               values_to = "dna_concentration") %>% 
  filter (!dna_concentration < 0)

plot_scatter <- 
ggplot (quant_data,
        aes (x = ng_ul_linear, y = ng_ul_power)) + 
  geom_point() +
  geom_abline(intercept = c(0,0), lty = 2, col = "black") +
  theme_bw() +
  labs (x = "DNA concentration (ng/uL)\nLinear Model",
        y = "DNA concentration (ng/uL)\nPower Model")
plot_scatter

plot_box <- quant_data_pivot %>% 
  ggplot (aes ( x = model, y = dna_concentration)) + 
  geom_boxplot() + 
  geom_jitter(col = "gray")+
  labs (x = "Regression Model", y = "DNA concentration (ng/uL)") +
  theme_bw() +
  stat_compare_means(method = "wilcox.test", paired = T)
plot_box

# Plot differences of DNA concentration calculated with linear and power model



#### VISUALIZE DATASET ####
plot <- plot_grid (
  plot_linear + labs(subtitle = "Linear regression of standards"),
  plot_power + labs (subtitle = "Power regression of standards"),
  plot_scatter + labs (subtitle = "DNA concentration between linear and\npower regression models"), 
  plot_box + labs (subtitle = "Comparison of calculated DNA concentration\nusing paired Wilcoxon Test"),
  labels = "AUTO",
  align = "hv",
  nrow = 2)

plot


#### GENERATE OUTPUT ####

# Save plot
ggsave (plot = plot,
        file = paste0(outdir, "quant-plot_", Sys.Date(), ".png"),
        height = 7.5,
        width = 7.5,
        dpi = 330,
        units = "in")
                      
        
# Save data file
write.csv (quant_data,
           file = paste0(outdir, "quant-data_", Sys.Date(),".csv")
)




