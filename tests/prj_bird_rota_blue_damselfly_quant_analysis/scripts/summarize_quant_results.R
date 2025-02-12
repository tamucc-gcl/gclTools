#### calculate_mean_quant_plate04.R ####
## Kevin Labrador
## 2025-02-06

# Set-up the working directory in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Clear global environment.
rm(list = ls())

#### PROGRAM DEFINED VARIABLES

#### Load Libraries ####
pacman::p_load(
  janitor,
  readxl,
  tidyverse
)

#### USER DEFINED ####
# Set infiles
quant_files_to_summarize <-
  list.files(
    "../data_processed/",
    pattern = "_quant_report_replicates.csv",
    full.names = T
  ) %>%
  map_df (~.x %>% read_csv)


# Assign output file paths
path_quant_report_summarized <-
  paste0(
    "../data/processed/",
    "rbd_edna_quant_report_merged_summarized.csv"
  )


path_quant_plot <-
  paste0(
    "../results/",
    "rbd_quant_sample_assessment.png"
  )

#### Summarize quant results ####

quant_files_summarized <-
  quant_files_to_summarize %>%
  group_by(
    plate_id,
    plate_column,
    plate_row,
    sample_type,
    dna_extract_tube_id
  ) %>%
  summarize(
    # Number of replicates
    n_reps = n(),

    # Means (na.rm = TRUE ensures that NA values won't break the calculation)
    sample_volume_mean = mean(sample_volume, na.rm = TRUE),
    rfu_mean           = mean(rfu, na.rm = TRUE),
    dna_per_well_mean  = mean(dna_per_well, na.rm = TRUE),
    ng_per_ul_mean     = mean(ng_per_ul, na.rm = TRUE),

    # Standard deviation of ng_per_ul
    ng_per_ul_sd = sd(ng_per_ul, na.rm = TRUE),

    # CV of ng_per_ul (SD / Mean * 100)
    ng_per_ul_cv = (sd(ng_per_ul, na.rm = TRUE) / mean(ng_per_ul, na.rm = TRUE)) * 100,

    # 95% CI for ng_per_ul (using a t-based approach, suitable for small sample sizes)
    ng_per_ul_ci_95 = qt(0.975, df = n() - 1) * (sd(ng_per_ul, na.rm = TRUE) / sqrt(n())),

    # % of replicates that are "within_lod"
    pct_reps_within_lod = (sum(quant_status == "within_lod", na.rm = TRUE) / n()) * 100,

    pct_reps_outside_lod = 100 - pct_reps_within_lod
  ) %>%
  ungroup() %>%
  mutate (
    sample_well_id =
      as.factor(
        sprintf("%s%02d",
                plate_row,
                plate_column
        )
      )
  )


#### Visualize quant results ####

(quant_plot <-
   quant_files_summarized  %>%
   ggplot (
     aes (
       x = fct_reorder (
         sample_well_id,
         desc(sample_well_id)
       ),
       y = ng_per_ul_mean,
       shape = sample_type,
       col = pct_reps_outside_lod
     )
   )  +
   geom_point(size = 1.5) +
   scale_color_gradientn(
     colors = c("gray50", "pink", "red3"),
     values = c(0, 0.5, 1),
     limits = c(0, 100),
     breaks = c(0, 50, 100)
   )+
   geom_errorbar(
     aes(ymin = ng_per_ul_mean - ng_per_ul_ci_95,
         ymax = ng_per_ul_mean + ng_per_ul_ci_95
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
   coord_flip() +
   facet_wrap (~plate_id,
               ncol = 5,
               scales =  "free_y") +
   theme_bw() +
   labs (
     y = "DNA Concentration (ng/uL) \u00b1 95% CI",
     x = "Sample Well ID",
     color = "Percent replicates outside\nlimit of detection",
     shape = "Sample Type"
   ) +
   theme (axis.text.x = element_text (size = 7))
)


#### OUTPUT RESULTS ####
ggsave(
  quant_plot,
  file = path_quant_plot,
  width = 12,
  height = 10,
  units = "in",
  dpi = 330)


# Export data frame
quant_files_summarized %>%
  select(
    -sample_well_id
    -contains("_cv"),
    -contains("_rfu"),
    -contains("per_well"),
    -contains("_sd")
  ) %>%
  arrange (plate_id, plate_row, plate_column) %>%
  write_csv(
    file = path_quant_report_summarized
  )
