#### troubleshooting gclTools repo changes ####

# from gclTools repo, run tamuccGCL.Rproj
# or, if already in gclTools project, just reload project
# from within r studio (upper right, gclTools pull down, select gclTools)
require(devtools)
load_all()

#### INITIALIZE  ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### USER DEFINED VARIABLES ####
plate_map_path = "./rbd_edna-plate01_AccuClear_2024-05-14_plate-map.csv"
quant_data_path = "./rbd_edna-plate01_AccuClear_2024-05-14_raw-data.csv"

#### LOAD PACKAGES  ####
# if troubleshooting gclTools repo, don't run this, otherwise, run this
# require (gclTools)

#### AUTO GEN VARIABLES ####

formatted_plate_map_path <-
  str_c(
    basename(plate_map_path),
    "_formatted",
    ".csv"
  )


# formatted_quant_data_path <-
#   str_c(
#     basename(quant_data_path),
#     "_formatted",
#     ".csv"
#   )

#### WRANGLE DATA INTO GCLTOOLS FORMAT ####
# Import data files

# raw_data <-
#   read_csv(quant_data_path) %>%
#   select(where(~ !all(is.na(.)))) %>%
#
#   mutate(Sample = if_else(row_number() == 1, 1, NA_integer_)) %>%
#   mutate(
#     R = "R",
#     .after = "Value"
#   )

plate_map <-
  read_csv(plate_map_path) %>%
  mutate(sample_id =
           str_remove(
             sample_id,
             "std\\-"
           )) %>%
  select(where(~ !all(is.na(.))))

# write_csv(
#   raw_data,
#   formatted_quant_data_path
# )

write_csv(
  plate_map,
  formatted_plate_map_path
)

#### READ IN DATA ####

quant_data <-
  load_quant_files(
    quant_data_path,
    formatted_plate_map_path
  )

#### RUN & Visualize MODELS ####

# kit_reference
linear_zeroint_zerostd_zeroedrfu_model <-
  train_quant_standards(
    quant_data,
    model = "kit_reference",
    remove_standard = FALSE
  )

linear_zeroint_zerostd_zeroedrfu_model

linear_zeroint_zerostd_zeroedrfu_model$standard %>%
  filter(dna_per_ul != 0) %>%
  ggplot() +
  aes(
    x = log10(zeroed_rfu),
    y = log10(dna_per_well)
  ) +
  geom_point() +
  geom_line(
    aes(y=log10(fit))
  ) +
  theme_classic()

# linear
linear_freeint_zerostd_zeroedrfu_model <-
  train_quant_standards(
    quant_data,
    model = "linear",
    remove_standard = FALSE
  )
linear_freeint_zerostd_zeroedrfu_model

linear_freeint_zerostd_zeroedrfu_model$standard %>%
  filter(dna_per_ul != 0) %>%
  ggplot() +
  aes(
    x = log10(zeroed_rfu),
    y = log10(dna_per_well)
  ) +
  geom_point() +
  geom_line(
    aes(y=log10(fit))
  ) +
  theme_classic()

# linear_nozerostd_zeroedrfu
linear_freeint_nozerostd_zeroedrfu_model <-
  train_quant_standards(
    quant_data,
    model = "linear_nozerostd_zeroedrfu",
    remove_standard = FALSE  # FALSE or a vector of the std sample_id's to exclude
  )
linear_freeint_nozerostd_zeroedrfu_model

linear_freeint_nozerostd_zeroedrfu_model$standard %>%
  filter(dna_per_ul != 0) %>%
  ggplot() +
  aes(
    x = log10(zeroed_rfu),
    y = log10(dna_per_well)
  ) +
  geom_point() +
  geom_line(
    aes(y=log10(fit))
  ) +
  theme_classic()

# power
power_freeint_nozerostd_zeroedrfu_model <-
  train_quant_standards(
    quant_data,
    model = "power",
    remove_standard = FALSE
  )
power_freeint_nozerostd_zeroedrfu_model

power_freeint_nozerostd_zeroedrfu_model$standard %>%
  filter(dna_per_ul != 0) %>%
  ggplot() +
  aes(
    x = log(zeroed_rfu),
    y = log(dna_per_well)
  ) +
  geom_point() +
  geom_line(
    aes(y=log(fit))
  ) +
  theme_classic()

power_freeint_nozerostd_zeroedrfu_model$standard %>%
  filter(dna_per_ul != 0) %>%
  ggplot() +
  aes(
    x = zeroed_rfu,
    y = dna_per_well
  ) +
  geom_point() +
  geom_line(
    aes(y=fit)
  ) +
  theme_classic()

# power_rawrfu
power_freeint_nozerostd_rawrfu_model <-
  train_quant_standards(
    quant_data,
    model = "power_rawrfu",
    remove_standard = FALSE
  )
power_freeint_nozerostd_rawrfu_model

power_freeint_nozerostd_rawrfu_model$standard %>%
  filter(dna_per_ul != 0) %>%
  ggplot() +
  aes(
    x = log(zeroed_rfu),
    y = log(dna_per_well)
  ) +
  geom_point() +
  geom_line(
    aes(y=log(fit))
  ) +
  theme_classic()

power_freeint_nozerostd_rawrfu_model$standard %>%
  filter(dna_per_ul != 0) %>%
  ggplot() +
  aes(
    x = zeroed_rfu,
    y = dna_per_well
  ) +
  geom_point() +
  geom_line(
    aes(y=fit)
  ) +
  theme_classic()

#### Create and Print Report Output ####
linear_zeroint_zerostd_zeroedrfu_report <-
  quant_dna(
    quant_data = quant_data,
    trained_model = linear_zeroint_zerostd_zeroedrfu_model
  )
linear_zeroint_zerostd_zeroedrfu_report

linear_freeint_zerostd_zeroedrfu_report <-
  quant_dna(
    quant_data = quant_data,
    trained_model = linear_freeint_zerostd_zeroedrfu_model
  )
linear_freeint_zerostd_zeroedrfu_report

power_freeint_nozerostd_zeroedrfu_report <-
  quant_dna(
    quant_data = quant_data,
    trained_model = power_freeint_nozerostd_zeroedrfu_model
  )
power_freeint_nozerostd_zeroedrfu_report

power_freeint_nozerostd_rawrfu_report <-
  quant_dna(
    quant_data = quant_data,
    trained_model = power_freeint_nozerostd_rawrfu_model
  )
power_freeint_nozerostd_rawrfu_report

## Not run:
export_quant_report (quant_report)

## End(Not run)
