source("R/read_multichannel.R")
source("R/helpers.R")

library(dplyr)
library(tidyr)

# Processes a single TERT-146 channel_based_multichannel CSV.
# Returns a list with two data frames:
#   $results  — one row per patient sample with validity flags and mutation call
#   $controls — NTC and PC rows with pass/fail flags
process_tert <- function(file_path) {
  df_raw <- read_multichannel(file_path)

  samples  <- df_raw |> is_patient()
  controls <- df_raw |> filter(Sample %in% c("NTC", "PC"))

  results <- samples |>
    mutate(
      valid_droplets    = Total >= 20000,
      valid_ic          = VIC_Positives >= 10,
      valid             = valid_droplets & valid_ic,
      FA                = (FAM_Positives / (FAM_Positives + VIC_Positives)) * 100,
      mutation_positive = FAM_Positives >= 5 & FA >= 0.5,
      result            = case_when(
        !valid            ~ "Invalid",
        mutation_positive ~ "Positive",
        TRUE              ~ "Negative"
      )
    ) |>
    select(`Run name`, Sample, Well, Total,
           VIC_Positives, FAM_Positives,
           FA, valid_droplets, valid_ic, valid, mutation_positive, result)

  ctrl_results <- controls |>
    mutate(
      valid_droplets = Total >= 20000,
      valid_ic       = VIC_Positives >= 10,
      control_pass   = case_when(
        Sample == "NTC" ~ valid_droplets & FAM_Positives == 0,
        Sample == "PC"  ~ valid_droplets & valid_ic
      )
    ) |>
    select(`Run name`, Sample, Well, Total,
           VIC_Positives, FAM_Positives, valid_droplets, valid_ic, control_pass)

  list(results = results, controls = ctrl_results)
}
