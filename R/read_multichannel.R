library(readr)
library(dplyr)
library(tidyr)

# Reads a single channel_based_multichannel CSV and returns a wide data frame
# with one row per sample and columns: Run name, Sample, Well, Total,
# FAM_Positives, VIC_Positives.
read_multichannel <- function(file_path) {
  raw <- read_csv(file_path, show_col_types = FALSE) |>
    rename(`Run name` = Run)

  totals <- raw |>
    filter(Channels == "FAM") |>
    select(`Run name`, Sample, Well, Total)

  counts <- raw |>
    filter(Channels %in% c("FAM", "VIC")) |>
    select(`Run name`, Sample, Well, Channels, Count) |>
    pivot_wider(
      names_from  = Channels,
      values_from = Count,
      names_glue  = "{Channels}_Positives"
    )

  left_join(totals, counts, by = c("Run name", "Sample", "Well")) |>
    mutate(
      Total         = as.numeric(Total),
      FAM_Positives = as.numeric(FAM_Positives),
      VIC_Positives = as.numeric(VIC_Positives)
    )
}
