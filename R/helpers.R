library(dplyr)

# Returns only patient rows — excludes controls, blanks, and placeholder IDs.
is_patient <- function(df) {
  df |>
    filter(
      !is.na(Sample), Sample != "",
      !Sample %in% c("Average", "NTC", "PC"),
      !grepl("^(0|1|UM_0)$", Sample)
    )
}


