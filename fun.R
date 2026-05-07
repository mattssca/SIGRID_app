library(readr)
library(dplyr)
library(tidyr)

process_uromonitor <- function(file_path, assay = c("TERT", "FGFR3")) {
  assay <- match.arg(assay)
  
  # ── 1. Read ───────────────────────────────────────────────────────────────
  no_prefix   <- c("MeanQC", "CV_QC", "QC Message")
  channel_row <- read_csv(file_path, n_max = 1, col_names = FALSE,
                          show_col_types = FALSE) |> unlist() |> as.character()
  col_row     <- read_csv(file_path, skip = 1, n_max = 1, col_names = FALSE,
                          show_col_types = FALSE) |> unlist() |> as.character()
  channel <- ""
  col_names <- mapply(function(ch, base, i) {
    if (!is.na(ch) && ch != "") channel <<- ch
    if (is.na(base) || base == "") base <- paste0("X", i)
    if (channel %in% c("FAM", "VIC") && !base %in% no_prefix)
      paste0(channel, "_", base)
    else base
  }, channel_row, col_row, seq_along(col_row))
  
  df_raw <- read_csv(file_path, skip = 2, col_names = col_names,
                     col_types = cols(.default = "c"), show_col_types = FALSE) |>
    fill(Group, .direction = "down")
  
  # ── 2. Coerce numeric columns that exist in this file ─────────────────────
  coerce_df <- function(df) {
    num_cols <- intersect(
      c("Total", "FAM_Positives", "FAM_Conc.", "FAM_Lambda (cp/Rxn)",
        "VIC_Positives", "VIC_Conc.", "VIC_Lambda (cp/Rxn)", "MeanQC", "CV_QC"),
      names(df)
    )
    df |> mutate(across(all_of(num_cols), as.numeric))
  }
  
  # ── 3. Split raw into samples / controls ──────────────────────────────────
  is_patient <- function(df) df |>
    filter(!is.na(Sample), Sample != "",
           !Sample %in% c("Average", "NTC", "PC"),
           !grepl("^(0|1|UM_0)$", Sample))
  
  # ── 4. Results: TERT ──────────────────────────────────────────────────────
  if (assay == "TERT") {
    samples  <- df_raw |> is_patient() |> coerce_df()
    controls <- df_raw |> filter(Sample %in% c("NTC", "PC")) |> coerce_df()
    
    ntc_bg <- df_raw |>
      filter(Sample == "NTC", !is.na(`Run name`), `Run name` != "") |>
      mutate(FAM_Positives = as.numeric(FAM_Positives)) |>
      select(`Run name`, ntc_fam = FAM_Positives)
    
    results <- samples |>
      left_join(ntc_bg, by = "Run name") |>
      mutate(
        valid_droplets     = Total >= 20000,
        valid_ic           = VIC_Positives >= 10,
        valid              = valid_droplets & valid_ic,
        FAM_Positives_norm = pmax(0, FAM_Positives - ntc_fam),
        FA                 = (FAM_Positives / (FAM_Positives + VIC_Positives)) * 100,
        mutation_positive  = FAM_Positives_norm >= 5 & FA >= 0.5,
        result             = case_when(
          !valid            ~ "Invalid",
          mutation_positive ~ "Positive",
          TRUE              ~ "Negative"
        )
      ) |>
      select(`Run name`, Sample, Group, Well, Total,
             VIC_Positives, FAM_Positives, FAM_Positives_norm,
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
      select(`Run name`, Sample, Group, Well, Total,
             VIC_Positives, FAM_Positives, valid_droplets, valid_ic, control_pass)
    
    # ── 5. Results: FGFR3 ─────────────────────────────────────────────────────
  } else {
    is_ic_plate  <- function(x) grepl("IC",  x, ignore.case = TRUE)
    is_mut_plate <- function(x) !grepl("IC", x, ignore.case = TRUE)
    
    extract_plate_num <- function(x) {
      ifelse(grepl("[Pp]late\\d+", x),
             sub(".*[Pp]late(\\d+).*", "\\1", x),
             NA_character_)
    }
    
    ic_rows <- df_raw |>
      filter(is_ic_plate(`Run name`),
             !Sample %in% c("NTC", "PC"),
             !is.na(Sample), Sample != "",
             !grepl("^(0|1|UM_0)$", Sample)) |>
      mutate(FAM_Positives = as.numeric(FAM_Positives),
             plate_num     = extract_plate_num(`Run name`))
    
    ic_by_plate <- ic_rows |>
      filter(!is.na(plate_num)) |>
      group_by(Sample, plate_num) |>
      summarise(IC_Positives = max(FAM_Positives, na.rm = TRUE), .groups = "drop")
    
    ic_reruns <- ic_rows |>
      filter(is.na(plate_num)) |>
      group_by(Sample) |>
      summarise(IC_Positives_reruns = max(FAM_Positives, na.rm = TRUE), .groups = "drop")
    
    ntc_bg <- df_raw |>
      filter(Sample == "NTC", is_mut_plate(`Run name`),
             !is.na(`Run name`), `Run name` != "") |>
      mutate(FAM_Positives = as.numeric(FAM_Positives)) |>
      select(`Run name`, ntc_fam = FAM_Positives)
    
    results <- df_raw |>
      filter(is_mut_plate(`Run name`)) |>
      is_patient() |>
      coerce_df() |>
      mutate(
        mut_target = case_when(
          grepl("248", `Run name`) ~ "FGFR3_248",
          grepl("249", `Run name`) ~ "FGFR3_249",
          TRUE                     ~ NA_character_
        ),
        plate_num = extract_plate_num(`Run name`)
      ) |>
      left_join(ntc_bg,      by = "Run name") |>
      left_join(ic_by_plate, by = c("Sample", "plate_num")) |>
      left_join(ic_reruns,   by = "Sample") |>
      mutate(
        IC_Positives       = coalesce(IC_Positives, IC_Positives_reruns),
        ntc_fam            = replace_na(ntc_fam, 0),
        FAM_Positives_norm = pmax(0, FAM_Positives - ntc_fam),
        valid_droplets     = Total >= 20000,
        valid_ic           = !is.na(IC_Positives) & IC_Positives >= 10,
        valid              = valid_droplets & valid_ic,
        FA                 = (FAM_Positives / (FAM_Positives + IC_Positives)) * 100,
        mutation_positive  = FAM_Positives_norm >= 5 & FA >= 0.5,
        result             = case_when(
          !valid            ~ "Invalid",
          mutation_positive ~ "Positive",
          TRUE              ~ "Negative"
        )
      ) |>
      select(`Run name`, Sample, Group, Well, Total, mut_target,
             IC_Positives, FAM_Positives, FAM_Positives_norm,
             FA, valid_droplets, valid_ic, valid, mutation_positive, result)
    
    ctrl_results <- df_raw |>
      filter(is_mut_plate(`Run name`), Sample %in% c("NTC", "PC")) |>
      coerce_df() |>
      mutate(
        valid_droplets = Total >= 20000,
        control_pass   = case_when(
          Sample == "NTC" ~ valid_droplets & FAM_Positives == 0,
          Sample == "PC"  ~ valid_droplets & FAM_Positives >= 5
        )
      ) |>
      select(`Run name`, Sample, Group, Well, Total,
             FAM_Positives, valid_droplets, control_pass)
  }
  
  list(results = results, controls = ctrl_results)
}
