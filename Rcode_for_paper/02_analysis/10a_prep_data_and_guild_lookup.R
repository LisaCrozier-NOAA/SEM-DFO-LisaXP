# ==============================================================================
# Script: 10a_prep_data_and_guild_lookup.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Single-execution data preparation, 19yr/24yr filtering, Min-Max [0, 1]
#          scaling, CSV export for Excel, and dynamic Guild Lookup Table creation.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------------------------
# STEP 0: CONFIGURATION SWITCH & PATHS
# ------------------------------------------------------------------------------
#years_included <- "19yr" # Switch: "19yr" or "24yr"
years_included <- "24yr" # Switch: "19yr" or "24yr"

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# STEP 1: APPLY 19-YEAR VS 24-YEAR COLUMN FILTERING
# ------------------------------------------------------------------------------
if (years_included == "24yr") {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
} else {
  # 19yr: Keep columns with at least 19 years of observations
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) >= 19]
}

filtered_df <- raw_df %>% 
  select(all_of(valid_cols)) %>% 
  filter(complete.cases(.))

cat(sprintf("[%s Selection] Retained %d columns across %d complete rows.\n", 
            years_included, ncol(filtered_df), nrow(filtered_df)))

# ------------------------------------------------------------------------------
# STEP 2: MIN-MAX SCALING [0, 1] & EXCEL-FRIENDLY CSV EXPORT
# ------------------------------------------------------------------------------
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

minmax_complete_data <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), min_max_scale))

# Export directly to CSV for easy inspection in Excel
write_csv(minmax_complete_data, file.path(data_out_dir, paste0("minmax_complete_data_", years_included, ".csv")))

# ------------------------------------------------------------------------------
# STEP 3: DEFINE FUNCTIONAL PREY POOLS & SPATIAL OVERLAP BY GUILD
# ------------------------------------------------------------------------------
all_cols <- names(minmax_complete_data)

get_prey_pool <- function(guild_num, stage_suffix) {
  sfx <- if (stage_suffix != "") paste0("_", stage_suffix) else ""
  
  m_squid     <- paste0("X04_marketsquid_GAM", sfx)
  anchovy     <- paste0("X05_anchovy_GAM", sfx)
  sardine     <- paste0("X05_DFA_abundSardine", sfx)
  herring_ncc <- paste0("X05_herring_GAM", sfx)
  eulachon    <- paste0("X05_eulachon_during_chinook", sfx)
  hake        <- paste0("X09_DFA_HakeAge5Plus", sfx)
  
  pollock     <- paste0("X13_pollock_age1plus", sfx)
  herring_ak  <- paste0("X13_sitkaHerring_EGoA", sfx)
  capelin     <- paste0("X13_capelin_WGoA", sfx)
  krill_ak    <- paste0("X12_egoa_krill", sfx)
  pink_ak     <- paste0("X14_pinkSalmonNorthAmerica", sfx)
  
  prey_list <- switch(
    as.character(guild_num),
    "1" = c(m_squid, anchovy, sardine, herring_ncc, eulachon),
    "2" = c(hake, m_squid, anchovy, sardine, herring_ncc, eulachon),
    "3" = c(pollock, herring_ak, capelin, krill_ak),
    "4" = c(pollock, herring_ak, capelin),
    "5" = c(pink_ak, pollock, herring_ak, capelin),
    c()
  )
  
  intersect(prey_list, all_cols)
}

get_spatial_overlap_col <- function(pred_col, region, stage_suffix) {
  sfx <- if (stage_suffix != "") paste0("_", stage_suffix) else ""
  
  if (grepl("X09_DFA_HakeAge5Plus", pred_col)) {
    col <- paste0("X12_DFA_biomassEuphShelfSum", sfx)
  } else if (region == "NCC") {
    col <- paste0("X01_habCompInd", sfx)
  } else {
    col <- paste0("X21_sst_egoa_junjulaug", sfx)
  }
  
  if (col %in% all_cols) return(col) else return(NA_character_)
}

# ------------------------------------------------------------------------------
# STEP 4: DYNAMICALLY BUILD GUILD LOOKUP TABLE
# ------------------------------------------------------------------------------
predator_cols <- all_cols[grepl("^(X08_|X09_|X10_|X11_|X15_)", all_cols)]

blacklist <- c("canaryRockfish", "chilipepper", "Harbor_seal_CR",
               "X09_DFA_ChinAbundSnakeFall_adultyr", "X09_DFA_HakeAge5Plus_adultyr",
               "sharkCatchGoA", "X15_spinyDogfishGoA_predAK_adultyr")
predator_cols <- predator_cols[!sapply(predator_cols, function(x) any(sapply(blacklist, grepl, x)))]

lookup_rows <- list()

for (pred in predator_cols) {
  region <- if (grepl("^(X08_|X09_|X10_)", pred) || pred == "X11_DFA_Harbour_p_WS") "NCC" else "AK"
  
  stage_suffix <- if (grepl("smoltyr$", pred)) {
    "smoltyr"
  } else if (grepl("adultyr$", pred)) {
    "adultyr"
  } else {
    "smoltyr"
  }
  
  guild_num <- case_when(
    region == "NCC" & grepl("^(X08_|X11_DFA_Harbour_p_WS)", pred) ~ 1,
    region == "NCC" ~ 2,
    region == "AK" & grepl("(ssl|sablefish|spinyDogfish)", pred, ignore.case = TRUE) ~ 3,
    region == "AK" & grepl("(halibut|PacificCod|Arrowtooth)", pred, ignore.case = TRUE) ~ 4,
    region == "AK" & grepl("shark", pred, ignore.case = TRUE) ~ 5,
    TRUE ~ 3
  )
  
  prey_candidates <- get_prey_pool(guild_num, stage_suffix)
  spatial_col     <- get_spatial_overlap_col(pred, region, stage_suffix)
  
  lookup_rows[[length(lookup_rows) + 1]] <- tibble(
    region              = region,
    guild_num           = guild_num,
    pred_data_col       = pred,
    stage_suffix        = stage_suffix,
    spatial_overlap_col = spatial_col,
    all_candidates_str  = paste(c(spatial_col, prey_candidates)[!is.na(c(spatial_col, prey_candidates))], collapse = ", ")
  )
}

guild_pred_lookup <- bind_rows(lookup_rows)
write_csv(guild_pred_lookup, file.path(data_out_dir, paste0("guild_predator_lookup_table_", years_included, ".csv")))

cat("Data preparation, CSV export, and lookup table generation completed successfully.\n")
