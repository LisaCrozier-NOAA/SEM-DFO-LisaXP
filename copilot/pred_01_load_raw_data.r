

#load raw data
#create ak_yr, sem_yr, shark_yr 

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(readr)
  library(stringr)
})

# =========================================================
# Clean get-data script
# Goal: create annual source tables:
#   - ak_yr
#   - sem_yr
#   - shark_yr
# from raw inputs, with minimal processing.
# =========================================================

# ---------------------------
# 0) Paths
# ---------------------------
path_sem   <- "outputs_4/sem_master_data.csv"
path_shark <- "data_Lisa/AKshark.table19.3_2022assess.csv"

path_egoa  <- "Ferris-DFAIndicators-goa/data/EGOA_EcoState_Data_Jan2023.csv"
path_wgoa  <- "Ferris-DFAIndicators-goa/data/WGOA_EcoState_Data_Jan2023.csv"

out_dir    <- "copilot/outputs_2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# 1) Helpers
# ---------------------------
guess_year_col <- function(df) {
  nms <- names(df)
  hit <- nms[str_detect(nms, regex("^(year|yr)$|year", ignore_case = TRUE))]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

to_annual_numeric <- function(df, df_name = "data") {
  yc <- guess_year_col(df)
  if (is.na(yc)) stop(df_name, ": could not infer a year column.")
  
  df %>%
    rename(year = all_of(yc)) %>%
    mutate(year = suppressWarnings(as.integer(year))) %>%
    filter(!is.na(year)) %>%
    group_by(year) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    arrange(year)
}

safe_read_csv_clean <- function(path, name) {
  if (!file.exists(path)) stop("Missing file: ", path, " (", name, ")")
  read_csv(path, show_col_types = FALSE) %>% clean_names()
}

# ---------------------------
# 2) Read raw data
# ---------------------------
sem_raw   <- safe_read_csv_clean(path_sem, "SEM master")
shark_raw <- safe_read_csv_clean(path_shark, "Shark table")
egoa_raw  <- safe_read_csv_clean(path_egoa, "EGOA indicators")
wgoa_raw  <- safe_read_csv_clean(path_wgoa, "WGOA indicators")

# ---------------------------
# 3) Build sem_yr
# ---------------------------
sem_yr <- to_annual_numeric(sem_raw, "sem_raw")

# ---------------------------
# 4) Build ak_yr (EGOA + WGOA combined, annual)
# ---------------------------
egoa2 <- egoa_raw %>% mutate(region = "egoa")
wgoa2 <- wgoa_raw %>% mutate(region = "wgoa")

ak_raw <- bind_rows(egoa2, wgoa2)
ak_yr  <- to_annual_numeric(ak_raw, "ak_raw")

# ---------------------------
# 5) Build shark_yr (GOA + BSAI side-by-side by year)
# ---------------------------
if (!"fmp" %in% names(shark_raw)) {
  stop("shark_raw is missing 'fmp' column needed to split GOA/BSAI.")
}

shark_goa  <- shark_raw %>% filter(str_to_upper(fmp) == "GOA")
shark_bsai <- shark_raw %>% filter(str_to_upper(fmp) == "BSAI")

# infer year col before join
yg <- guess_year_col(shark_goa)
yb <- guess_year_col(shark_bsai)
if (is.na(yg) || is.na(yb)) stop("Could not infer year column in shark GOA/BSAI subsets.")

shark_joined <- shark_goa %>%
  rename(year = all_of(yg)) %>%
  left_join(
    shark_bsai %>% rename(year = all_of(yb)),
    by = "year",
    suffix = c("_goa", "_bsai")
  ) %>%
  clean_names()

shark_yr <- to_annual_numeric(shark_joined, "shark_joined")

# ---------------------------
# 6) Write outputs
# ---------------------------
write_csv(ak_yr,    file.path(out_dir, "ak_yr.csv"))
write_csv(sem_yr,   file.path(out_dir, "sem_yr.csv"))
write_csv(shark_yr, file.path(out_dir, "shark_yr.csv"))

# Also save RDS for fast reload
saveRDS(ak_yr,    file.path(out_dir, "ak_yr.rds"))
saveRDS(sem_yr,   file.path(out_dir, "sem_yr.rds"))
saveRDS(shark_yr, file.path(out_dir, "shark_yr.rds"))

message("Created annual tables:")
message(" - ak_yr:    ", nrow(ak_yr), " years, ", ncol(ak_yr), " columns")
message(" - sem_yr:   ", nrow(sem_yr), " years, ", ncol(sem_yr), " columns")
message(" - shark_yr: ", nrow(shark_yr), " years, ", ncol(shark_yr), " columns")
message("Wrote outputs to: ", out_dir)

# Keep in environment for downstream scripts:
# ak_yr, sem_yr, shark_yr


#Reduce to relevant columns===============
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(readr)
})

# Assumes these exist from your getdata workflow:
# ak_yr, sem_yr, shark_yr
# (or sem if you prefer raw sem, but sem_yr is cleaner for annual modeling)

required <- c("ak_yr", "sem_yr", "shark_yr")
miss <- required[!vapply(required, exists, logical(1))]
if (length(miss) > 0) {
  stop("Missing objects: ", paste(miss, collapse = ", "),
       ". Run your getdata script first.")
}

if (!exists("out_dir")) out_dir <- "copilot/outputs_2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot("year" %in% names(ak_yr), "year" %in% names(sem_yr), "year" %in% names(shark_yr))

# ---------------------------
# 1) Reconstruct candidate pools exactly like your workflow intent
# ---------------------------
find_col <- function(df, patterns) {
  nm <- names(df)
  nm[str_detect(nm, regex(paste(patterns, collapse="|"), ignore_case = TRUE))]
}

# AK candidates
sst_cands <- find_col(ak_yr, c("sst","spring_sst","summer_sst","sea_surface","enso","pdo","swln_temp"))
ssl_cands <- find_col(ak_yr, c("steller","sea_lion","ssl"))
forage_ak <- find_col(ak_yr, c("capelin","sand_lance","sand lance","ammod","herring","herr","forage"))

# SEM candidates
forage_sem <- find_col(sem_yr %>% select(!contains("x05")), c("capelin","sand_lance","sand lance","ammod","herring","forage"))
salmon_cands <- find_col(sem_yr, c("^x07_","^x16_"))
shark_sem <- find_col(sem_yr, c("shark","salmon_shark"))

# Shark table candidates
shark_cands <- find_col(shark_yr, c("shark","biomass","cpue","salmon_shark","sleeper"))

# Union = all tested families
candidate_cols <- unique(c(
  sst_cands, ssl_cands,
  forage_ak, forage_sem,
  salmon_cands,
  shark_sem, shark_cands
))

# Keep only columns that actually exist in each source
ak_keep <- intersect(candidate_cols, names(ak_yr))
sem_keep <- intersect(candidate_cols, names(sem_yr))
shark_keep <- intersect(candidate_cols, names(shark_yr))

# ---------------------------
# 2) Build annual master by year
# ---------------------------
years <- sort(unique(c(ak_yr$year, sem_yr$year, shark_yr$year)))

master <- tibble(year = years) %>%
  left_join(ak_yr %>% select(any_of(c("year", ak_keep))), by = "year", suffix = c("", "_ak")) %>%
  left_join(sem_yr %>% select(any_of(c("year", sem_keep))), by = "year", suffix = c("", "_sem")) %>%
  left_join(shark_yr %>% select(any_of(c("year", shark_keep))), by = "year", suffix = c("", "_shark"))

# ---------------------------
# 3) Coalesce duplicate names created by joins (x / x_sem / x_shark)
# ---------------------------
# For columns present in >1 table, create one unified column by first non-NA
all_names <- names(master)
base_names <- unique(gsub("(_ak|_sem|_shark)$", "", all_names))

for (b in base_names) {
  variants <- intersect(c(b, paste0(b, "_ak"), paste0(b, "_sem"), paste0(b, "_shark")), names(master))
  if (length(variants) <= 1) next
  
  vals <- master[variants]
  master[[b]] <- apply(vals, 1, function(r) {
    z <- r[!is.na(r)]
    if (length(z) == 0) NA else z[[1]]
  })
}

# Keep clean output: year + unified candidate columns only
master_clean <- master %>%
  select(year, any_of(candidate_cols)) %>%
  arrange(year)

# ---------------------------
# 4) Write outputs
# ---------------------------
write_csv(master_clean, file.path(out_dir, "data_all_tested_columns_annual.csv"))

dict <- tibble(
  column = setdiff(names(master_clean), "year"),
  in_ak_yr = column %in% names(ak_yr),
  in_sem_yr = column %in% names(sem_yr),
  in_shark_yr = column %in% names(shark_yr),
  candidate_group = case_when(
    column %in% sst_cands ~ "sst_or_ocean_temp",
    column %in% ssl_cands ~ "ssl",
    column %in% forage_ak | column %in% forage_sem ~ "forage",
    column %in% salmon_cands ~ "salmon_response_or_covariate",
    column %in% shark_cands | column %in% shark_sem ~ "shark",
    TRUE ~ "other"
  )
)

write_csv(dict, file.path(out_dir, "data_meta_all_tested_columns_dictionary.csv"))

message("Wrote:")
message(" - ", file.path(out_dir, "data_all_tested_columns_annual.csv"))
message(" - ", file.path(out_dir, "data_meta_all_tested_columns_dictionary.csv"))
message("Columns in master: ", ncol(master_clean) - 1)