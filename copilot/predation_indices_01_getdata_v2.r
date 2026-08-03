
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(readr)
  library(stringr)
  library(glue)
})

# ---------------------------
# 0) Paths---------
# ---------------------------
path_guilds      <- "Doug_code/guildfile.excludecol.csv"
path_indicators  <- "Doug_code/indicators.csv"
path_lookup      <- "outputs_4/var_lookup_NCC_AK.SAR.csv"
path_sem         <- "outputs_4/sem_master_data.csv"
path_shark       <- "data_Lisa/AKshark.table19.3_2022assess.csv"

path_egoa        <- "Ferris-DFAIndicators-goa/data/EGOA_EcoState_Data_Jan2023.csv"
path_wgoa        <- "Ferris-DFAIndicators-goa/data/WGOA_EcoState_Data_Jan2023.csv"
path_egoa_meta   <- "Ferris-DFAIndicators-goa/data/EGOA_metadata.csv"
path_wgoa_meta   <- "Ferris-DFAIndicators-goa/data/WGOA_metadata.csv"
path_clim_trend  <- "Ferris-DFAIndicators-goa/MARSS results/climate_trends.csv"

out_dir          <- "copilot/outputs_2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# 1) Read + clean--------
# ---------------------------
guilds <- read_csv(path_guilds, show_col_types = FALSE) %>% clean_names()
inds   <- read_csv(path_indicators, show_col_types = FALSE) %>% clean_names()
lookup <- read_csv(path_lookup, show_col_types = FALSE) %>% clean_names()
sem    <- read_csv(path_sem, show_col_types = FALSE) %>% clean_names()

shark_raw <- read_csv(path_shark, show_col_types = FALSE) %>% clean_names()
shark_goa <- shark_raw %>% filter(fmp == "GOA")
shark_bsai <- shark_raw %>% filter(fmp == "BSAI")
shark <- left_join(shark_goa, shark_bsai, by = join_by(year), suffix = c("_goa","_bsai"))

egoa   <- read_csv(path_egoa, show_col_types = FALSE) %>% clean_names()
wgoa   <- read_csv(path_wgoa, show_col_types = FALSE) %>% clean_names()

# Metadata loaded for context/possible future use (not required in core math below)
egoa_m <- read_csv(path_egoa_meta, show_col_types = FALSE) %>% clean_names()
wgoa_m <- read_csv(path_wgoa_meta, show_col_types = FALSE) %>% clean_names()

clim_tr <- if (file.exists(path_clim_trend)) {
  read_csv(path_clim_trend, show_col_types = FALSE) %>% clean_names()
} else {
  tibble()
}

# ---------------------------
# 2) Helpers----
# ---------------------------
stop_if_missing <- function(df, nm, required_cols) {
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0) stop(glue("{nm} missing required columns: {paste(miss, collapse=', ')}"))
}

guess_year_col <- function(df) {
  ycols <- names(df)[str_detect(names(df), "^(year|yr)$|year")]
  if (length(ycols) == 0) return(NA_character_)
  ycols[1]
}

find_col <- function(df, patterns) {
  nm <- names(df)
  nm[str_detect(nm, regex(paste(patterns, collapse="|"), ignore_case = TRUE))]
}

pick1 <- function(x) ifelse(length(x) > 0, x[1], NA_character_)

z <- function(x) as.numeric(scale(x))

safe_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

# ---------------------------
# 4) Build annual source tables-------
# ---------------------------
# GOA
egoa2 <- egoa %>% mutate(region = "EGOA")
wgoa2 <- wgoa %>% mutate(region = "WGOA")
ak_raw <- bind_rows(egoa2, wgoa2)

year_col_ak <- guess_year_col(ak_raw)
if (is.na(year_col_ak)) stop("Could not infer year column in GOA raw data.")
ak_yr <- ak_raw %>%
  rename(year = all_of(year_col_ak)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# SEM
year_col_sem <- guess_year_col(sem)
if (is.na(year_col_sem)) stop("Could not infer year column in sem_master_data.csv")
sem_yr <- sem %>%
  rename(year = all_of(year_col_sem)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# Shark
year_col_shark <- guess_year_col(shark)
if (is.na(year_col_shark)) stop("Could not infer year column in shark table.")
shark_yr <- shark %>%
  rename(year = all_of(year_col_shark)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# Climate trend optional
if (nrow(clim_tr) > 0) {
  yc <- guess_year_col(clim_tr)
  if (!is.na(yc)) clim_tr <- clim_tr %>% rename(year = all_of(yc))
}

# ---------------------------
# 5) Candidate columns (auto)---------
# ---------------------------
sst_cands   <- find_col(ak_yr, c("sst","spring_sst","summer_sst","sea_surface"))
ssl_cands   <- find_col(ak_yr, c("steller","sea_lion","ssl"))

# forage from GOA + SEM (user-vetted logic)
forage_ak <- find_col(ak_raw, c("capelin", "sand_lance", "sand lance", "ammod", "herring", "herr"))
forage_sem <- find_col(sem %>% select(!contains("x05")), c("capelin", "sand_lance", "sand lance", "ammod", "herring", "forage"))
forage_cands  <- unique(c(forage_ak, forage_sem))


shark_cands_shark <- find_col(shark_yr, c("shark","biomass","cpue","salmon_shark"))
shark_cands_sem   <- find_col(sem_yr, c("shark","salmon_shark"))
shark_cands <- unique(c(shark_cands_shark, shark_cands_sem))

# Strict salmon candidates per user: x07_* and x16_*
salmon_cands <- find_col(sem_yr, c("^x07_", "^x16_"))

message("Candidates by source:\n",
        "SST (ak_yr): ", paste(sst_cands, collapse=", "), "\n",
        "SSL (ak_yr): ", paste(ssl_cands, collapse=", "), "\n",
        "Forage (ak+sem): ", paste(forage_cands, collapse=", "), "\n",
        "Shark (shark+sem): ", paste(shark_cands, collapse=", "), "\n",
        "Salmon (sem): ", paste(salmon_cands, collapse=", "))

# # ---------------------------
# # 6) Manual picks (edit here if needed)
# # ---------------------------
# col_sst    <- pick1(sst_cands)            # from ak_yr
# col_ssl    <- pick1(ssl_cands)            # from ak_yr
# col_shark  <- pick1(shark_cands_shark)    # prefer shark table

# User-specified salmon columns
col_salmon_juv   <- "x07_dfa_cpue_int_spr_jun_hw"
col_salmon_adult <- "x16_sar"

# Forage columns from both sources
cols_forage_ak  <- intersect(forage_ak,  names(ak_yr))
cols_forage_sem <- intersect(forage_sem, names(sem_yr))

if (any(is.na(c(col_sst, col_ssl, col_shark)))) {
  stop("Set col_sst/col_ssl/col_shark manually based on candidate printout.")
}
if (!all(c(col_salmon_juv, col_salmon_adult) %in% names(sem_yr))) {
  stop("Juvenile/adult salmon columns not found in sem_yr.")
}
if (length(cols_forage_ak) + length(cols_forage_sem) == 0) {
  stop("No forage columns selected from AK or SEM sources.")
}

# ---------------------------
# 7) Unified annual design table
# ---------------------------
design <- tibble(year = sort(unique(c(ak_yr$year, sem_yr$year, shark_yr$year)))) %>%
  left_join(
    ak_yr %>%
      select(
        year,
        SST_raw = all_of(col_sst),
        SSL_raw = all_of(col_ssl),
        all_of(cols_forage_ak)
      ),
    by = "year"
  ) %>%
  left_join(
    sem_yr %>%
      select(
        year,
        all_of(col_salmon_juv),
        all_of(col_salmon_adult),
        all_of(cols_forage_sem)
      ),
    by = "year"
  ) %>%
  left_join(
    shark_yr %>%
      select(
        year,
        Shark_raw = all_of(col_shark)
      ),
    by = "year"
  ) %>%
  mutate(
    Salmon_raw = rowMeans(cbind(.data[[col_salmon_juv]], .data[[col_salmon_adult]]), na.rm = TRUE)
  )

# Build forage composite F_raw
forage_cols_all <- c(cols_forage_ak, cols_forage_sem)
forage_cols_all <- forage_cols_all[forage_cols_all %in% names(design)]

if (length(forage_cols_all) == 0) {
  stop("No forage columns present in design after joins.")
}

forage_scaled <- design %>%
  select(all_of(forage_cols_all)) %>%
  mutate(across(everything(), safe_scale))

design <- design %>%
  mutate(
    F_raw = rowMeans(as.matrix(forage_scaled), na.rm = TRUE)
  )

# Optional climate trend
if (nrow(clim_tr) > 0 && "year" %in% names(clim_tr)) {
  num_cols <- names(clim_tr)[sapply(clim_tr, is.numeric)]
  num_cols <- setdiff(num_cols, "year")
  if (length(num_cols) > 0) {
    col_clim <- num_cols[1]
    design <- design %>%
      left_join(clim_tr %>% select(year, clim_raw = all_of(col_clim)), by = "year") %>%
      mutate(W_clim = safe_scale(clim_raw))
  } else {
    design$W_clim <- NA_real_
  }
} else {
  design$W_clim <- NA_real_
}

