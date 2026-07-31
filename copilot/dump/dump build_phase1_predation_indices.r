# Phase 1 Alaska predator/prey mechanistic indices
# Author: Copilot draft for LisaCrozier-NOAA
# Purpose:
#   1) Harmonize indicator naming across mapping files
#   2) Build warm-state index W_t (SST and optional climate DFA trend)
#   3) Build SSL Type III prey-switching index
#   4) Build salmon shark metabolic predation index
#   5) Export SEM-ready candidate AK predation index table
#
# NOTE:
# - This script is intentionally transparent and parameterized.
# - You will likely tweak column selections once exact names are confirmed.

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(readr)
  library(stringr)
})

# ---------------------------
# 0) Paths
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

out_dir          <- "copilot/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# 1) Read + clean
# ---------------------------
guilds <- read_csv(path_guilds, show_col_types = FALSE) %>% clean_names()
inds   <- read_csv(path_indicators, show_col_types = FALSE) %>% clean_names()
lookup <- read_csv(path_lookup, show_col_types = FALSE) %>% clean_names()
sem    <- read_csv(path_sem, show_col_types = FALSE) %>% clean_names()
sharkGOA    <- read_csv(path_shark, show_col_types = FALSE) %>% clean_names() %>%   filter(fmp=="GOA") 
sharkBSAI    <- read_csv(path_shark, show_col_types = FALSE) %>% clean_names() %>%   filter(fmp=="BSAI")
shark <- left_join(sharkGOA,sharkBSAI, join_by(year),suffix = c("_goa","_bsai"))


egoa   <- read_csv(path_egoa, show_col_types = FALSE) %>% clean_names()
wgoa   <- read_csv(path_wgoa, show_col_types = FALSE) %>% clean_names()

egoa_m <- read_csv(path_egoa_meta, show_col_types = FALSE) %>% clean_names()
wgoa_m <- read_csv(path_wgoa_meta, show_col_types = FALSE) %>% clean_names()

clim_tr <- if (file.exists(path_clim_trend)) {
  read_csv(path_clim_trend, show_col_types = FALSE) %>% clean_names()
} else {
  tibble()
}

# ---------------------------
# 2) Basic validation helpers
# ---------------------------
stop_if_missing <- function(df, nm, required_cols) {
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0) {
    stop(glue::glue("{nm} missing required columns: {paste(miss, collapse=', ')}"))
  }
}

# We only strictly require short_name for joins.
# If your file uses shortName and janitor made short_name, this should work.
stop_if_missing(guilds, "guilds", c("short_name"))
stop_if_missing(inds,   "indicators", c("short_name"))

# Try to infer year column from each frame
guess_year_col <- function(df) {
  ycols <- names(df)[str_detect(names(df), "^(year|yr)$|year")]
  if (length(ycols) == 0) return(NA_character_)
  ycols[1]
}

# ---------------------------
# 3) Build canonical indicator dictionary
# ---------------------------
# Keep key SEM structure info from guilds
dict <- guilds %>%
  select(any_of(c("short_name", "semnode", "guild"))) %>%
  distinct()

# Add indicator metadata
dict <- dict %>%
  left_join(inds, by = "short_name", suffix = c("_guild", "_ind"))

# Add alias table if it has a recognizable short_name or alias field
# (Flexible join attempt)
lookup_cols <- names(lookup)
if (!"short_name" %in% lookup_cols) {
  # try common candidates
  cand <- c("shortname", "short_name_model", "short_name_sem")
  hit <- cand[cand %in% lookup_cols]
  if (length(hit) > 0) {
    lookup <- lookup %>% rename(short_name = all_of(hit[1]))
  }
}
if ("short_name" %in% names(lookup)) {
  dict <- dict %>% left_join(lookup, by = "short_name")
}

write_csv(dict, file.path(out_dir, "phase1_indicator_dictionary.csv"))

# ---------------------------
# 4) Build raw AK annual table
# ---------------------------
# Add region tags, bind
egoa2 <- egoa %>% mutate(region = "EGOA")
wgoa2 <- wgoa %>% mutate(region = "WGOA")
ak_raw <- bind_rows(egoa2, wgoa2)

year_col_ak <- guess_year_col(ak_raw)
if (is.na(year_col_ak)) stop("Could not infer year column in GOA raw data.")

ak_raw <- ak_raw %>% rename(year = all_of(year_col_ak))

# ---------------------------
# 5) Identify candidate columns by pattern (edit as needed)
# ---------------------------
find_col <- function(df, patterns) {
  nm <- names(df)
  hits <- nm[str_detect(nm, regex(paste(patterns, collapse="|"), ignore_case = TRUE))]
  hits
}

# SST candidates
sst_cands <- find_col(ak_raw, c("sst", "sea_surface_temp", "spring_sst", "summer_sst"))

# Forage (capelin, sand lance, herring as candidate preferred forage components)
forage_cands1 <- find_col(ak_raw, c("capelin", "sand_lance", "sand lance","ammod", "herring", "herr"))
forage_cands2 <- find_col(sem %>% select(!contains("x05")), c("capelin", "sand_lance", "sand lance","ammod", "herring", "forage")) 
forage_cands <- c(forage_cands1,forage_cands2)

# SSL abundance / trend
ssl_cands <- find_col(ak_raw, c("steller", "sea_lion", "sea lion", "ssl"))

# Salmon shark abundance proxy
shark_cands1 <- find_col(shark, c("salmon_shark", "salmon shark", "shark"))
shark_cands2 <- find_col(sem, c("salmon_shark", "salmon shark", "shark"))
shark_cands <-c(shark_cands1,shark_cands2)

# Salmon availability proxy in AK context
salmon_cands <- find_col(sem, c("x07", "x16"))

message("Candidate columns:\n",
        "SST: ", paste(sst_cands, collapse=", "), "\n",
        "Forage: ", paste(forage_cands, collapse=", "), "\n",
        "SSL: ", paste(ssl_cands, collapse=", "), "\n",
        "Shark: ", paste(shark_cands1, collapse=", "), "\n",
        "Salmon: ", paste(salmon_cands, collapse=", "))

# ---------------------------
# 6) USER-EDITABLE COLUMN SELECTION
# ---------------------------
# Pick one or more exact column names after first run messages.
# For now choose first hit as default.
pick1 <- function(x) ifelse(length(x) > 0, x[1], NA_character_)

col_sst      <- pick1(sst_cands)
col_ssl      <- pick1(ssl_cands)
col_shark    <- pick1(shark_cands)
col_salmon   <- pick1(salmon_cands)

# forage: allow multiple and average after z-scaling
cols_forage  <- forage_cands

required_selected <- c(col_sst, col_ssl, col_shark, col_salmon)
if (any(is.na(required_selected)) || length(cols_forage) == 0) {
  stop("Need to manually set selected columns (col_sst, col_ssl, col_shark, col_salmon, cols_forage).")
}

# ---------------------------
# 7) Standardization utilities
# ---------------------------
z <- function(x) as.numeric(scale(x))

# Optional 3-pt running mean (centered)
smooth3 <- function(x) {
  stats::filter(x, rep(1/3, 3), sides = 2) %>% as.numeric()
}

# ---------------------------
# 8) Construct annual mechanistic components
# ---------------------------
ak <- ak_raw %>%
  group_by(year) %>%
  summarise(
    SST_raw      = mean(.data[[col_sst]], na.rm = TRUE),
    SSL_raw      = mean(.data[[col_ssl]], na.rm = TRUE),
    Shark_raw    = mean(.data[[col_shark]], na.rm = TRUE),
    Salmon_raw   = mean(.data[[col_salmon]], na.rm = TRUE),
    across(all_of(cols_forage), ~mean(.x, na.rm = TRUE), .names = "{.col}"),
    .groups = "drop"
  )

# Preferred forage index: mean z-score across selected forage components
forage_mat <- ak %>% select(all_of(cols_forage)) %>% as.matrix()
F_raw <- rowMeans(scale(forage_mat), na.rm = TRUE)

ak <- ak %>%
  mutate(
    F_raw = as.numeric(F_raw),
    
    # Warm-state index W_t (higher = warmer)
    W_sst = z(SST_raw),
    
    # Optional smoothing pass (disabled by default for direct transparency)
    W_sst_sm = W_sst,
    F_sm     = z(F_raw),
    SSL_sm   = z(SSL_raw),
    Shark_sm = z(Shark_raw),
    Salmon_sm= z(Salmon_raw)
  )

# Optional climate trend integration
if (nrow(clim_tr) > 0) {
  yc <- guess_year_col(clim_tr)
  if (!is.na(yc)) {
    clim_tr2 <- clim_tr %>% rename(year = all_of(yc))
    # pick first numeric trend-like column (manual override may be needed)
    num_cols <- names(clim_tr2)[sapply(clim_tr2, is.numeric)]
    num_cols <- setdiff(num_cols, "year")
    if (length(num_cols) > 0) {
      col_clim <- num_cols[1]
      ak <- ak %>%
        left_join(clim_tr2 %>% select(year, clim_raw = all_of(col_clim)), by = "year") %>%
        mutate(W_clim = z(clim_raw))
    }
  }
}

# ---------------------------
# 9) Mechanistic index formulas
# ---------------------------
# Parameter set (edit in sensitivity runs)
m       <- 2.0
alpha   <- 1.0
beta    <- 0.8
Q10     <- 2.0
T_ref   <- mean(ak$SST_raw, na.rm = TRUE)

# Use SST-based warm state as primary
W_use <- ak$W_sst_sm

# Compression/overlap placeholders (can replace with better proxies later)
C_t <- rep(1, nrow(ak))                  # SSL compression modifier
O_t <- pmax(0.2, 1 + 0.3 * W_use)        # Shark overlap modifier, bounded

# SSL Type III
Nsalmon <- pmax(1e-6, ak$Salmon_sm - min(ak$Salmon_sm, na.rm = TRUE) + 1)
F_t     <- pmax(1e-6, ak$F_sm     - min(ak$F_sm, na.rm = TRUE) + 1)

p_salmon <- (Nsalmon^m) / ((Nsalmon^m) + (alpha * F_t * exp(-beta * W_use))^m)
I_SSL    <- ak$SSL_sm * p_salmon * C_t

# Shark metabolic index
M_t      <- Q10^((ak$SST_raw - T_ref)/10)
I_Shark  <- ak$Shark_sm * M_t * O_t

# Composite PredAK candidate
I_PredAK <- z(I_SSL) + z(I_Shark)

ak_idx <- ak %>%
  transmute(
    year,
    W_sst,
    W_sst_sm,
    W_clim = if ("W_clim" %in% names(.)) W_clim else NA_real_,
    F_raw,
    SSL_raw,
    Shark_raw,
    Salmon_raw,
    p_salmon_ssl = p_salmon,
    I_SSL = I_SSL,
    M_shark = M_t,
    I_Shark = I_Shark,
    I_PredAK = I_PredAK
  )

write_csv(ak_idx, file.path(out_dir, "phase1_AK_mechanistic_indices.csv"))

# ---------------------------
# 10) Join with sem_master_data for modeling handoff
# ---------------------------
year_col_sem <- guess_year_col(sem)
if (is.na(year_col_sem)) stop("Could not infer year column in sem_master_data.csv")
sem2 <- sem %>% rename(year = all_of(year_col_sem))

sem_plus <- sem2 %>%
  left_join(ak_idx, by = "year")

write_csv(sem_plus, file.path(out_dir, "sem_master_data_plus_phase1_indices.csv"))

# ---------------------------
# 11) Minimal diagnostics
# ---------------------------
diag_long <- ak_idx %>%
  select(year, W_sst, p_salmon_ssl, I_SSL, I_Shark, I_PredAK) %>%
  pivot_longer(-year, names_to = "series", values_to = "value")

p <- ggplot(diag_long, aes(year, value, color = series)) +
  geom_line(linewidth = 0.8) +
  theme_bw() +
  labs(title = "Phase 1 Mechanistic Index Diagnostics", x = "Year", y = "Value")

ggsave(file.path(out_dir, "phase1_index_diagnostics.png"), p, width = 10, height = 6, dpi = 150)

message("Done. Wrote:\n",
        "- ", file.path(out_dir, "phase1_indicator_dictionary.csv"), "\n",
        "- ", file.path(out_dir, "phase1_AK_mechanistic_indices.csv"), "\n",
        "- ", file.path(out_dir, "sem_master_data_plus_phase1_indices.csv"), "\n",
        "- ", file.path(out_dir, "phase1_index_diagnostics.png"))