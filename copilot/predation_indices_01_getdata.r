# Phase 1 Alaska predator/prey mechanistic indices
# Repo: LisaCrozier-NOAA/SEM-DFO-LisaXP
# Date: 2026-07-30
#
# Purpose
#   1) Harmonize indicator naming across mapping files
#   2) Build annual hybrid-source design table (GOA raw + SEM + shark assessment + climate trend)
#   3) Construct mechanistic indices:
#        - SSL Type III prey-switching index
#        - Salmon shark metabolic predation index
#        - Composite PredAK index
#   4) Export SEM-ready outputs for downstream modeling
#
# Notes
#   - Uses strict salmon proxies from SEM: x07_* (juvenile), x16_* (adult/SAR proxy provided by user)
#   - Keeps assumptions transparent and parameterized
#   - If candidate auto-detection is imperfect, set manual picks explicitly

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

out_dir          <- "copilot/outputs"
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
# 3) Build canonical indicator dictionary-----
# ---------------------------
stop_if_missing(guilds, "guilds", c("short_name"))
stop_if_missing(inds,   "indicators", c("short_name"))

dict <- guilds %>%
  select(any_of(c("short_name", "semnode", "guild"))) %>%
  distinct() %>%
  left_join(inds, by = "short_name", suffix = c("_guild", "_ind"))

lookup_cols <- names(lookup)
if (!"short_name" %in% lookup_cols) {
  cand <- c("shortname", "short_name_model", "short_name_sem")
  hit <- cand[cand %in% lookup_cols]
  if (length(hit) > 0) lookup <- lookup %>% rename(short_name = all_of(hit[1]))
}
if ("short_name" %in% names(lookup)) {
  dict <- dict %>% left_join(lookup, by = "short_name")
}

write_csv(dict, file.path(out_dir, "phase1_indicator_dictionary.csv"))

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

# ---------------------------
# 6) Manual picks (edit here if needed)
# ---------------------------
col_sst    <- pick1(sst_cands)            # from ak_yr
col_ssl    <- pick1(ssl_cands)            # from ak_yr
col_shark  <- pick1(shark_cands_shark)    # prefer shark table

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

# ---------------------------
# 8) Mechanistic indices
# ---------------------------
# Parameters (tune in sensitivity runs)
m     <- 2.0
alpha <- 1.0
beta  <- 0.8
Q10   <- 2.0
T_ref <- mean(design$SST_raw, na.rm = TRUE)

design <- design %>%
  mutate(
    W_sst   = safe_scale(SST_raw),
    W_use   = W_sst,   # swap to W_clim for alternate forcing if desired
    SSL_z   = safe_scale(SSL_raw),
    Shark_z = safe_scale(Shark_raw),
    Sal_z   = safe_scale(Salmon_raw),
    F_z     = safe_scale(F_raw)
  )

# Shift scaled salmon/forage to positive support for Hill function
Nsalmon <- pmax(1e-6, design$Sal_z - min(design$Sal_z, na.rm = TRUE) + 1)
F_t     <- pmax(1e-6, design$F_z   - min(design$F_z, na.rm = TRUE) + 1)

# Modifiers (simple placeholders for now)
C_t <- rep(1, nrow(design))                  # SSL compression modifier
O_t <- pmax(0.2, 1 + 0.3 * design$W_use)     # shark overlap modifier

# SSL Type III switching index
p_salmon <- (Nsalmon^m) / ((Nsalmon^m) + (alpha * F_t * exp(-beta * design$W_use))^m)
I_SSL    <- design$SSL_z * p_salmon * C_t

# Shark metabolic predation index
M_t      <- Q10^((design$SST_raw - T_ref)/10)
I_Shark  <- design$Shark_z * M_t * O_t

# Composite PredAK candidate
I_PredAK <- safe_scale(I_SSL) + safe_scale(I_Shark)

ak_idx <- design %>%
  transmute(
    year,
    SST_raw,
    W_sst,
    W_clim,
    W_use,
    F_raw,
    SSL_raw,
    Shark_raw,
    Salmon_raw,
    p_salmon_ssl = p_salmon,
    I_SSL,
    M_shark = M_t,
    I_Shark,
    I_PredAK
  )

# ---------------------------
# 9) Save outputs
# ---------------------------
write_csv(ak_idx, file.path(out_dir, "phase1_AK_mechanistic_indices.csv"))

sem_plus <- sem_yr %>% left_join(ak_idx, by = "year")
write_csv(sem_plus, file.path(out_dir, "sem_master_data_plus_phase1_indices.csv"))

# NA diagnostics
na_diag <- tibble(
  column = names(ak_idx),
  n_na   = map_int(ak_idx, ~sum(is.na(.x))),
  pct_na = round(100 * n_na / nrow(ak_idx), 2)
)
write_csv(na_diag, file.path(out_dir, "phase1_index_na_diagnostics.csv"))

# simple plot
diag_long <- ak_idx %>%
  select(year, W_sst, p_salmon_ssl, I_SSL, I_Shark, I_PredAK) %>%
  pivot_longer(-year, names_to = "series", values_to = "value")

p <- ggplot(diag_long, aes(year, value, color = series)) +
  geom_line(linewidth = 0.8) +
  theme_bw() +
  labs(title = "Phase 1 Mechanistic Index Diagnostics", x = "Year", y = "Value")

ggsave(file.path(out_dir, "phase1_index_diagnostics.png"), p, width = 10, height = 6, dpi = 150)

message("Done. Files written to: ", out_dir, "\n",
        "- phase1_indicator_dictionary.csv\n",
        "- phase1_AK_mechanistic_indices.csv\n",
        "- sem_master_data_plus_phase1_indices.csv\n",
        "- phase1_index_na_diagnostics.csv\n",
        "- phase1_index_diagnostics.png")
