

# Phase 1 fitting script (1998-2021, lag 0)
# Uses direct model:
#   x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + I_SSL + I_Shark
#
# Assumes these objects/functions already exist from your build script:
#   ak_raw, sem, shark, clim_tr (optional), out_dir
#   guess_year_col(), find_col(), safe_scale()

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(glue)
})

if (!exists("out_dir")) out_dir <- "copilot/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- safety helpers ----------
if (!exists("safe_scale")) {
  safe_scale <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    as.numeric(scale(x))
  }
}
if (!exists("guess_year_col")) {
  guess_year_col <- function(df) {
    ycols <- names(df)[stringr::str_detect(names(df), "^(year|yr)$|year")]
    if (length(ycols) == 0) return(NA_character_)
    ycols[1]
  }
}
if (!exists("find_col")) {
  find_col <- function(df, patterns) {
    nm <- names(df)
    nm[stringr::str_detect(nm, stringr::regex(paste(patterns, collapse="|"), ignore_case = TRUE))]
  }
}

# ---------------------------
# 1) Build annual source tables--------
# ---------------------------
year_col_ak <- guess_year_col(ak_raw)
if (is.na(year_col_ak)) stop("Could not infer year in ak_raw.")
ak_yr <- ak_raw %>%
  rename(year = all_of(year_col_ak)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

year_col_sem <- guess_year_col(sem)
if (is.na(year_col_sem)) stop("Could not infer year in sem.")
sem_yr <- sem %>%
  rename(year = all_of(year_col_sem)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

year_col_shark <- guess_year_col(shark)
if (is.na(year_col_shark)) stop("Could not infer year in shark.")
shark_yr <- shark %>%
  rename(year = all_of(year_col_shark)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

if (exists("clim_tr") && nrow(clim_tr) > 0) {
  yc <- guess_year_col(clim_tr)
  if (!is.na(yc)) clim_tr <- clim_tr %>% rename(year = all_of(yc))
}

# ---------------------------
# 2) Candidate sets-----------
# ---------------------------
sst_cands <- find_col(ak_yr, c("sst","spring_sst","summer_sst","sea_surface"))
ssl_cands <- find_col(ak_yr, c("steller","sea_lion","ssl"))

# user-vetted forage logic
forage_ak  <- find_col(ak_raw, c("capelin", "sand_lance", "sand lance", "ammod", "herring", "herr"))
forage_sem <- find_col(sem %>% select(!contains("x05")), c("capelin", "sand_lance", "sand lance", "ammod", "herring", "forage"))

# keep only those that exist in annual tables
forage_ak_use  <- intersect(forage_ak, names(ak_yr))
forage_sem_use <- intersect(forage_sem, names(sem_yr))

shark_cands <- find_col(shark_yr, c("shark","biomass","cpue","salmon_shark"))
if (length(shark_cands) == 0) shark_cands <- find_col(sem_yr, c("shark","salmon_shark"))

# strict salmon fields
col_salmon_juv   <- "x07_dfa_cpue_int_spr_jun_hw"
col_salmon_adult <- "x16_sar"

if (!all(c(col_salmon_juv, col_salmon_adult) %in% names(sem_yr))) {
  stop("Required salmon columns not found in sem_yr: x07_dfa_cpue_int_spr_jun_hw and/or x16_sar.")
}
if (length(sst_cands)==0 || length(ssl_cands)==0 || length(shark_cands)==0) {
  stop("One or more candidate sets are empty (sst/ssl/shark). Check naming patterns.")
}
if (length(forage_ak_use) + length(forage_sem_use) == 0) {
  stop("No forage columns available after forage filtering.")
}

message("Candidate counts: ",
        "SST=", length(sst_cands),
        ", SSL=", length(ssl_cands),
        ", Shark=", length(shark_cands),
        ", ForageAK=", length(forage_ak_use),
        ", ForageSEM=", length(forage_sem_use))

# ---------------------------
# 3) Core fit function----------
# ---------------------------
fit_one_combo <- function(col_sst, col_ssl, col_shark,
                          k0 = 0.0,
                          kT = 1.0,  # warm effect on switching (+)
                          kF = 1.0,  # forage reduces switching (used as -kF*F_t)
                          Q10 = 2.0,
                          use_climate = FALSE) {
  
  # build selected data
  data_selected <- tibble(year = sort(unique(c(ak_yr$year, sem_yr$year, shark_yr$year)))) %>%
    left_join(
      ak_yr %>%
        select(
          year,
          SST_raw = all_of(col_sst),
          SSL_raw = all_of(col_ssl),
          all_of(forage_ak_use)
        ),
      by = "year"
    ) %>%
    left_join(
      sem_yr %>%
        select(
          year,
          all_of(col_salmon_juv),   # NCC baseline covariate
          all_of(col_salmon_adult), # SAR response
          all_of(forage_sem_use)
        ),
      by = "year"
    ) %>%
    left_join(
      shark_yr %>%
        select(year, Shark_raw = all_of(col_shark)),
      by = "year"
    )
  
  # forage composite
  forage_cols_all <- c(forage_ak_use, forage_sem_use)
  forage_cols_all <- forage_cols_all[forage_cols_all %in% names(data_selected)]
  
  if (length(forage_cols_all) == 0) {
    return(tibble(
      sst_col = col_sst, ssl_col = col_ssl, shark_col = col_shark,
      forcing = ifelse(use_climate, "T_clim", "T_sst"),
      n = 0,
      aic = NA_real_, r2 = NA_real_,
      b_cpue = NA_real_, b_ssl = NA_real_, b_shark = NA_real_,
      p_cpue = NA_real_, p_ssl = NA_real_, p_shark = NA_real_
    ))
  }
  
  forage_scaled <- data_selected %>%
    select(all_of(forage_cols_all)) %>%
    mutate(across(everything(), safe_scale))
  
  data_selected <- data_selected %>%
    mutate(F_raw = rowMeans(as.matrix(forage_scaled), na.rm = TRUE))
  
  # optional climate forcing
  data_selected$T_clim <- NA_real_
  if (use_climate && exists("clim_tr") && nrow(clim_tr) > 0 && "year" %in% names(clim_tr)) {
    num_cols <- names(clim_tr)[sapply(clim_tr, is.numeric)]
    num_cols <- setdiff(num_cols, "year")
    if (length(num_cols) > 0) {
      data_selected <- data_selected %>%
        left_join(clim_tr %>% select(year, clim_raw = all_of(num_cols[1])), by = "year") %>%
        mutate(T_clim = safe_scale(clim_raw))
    }
  }
  
  # mechanistic construction
  T_ref <- mean(data_selected$SST_raw, na.rm = TRUE)
  
  data_selected <- data_selected %>%
    mutate(
      T_sst   = safe_scale(SST_raw),
      T_use   = if (use_climate && !all(is.na(T_clim))) T_clim else T_sst,
      SSL_z   = safe_scale(SSL_raw),
      Shark_z = safe_scale(Shark_raw),
      F_z     = safe_scale(F_raw)
    )
  
  # SSL switching: no salmon abundance term
  # p_switch increases with T_use and decreases with forage F_z
  p_switch <- plogis(k0 + kT * data_selected$T_use - kF * data_selected$F_z)
  
  C_t   <- rep(1, nrow(data_selected))
  I_SSL <- data_selected$SSL_z * p_switch * C_t
  
  # Shark metabolic index
  O_t     <- pmax(0.2, 1 + 0.3 * data_selected$T_use)
  M_t     <- Q10^((data_selected$SST_raw - T_ref) / 10)
  I_Shark <- data_selected$Shark_z * M_t * O_t
  
  # fit window
  fit_dat <- data_selected %>%
    mutate(
      I_SSL = I_SSL,
      I_Shark = I_Shark
    ) %>%
    filter(year >= 1998, year <= 2021) %>%
    filter(
      !is.na(.data[[col_salmon_adult]]),
      !is.na(.data[[col_salmon_juv]]),
      !is.na(I_SSL),
      !is.na(I_Shark)
    )
  
  if (nrow(fit_dat) < 12) {
    return(tibble(
      sst_col = col_sst, ssl_col = col_ssl, shark_col = col_shark,
      forcing = ifelse(use_climate, "T_clim", "T_sst"),
      n = nrow(fit_dat),
      aic = NA_real_, r2 = NA_real_,
      b_cpue = NA_real_, b_ssl = NA_real_, b_shark = NA_real_,
      p_cpue = NA_real_, p_ssl = NA_real_, p_shark = NA_real_
    ))
  }
  
  # direct model per user request
  mod <- lm(
    reformulate(c(col_salmon_juv, "I_SSL", "I_Shark"), response = col_salmon_adult),
    data = fit_dat
  )
  
  sm <- summary(mod)
  td <- broom::tidy(mod)
  
  get_term <- function(term, field) {
    out <- td %>% filter(.data$term == !!term) %>% pull(!!sym(field))
    if (length(out) == 0) NA_real_ else out[1]
  }
  
  tibble(
    sst_col = col_sst,
    ssl_col = col_ssl,
    shark_col = col_shark,
    forcing = ifelse(use_climate, "T_clim", "T_sst"),
    n = nrow(fit_dat),
    aic = AIC(mod),
    r2 = sm$r.squared,
    b_cpue  = get_term(col_salmon_juv, "estimate"),
    b_ssl   = get_term("I_SSL", "estimate"),
    b_shark = get_term("I_Shark", "estimate"),
    p_cpue  = get_term(col_salmon_juv, "p.value"),
    p_ssl   = get_term("I_SSL", "p.value"),
    p_shark = get_term("I_Shark", "p.value")
  )
}

# ---------------------------
# 4) Grid search over combos
# ---------------------------
grid <- tidyr::crossing(
  sst_col = sst_cands,
  ssl_col = ssl_cands,
  shark_col = shark_cands,
  forcing = c("T_sst", "T_clim")
)

results <- purrr::pmap_dfr(
  grid,
  function(sst_col, ssl_col, shark_col, forcing) {
    fit_one_combo(
      col_sst = sst_col,
      col_ssl = ssl_col,
      col_shark = shark_col,
      k0 = 0.0,
      kT = 1.0,
      kF = 1.0,
      Q10 = 2.0,
      use_climate = (forcing == "T_clim")
    )
  }
)

# ranking (customize sign expectations as desired)
results_ranked <- results %>%
  arrange(aic, desc(r2))

write_csv(results,        file.path(out_dir, "phase1_fit_all_combos.csv"))
write_csv(results_ranked, file.path(out_dir, "phase1_fit_ranked_combos.csv"))
write_csv(slice_head(results_ranked, n = 10), file.path(out_dir, "phase1_fit_top10.csv"))

print(slice_head(results_ranked, n = 10))

# ---------------------------
# 5) Rebuild and save best combo index series
# ---------------------------
best <- results_ranked %>% filter(!is.na(aic)) %>% slice(1)

if (nrow(best) == 1) {
  col_sst   <- best$sst_col[[1]]
  col_ssl   <- best$ssl_col[[1]]
  col_shark <- best$shark_col[[1]]
  use_climate <- best$forcing[[1]] == "T_clim"
  
  data_selected <- tibble(year = sort(unique(c(ak_yr$year, sem_yr$year, shark_yr$year)))) %>%
    left_join(
      ak_yr %>% select(year, SST_raw = all_of(col_sst), SSL_raw = all_of(col_ssl), all_of(forage_ak_use)),
      by = "year"
    ) %>%
    left_join(
      sem_yr %>% select(year, all_of(col_salmon_juv), all_of(col_salmon_adult), all_of(forage_sem_use)),
      by = "year"
    ) %>%
    left_join(
      shark_yr %>% select(year, Shark_raw = all_of(col_shark)),
      by = "year"
    )
  
  forage_cols_all <- c(forage_ak_use, forage_sem_use)
  forage_cols_all <- forage_cols_all[forage_cols_all %in% names(data_selected)]
  
  forage_scaled <- data_selected %>%
    select(all_of(forage_cols_all)) %>%
    mutate(across(everything(), safe_scale))
  
  data_selected <- data_selected %>%
    mutate(F_raw = rowMeans(as.matrix(forage_scaled), na.rm = TRUE))
  
  data_selected$T_clim <- NA_real_
  if (use_climate && exists("clim_tr") && nrow(clim_tr) > 0 && "year" %in% names(clim_tr)) {
    num_cols <- names(clim_tr)[sapply(clim_tr, is.numeric)]
    num_cols <- setdiff(num_cols, "year")
    if (length(num_cols) > 0) {
      data_selected <- data_selected %>%
        left_join(clim_tr %>% select(year, clim_raw = all_of(num_cols[1])), by = "year") %>%
        mutate(T_clim = safe_scale(clim_raw))
    }
  }
  
  # same parameterization as fit function
  k0 <- 0.0
  kT <- 1.0
  kF <- 1.0
  Q10 <- 2.0
  
  T_ref <- mean(data_selected$SST_raw, na.rm = TRUE)
  
  data_selected <- data_selected %>%
    mutate(
      T_sst   = safe_scale(SST_raw),
      T_use   = if (use_climate && !all(is.na(T_clim))) T_clim else T_sst,
      SSL_z   = safe_scale(SSL_raw),
      Shark_z = safe_scale(Shark_raw),
      F_z     = safe_scale(F_raw)
    )
  
  p_switch <- plogis(k0 + kT * data_selected$T_use - kF * data_selected$F_z)
  C_t <- rep(1, nrow(data_selected))
  I_SSL <- data_selected$SSL_z * p_switch * C_t
  
  O_t <- pmax(0.2, 1 + 0.3 * data_selected$T_use)
  M_t <- Q10^((data_selected$SST_raw - T_ref)/10)
  I_Shark <- data_selected$Shark_z * M_t * O_t
  
  I_PredAK <- safe_scale(I_SSL) + safe_scale(I_Shark)
  
  best_idx <- data_selected %>%
    transmute(
      year,
      SST_raw,
      T_sst,
      T_clim,
      T_use,
      F_raw,
      SSL_raw,
      Shark_raw,
      p_switch_ssl = p_switch,
      I_SSL,
      M_shark = M_t,
      I_Shark,
      I_PredAK
    )
  
  write_csv(best_idx, file.path(out_dir, "phase1_best_combo_indices.csv"))
  
  sem_plus_best <- sem_yr %>% left_join(best_idx, by = "year")
  write_csv(sem_plus_best, file.path(out_dir, "sem_master_data_plus_best_phase1_indices.csv"))
}

message("Fitting complete.\n",
        "- ", file.path(out_dir, "phase1_fit_all_combos.csv"), "\n",
        "- ", file.path(out_dir, "phase1_fit_ranked_combos.csv"), "\n",
        "- ", file.path(out_dir, "phase1_fit_top10.csv"), "\n",
        "- ", file.path(out_dir, "phase1_best_combo_indices.csv (if fit succeeded)\n",
                        "- ", file.path(out_dir, "sem_master_data_plus_best_phase1_indices.csv (if fit succeeded)"))