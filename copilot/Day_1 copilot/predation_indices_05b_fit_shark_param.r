#Explore best temperature for sharks and other model decisions

#results:
shark_model<-imp_params_restrict %>% filter(importance_rel_top1>0.8)
print(shark_model,n=Inf)

# model_group shark_col                  importance_aicw param           level                      shark_temp_col            Q10 overlap_form overlap_slope shark_transform  rank importance_rel_top1
# 1 long        pacific_sleeper_shark_goa            0.518 shark_col       pacific_sleeper_shark_goa  NA                         NA NA                    NA   NA                  1               1    
# 2 long        salmon_shark_goa                     0.482 shark_col       salmon_shark_goa           NA                         NA NA                    NA   NA                  2               0.932
# 3 short       pacific_sleeper_shark_bsai           0.541 shark_col       pacific_sleeper_shark_bsai NA                         NA NA                    NA   NA                  1               1    
# 4 short       salmon_shark_bsai                    0.459 shark_col       salmon_shark_bsai          NA                         NA NA                    NA   NA                  2               0.849
# 5 long        NA                                   0.366 shark_temp_col  enso_dj                    enso_dj                    NA NA                    NA   NA                  1               1    
# 6 long        NA                                   0.297 shark_temp_col  pdo_djf                    pdo_djf                    NA NA                    NA   NA                  3               0.813
# 7 long        NA                                   0.337 shark_temp_col  swln_temp_spr_176to226m    swln_temp_spr_176to226m    NA NA                    NA   NA                  2               0.922
# 8 short       NA                                   0.372 shark_temp_col  pdo_djf                    pdo_djf                    NA NA                    NA   NA                  1               1    
# 9 short       NA                                   0.331 shark_temp_col  swln_temp_spr_176to226m    swln_temp_spr_176to226m    NA NA                    NA   NA                  2               0.889
# 10 long        NA                                   1     Q10             2                          NA                          2 NA                    NA   NA                  1               1    
# 11 short       NA                                   1     Q10             2                          NA                          2 NA                    NA   NA                  1               1    
# 12 long        NA                                   0.472 overlap_form    linear                     NA                         NA linear                NA   NA                  1               1    
# 13 long        NA                                   0.418 overlap_form    logistic                   NA                         NA logistic              NA   NA                  2               0.885
# 14 short       NA                                   0.434 overlap_form    linear                     NA                         NA linear                NA   NA                  2               0.970
# 15 short       NA                                   0.447 overlap_form    logistic                   NA                         NA logistic              NA   NA                  1               1    
# 16 long        NA                                   0.209 overlap_slope   0.2                        NA                         NA NA                     0.2 NA                  4               0.865
# 17 long        NA                                   0.211 overlap_slope   0.4                        NA                         NA NA                     0.4 NA                  3               0.873
# 18 long        NA                                   0.228 overlap_slope   0.6                        NA                         NA NA                     0.6 NA                  2               0.943
# 19 long        NA                                   0.242 overlap_slope   0.8                        NA                         NA NA                     0.8 NA                  1               1    
# 20 short       NA                                   0.228 overlap_slope   0.2                        NA                         NA NA                     0.2 NA                  1               1    
# 21 short       NA                                   0.222 overlap_slope   0.4                        NA                         NA NA                     0.4 NA                  2               0.972
# 22 short       NA                                   0.217 overlap_slope   0.6                        NA                         NA NA                     0.6 NA                  3               0.953
# 23 short       NA                                   0.214 overlap_slope   0.8                        NA                         NA NA                     0.8 NA                  4               0.939
# 24 long        NA                                   0.359 shark_transform log1p_z                    NA                         NA NA                    NA   log1p_z             1               1    
# 25 long        NA                                   0.336 shark_transform z                          NA                         NA NA                    NA   z                   2               0.937
# 26 long        NA                                   0.305 shark_transform z_roll2                    NA                         NA NA                    NA   z_roll2             3               0.849
# 27 short       NA                                   0.360 shark_transform z                          NA                         NA NA                    NA   z                   2               0.963
# 28 short       NA                                   0.374 shark_transform z_roll2                    NA                         NA NA                    NA   z_roll2             1               1    


#===============================================


suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

if (!exists("out_dir")) out_dir <- "copilot/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Requires prebuilt annual tables in env:
# ak_yr, sem_yr, shark_yr, safe_scale, find_col

# ---------------------------
# Fixed SSL model choices
# ---------------------------
col_sst_ssl <- "sst_wgoa_coastwatch_junjulaug"
col_ssl <- "ssl_west_pup_pred"
forage_fixed <- c("stka_herr_matbiom", "mid_il_capelin")

col_salmon_juv <- "x07_dfa_cpue_int_spr_jun_hw"
col_salmon_adult <- "x16_sar"

# Fixed SSL params from your decision
k0 <- 0.0
kT <- 1.5
kF <- 2.0

# Shark-temperature candidates chosen by you
shark_temp_cands <- c(
  "sst_wgoa_coastwatch_junjulaug",
  "bts_wgoa_temp_sum_195to205m",
  "swln_temp_spr_176to226m",
  "ll_wgoa_temp_sum_246to255m",
  "pdo_djf",
  "enso_dj"
)
shark_temp_cands <- shark_temp_cands[shark_temp_cands %in% names(ak_yr)]

if (length(shark_temp_cands) == 0) stop("None of chosen shark temp candidates found in ak_yr.")

# Shark abundance candidates remain open
shark_cands <- find_col(shark_yr, c("shark","biomass","cpue","salmon_shark"))
if (length(shark_cands) == 0) shark_cands <- find_col(sem_yr, c("shark","salmon_shark"))
if (length(shark_cands) == 0) stop("No shark abundance candidates found.")

# Basic checks
stopifnot(
  col_sst_ssl %in% names(ak_yr),
  col_ssl %in% names(ak_yr),
  col_salmon_juv %in% names(sem_yr),
  col_salmon_adult %in% names(sem_yr)
)

forage_ak_use  <- intersect(forage_fixed, names(ak_yr))
forage_sem_use <- intersect(forage_fixed, names(sem_yr))
forage_cols_all <- unique(c(forage_ak_use, forage_sem_use))
if (length(forage_cols_all) == 0) stop("Fixed forage columns not found.")

# ---------------------------
# Helper: build base table once
# ---------------------------
data_base <- tibble(year = sort(unique(c(ak_yr$year, sem_yr$year, shark_yr$year)))) %>%
  left_join(
    ak_yr %>% select(year, SST_ssl = all_of(col_sst_ssl), SSL_raw = all_of(col_ssl), any_of(forage_ak_use), any_of(shark_temp_cands)),
    by = "year"
  ) %>%
  left_join(
    sem_yr %>% select(year, all_of(col_salmon_juv), all_of(col_salmon_adult), any_of(forage_sem_use)),
    by = "year"
  )

# fixed forage + SSL components
forage_scaled <- data_base %>%
  select(any_of(forage_cols_all)) %>%
  mutate(across(everything(), safe_scale))

data_base <- data_base %>%
  mutate(
    F_raw = rowMeans(as.matrix(forage_scaled), na.rm = TRUE),
    T_ssl = safe_scale(SST_ssl),
    SSL_z = safe_scale(SSL_raw),
    F_z   = safe_scale(F_raw)
  )

p_switch_ssl <- plogis(k0 + kT * data_base$T_ssl - kF * data_base$F_z)
I_SSL_fixed <- data_base$SSL_z * p_switch_ssl

# ---------------------------
# Shark model fit function
# ---------------------------
fit_shark_variant <- function(col_shark, shark_temp_col, Q10, overlap_form, overlap_slope, shark_transform) {
  
  d <- data_base %>%
    left_join(shark_yr %>% select(year, Shark_raw = all_of(col_shark)), by = "year")
  
  # shark abundance transform
  shark_vec <- d$Shark_raw
  if (shark_transform == "z") {
    Shark_use <- safe_scale(shark_vec)
  } else if (shark_transform == "log1p_z") {
    minv <- suppressWarnings(min(shark_vec, na.rm = TRUE))
    if (!is.finite(minv)) return(NULL)
    Shark_use <- safe_scale(log1p(shark_vec - minv))
  } else if (shark_transform == "z_roll2") {
    z0 <- safe_scale(shark_vec)
    zlag <- dplyr::lag(z0, 1)
    Shark_use <- rowMeans(cbind(z0, zlag), na.rm = TRUE)
  } else {
    return(NULL)
  }
  
  # shark temperature driver
  if (!shark_temp_col %in% names(d)) return(NULL)
  T_shark_raw <- d[[shark_temp_col]]
  T_shark <- safe_scale(T_shark_raw)
  
  Tref_shark <- mean(T_shark_raw, na.rm = TRUE)
  M_t <- Q10^((T_shark_raw - Tref_shark) / 10)
  
  # overlap forms
  O_t <- if (overlap_form == "constant") {
    rep(1, nrow(d))
  } else if (overlap_form == "linear") {
    pmax(0.1, 1 + overlap_slope * T_shark)
  } else if (overlap_form == "logistic") {
    plogis(overlap_slope * T_shark) * 2  # scaled roughly around 1
  } else {
    return(NULL)
  }
  
  I_Shark <- Shark_use * M_t * O_t
  
  fit_dat <- d %>%
    mutate(
      I_SSL = I_SSL_fixed,
      I_Shark = I_Shark
    ) %>%
    filter(year >= 1998, year <= 2021) %>%
    filter(
      !is.na(.data[[col_salmon_adult]]),
      !is.na(.data[[col_salmon_juv]]),
      !is.na(I_SSL),
      !is.na(I_Shark)
    )
  
  if (nrow(fit_dat) < 12) return(NULL)
  
  # Baseline and test for delta-AIC
  mod_base <- lm(reformulate(c(col_salmon_juv, "I_SSL"), response = col_salmon_adult), data = fit_dat)
  mod_test <- lm(reformulate(c(col_salmon_juv, "I_SSL", "I_Shark"), response = col_salmon_adult), data = fit_dat)
  
  td <- broom::tidy(mod_test)
  sm <- summary(mod_test)
  get_term <- function(term, field) {
    out <- td %>% filter(.data$term == term) %>% pull(!!sym(field))
    if (length(out) == 0) NA_real_ else out[1]
  }
  
  tibble(
    shark_col = col_shark,
    shark_temp_col = shark_temp_col,
    Q10 = Q10,
    overlap_form = overlap_form,
    overlap_slope = overlap_slope,
    shark_transform = shark_transform,
    n = nrow(fit_dat),
    
    aic_base = AIC(mod_base),
    aic_test = AIC(mod_test),
    delta_aic_shark = AIC(mod_test) - AIC(mod_base), # negative = shark helps
    r2_test = sm$r.squared,
    
    # coefficients/p-values for all model terms
    b_cpue  = get_term(col_salmon_juv, "estimate"),
    p_cpue  = get_term(col_salmon_juv, "p.value"),
    b_ssl   = get_term("I_SSL", "estimate"),
    p_ssl   = get_term("I_SSL", "p.value"),
    b_shark = get_term("I_Shark", "estimate"),
    p_shark = get_term("I_Shark", "p.value")
  )
  
}

# ---------------------------
# Sweep grid
# ---------------------------
grid <- tidyr::crossing(
  shark_col = shark_cands,
  shark_temp_col = shark_temp_cands,
  Q10 = c(1.5, 2.0, 2.5, 3.0),
  overlap_form = c("constant", "linear", "logistic"),
  overlap_slope = c(0.2, 0.4, 0.6, 0.8),
  shark_transform = c("z", "log1p_z", "z_roll2")
)

# constant overlap ignores slope; keep it but dedupe after
grid <- grid %>%
  mutate(overlap_slope = if_else(overlap_form == "constant", 0, overlap_slope)) %>%
  distinct()

res_shark <- purrr::pmap_dfr(
  grid,
  function(shark_col, shark_temp_col, Q10, overlap_form, overlap_slope, shark_transform) {
    fit <- fit_shark_variant(
      col_shark = shark_col,
      shark_temp_col = shark_temp_col,
      Q10 = Q10,
      overlap_form = overlap_form,
      overlap_slope = overlap_slope,
      shark_transform = shark_transform
    )
    if (is.null(fit)) tibble() else fit
  }
)

if (nrow(res_shark) == 0) stop("No shark variants fit.")

res_shark <- res_shark %>%
  mutate(
    model_group = case_when(n == 19 ~ "short", n == 24 ~ "long", TRUE ~ NA_character_),
    sig_shark = !is.na(p_shark) & p_shark < 0.05
  )

write_csv(res_shark, file.path(out_dir, "phase5_shark_sweep_all.csv"))

# rank primarily by shark improvement---------
ranked <- res_shark %>%
  mutate(
    sig_cpue  = !is.na(p_cpue)  & p_cpue  < 0.05,
    sig_ssl   = !is.na(p_ssl)   & p_ssl   < 0.05,
    sig_shark = !is.na(p_shark) & p_shark < 0.05,
    sig_all3  = sig_cpue & sig_ssl & sig_shark
  ) %>%
  arrange(
    desc(sig_cpue),      # prioritize models that keep CPUE significant
    desc(sig_ssl),       # then SSL significance
    delta_aic_shark,     # then shark improvement
    p_shark,             # then shark term significance
    desc(abs(b_shark))
  )
write_csv(ranked, file.path(out_dir, "phase5_shark_sweep_ranked.csv"))

# AIC-weighted importance within short/long using test-model AIC
w <- res_shark %>%
  filter(!is.na(model_group), !is.na(aic_test)) %>%
  group_by(model_group) %>%
  mutate(
    delta_aic = aic_test - min(aic_test, na.rm = TRUE),
    rel_like = exp(-0.5 * delta_aic),
    aic_weight = rel_like / sum(rel_like, na.rm = TRUE)
  ) %>%
  ungroup()

imp <- bind_rows(
  w %>% group_by(model_group, shark_col) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_col", level=shark_col),
  w %>% group_by(model_group, shark_temp_col) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_temp_col", level=shark_temp_col),
  w %>% group_by(model_group, Q10) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="Q10", level=as.character(Q10)),
  w %>% group_by(model_group, overlap_form) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="overlap_form", level=overlap_form),
  w %>% group_by(model_group, overlap_slope) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="overlap_slope", level=as.character(overlap_slope)),
  w %>% group_by(model_group, shark_transform) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_transform", level=shark_transform)
) %>%
  group_by(model_group, param) %>%
  mutate(rank = rank(-importance_aicw, ties.method="min"),
         importance_rel_top1 = importance_aicw / max(importance_aicw)) %>%
  ungroup()

write_csv(imp, file.path(out_dir, "phase5_shark_param_importance_short_long.csv"))

message("Saved shark sweep outputs in ", out_dir)

#recompute rankings without "other sharks"---------------
library(tidyverse)

# assumes res_shark already exists from phase5
res_shark_restrict <- res_shark %>%
  filter(
    str_detect(shark_col, "salmon_shark|pacific_sleeper_shark"),
    !str_detect(shark_col, "other_sharks")
  ) %>%
  filter(
    Q10==2
  ) %>%
  mutate(
    sig_cpue  = !is.na(p_cpue)  & p_cpue  < 0.05,
    sig_ssl   = !is.na(p_ssl)   & p_ssl   < 0.05,
    sig_shark = !is.na(p_shark) & p_shark < 0.05,
    sig_all3  = sig_cpue & sig_ssl & sig_shark,
    model_group = case_when(n == 19 ~ "short", n == 24 ~ "long", TRUE ~ NA_character_)
  )

# ranking with CPUE retained as priority
ranked_restrict <- res_shark_restrict %>%
  arrange(
    desc(sig_cpue),
    desc(sig_ssl),
    delta_aic_shark,
    p_shark,
    desc(abs(b_shark))
  )

write_csv(ranked_restrict, file.path(out_dir, "phase5_shark_sweep_ranked_no_other.csv"))
ranked_restrict %>% filter(Q10==2)

# recompute importance without other_sharks
imp_restrict <- res_shark_restrict %>%
  filter(!is.na(model_group), !is.na(aic_test)) %>%
  group_by(model_group) %>%
  mutate(
    delta_aic = aic_test - min(aic_test, na.rm = TRUE),
    rel_like = exp(-0.5 * delta_aic),
    aic_weight = rel_like / sum(rel_like, na.rm = TRUE)
  ) %>%
  ungroup()

imp_params_restrict <- bind_rows(
  imp_restrict %>% group_by(model_group, shark_col) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_col", level=shark_col),
  imp_restrict %>% group_by(model_group, shark_temp_col) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_temp_col", level=shark_temp_col),
  imp_restrict %>% group_by(model_group, Q10) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="Q10", level=as.character(Q10)),
  imp_restrict %>% group_by(model_group, overlap_form) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="overlap_form", level=overlap_form),
  imp_restrict %>% group_by(model_group, overlap_slope) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="overlap_slope", level=as.character(overlap_slope)),
  imp_restrict %>% group_by(model_group, shark_transform) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_transform", level=shark_transform)
) %>%
  group_by(model_group, param) %>%
  mutate(rank = rank(-importance_aicw, ties.method="min"),
         importance_rel_top1 = importance_aicw / max(importance_aicw)) %>%
  ungroup()

write_csv(imp_params_restrict, file.path(out_dir, "phase5_shark_param_importance_no_other.csv"))

print(imp_params_restrict, n = Inf)

imp_params_restrict %>% filter(importance_rel_top1<0.5)
imp_params_restrict %>% filter(importance_rel_top1>0.8)
