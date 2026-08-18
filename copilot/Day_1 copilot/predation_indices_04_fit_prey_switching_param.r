#Results: selected values for ssl prey switching function
# starting parameters
# k0 <- 0.0
# kW <- 1.0   # warm effect on switching (+)
# kF <- 1.0   # forage effect reducing switching (+ in formula as -kF*F_t)
# Q10 <- 2.0
# p_switch <- plogis(k0 + kW * design$W_use - kF * design$F_z)

# SSL switching: increases with T_use and decreases with forage F_z
p_switch <- plogis(k0 + kT * dat$sst - kF * dat$forage)

C_t   <- rep(1, nrow(data_selected))
I_SSL <- data_selected$SSL_z * p_switch * C_t


# For mainline model going forward:
#   
#   k0 = 0
#   kT = 1.5
#   kF = 2.0
#   Q10 = 2.0

#results: not sensitive at all to Q10. Let's just pick the standard in the literature.  
#k0 was -1 or -0.5 in the long model, and 0 or 0.5 in the short model. What does the negative vs positive mean? Should we just go with 0? I don't want to have different values for this in the short vs long models.  
#Let's go with 1.5 or 2 for kT -- is there a literature standard? 
#And 2 for kF

# param_importance %>% filter(importance_rel_top1>0.7)
# 
# model_group param     level             importance_aicw importance_rel_top1  rank
# <chr>       <chr>     <chr>                       <dbl>               <dbl> <int>
#   1 long        k0        -1                          0.486               1         1
# 2 long        k0        -0.5                        0.345               0.709     2
# 3 short       k0        0                           0.353               1         1
# 4 short       k0        0.5                         0.330               0.937     2
# 5 long        kT        2                           0.712               1         1
# 6 short       kT        1                           0.287               0.701     2
# 7 short       kT        1.5                         0.409               1         1
# 8 long        kF        1.5                         0.364               0.740     2
# 9 long        kF        2                           0.493               1         1
# 10 short       kF        2                           0.656               1         1
# 11 long        Q10       1.5                         0.341               1         1
# 12 long        Q10       2                           0.333               0.976     2
# 13 long        Q10       2.5                         0.326               0.957     3
# 14 short       Q10       1.5                         0.339               1         1
# 15 short       Q10       2                           0.333               0.983     2
# 16 short       Q10       2.5                         0.328               0.969     3
# 17 long        shark_col other_sharks_goa            0.964               1         1
# 18 short       shark_col other_sharks_bsai           0.574               1         1


# What the signs mean for k0
# In your switching function:
#   
#   p_switch = plogis(k0 + kT*T_use - kF*F_z)
# 
# k0 is baseline log-odds of switching when T_use = 0 and F_z = 0 (average conditions).
# k0 < 0 ⇒ lower baseline switching probability.
# k0 > 0 ⇒ higher baseline switching probability.
# So long-model preference for negative k0 means “low baseline switching unless conditions push it up,” while short-model preference near 0/0.5 means a moderate baseline.
# 
# Given you want one shared value across short/long, k0 = 0 is a very reasonable compromise (neutral baseline).
# 
# kT choice
# No universal “standard” scalar for this exact logistic switching form (because it depends on how covariates are scaled), so your empirical range matters more than external default.
# Given your support at 1.5–2, I’d pick:
#   
#   kT = 1.5 as a conservative primary value
# and run sensitivity at 2.0 (already basically your top support)
# kF choice
# Agree with your read: kF = 2 is strongly supported.
# 
# Q10
# Also agree. Since insensitive, use standard literature default: Q10 = 2.0.


#setup-------

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

if (!exists("out_dir")) out_dir <- "copilot/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Required objects/functions expected from prior scripts:
# ak_raw, sem, shark, clim_tr(optional), guess_year_col(), find_col(), safe_scale()

# ---------------------------
# 0) Fixed choices from Phase 2
# ---------------------------
col_sst <- "sst_wgoa_coastwatch_junjulaug"
col_ssl <- "ssl_west_pup_pred"
use_climate <- FALSE  # fixed to T_sst
forage_fixed <- c("stka_herr_matbiom", "mid_il_capelin")

col_salmon_juv   <- "x07_dfa_cpue_int_spr_jun_hw"
col_salmon_adult <- "x16_sar"

# ---------------------------
# 1) Build annual tables
# ---------------------------
year_col_ak <- guess_year_col(ak_raw)
year_col_sem <- guess_year_col(sem)
year_col_shark <- guess_year_col(shark)

ak_yr <- ak_raw %>%
  rename(year = all_of(year_col_ak)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

sem_yr <- sem %>%
  rename(year = all_of(year_col_sem)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

shark_yr <- shark %>%
  rename(year = all_of(year_col_shark)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# checks
stopifnot(col_sst %in% names(ak_yr))
stopifnot(col_ssl %in% names(ak_yr))
stopifnot(col_salmon_juv %in% names(sem_yr))
stopifnot(col_salmon_adult %in% names(sem_yr))

# shark candidates stay open
shark_cands <- find_col(shark_yr, c("shark","biomass","cpue","salmon_shark"))
if (length(shark_cands) == 0) shark_cands <- find_col(sem_yr, c("shark","salmon_shark"))
if (length(shark_cands) == 0) stop("No shark candidates found.")

# forage fixed, but can live in AK and/or SEM annual tables
forage_ak_use  <- intersect(forage_fixed, names(ak_yr))
forage_sem_use <- intersect(forage_fixed, names(sem_yr))
forage_cols_all <- unique(c(forage_ak_use, forage_sem_use))
if (length(forage_cols_all) == 0) stop("Fixed forage columns not found in annual tables.")

# ---------------------------
# 2) Fit function with locked SST/SSL/forage, variable shark + params
# ---------------------------
fit_shape_combo <- function(col_shark, k0, kT, kF, Q10) {
  
  data_selected <- tibble(year = sort(unique(c(ak_yr$year, sem_yr$year, shark_yr$year)))) %>%
    left_join(
      ak_yr %>%
        select(year, SST_raw = all_of(col_sst), SSL_raw = all_of(col_ssl), any_of(forage_ak_use)),
      by = "year"
    ) %>%
    left_join(
      sem_yr %>%
        select(year, all_of(col_salmon_juv), all_of(col_salmon_adult), any_of(forage_sem_use)),
      by = "year"
    ) %>%
    left_join(
      shark_yr %>% select(year, Shark_raw = all_of(col_shark)),
      by = "year"
    )
  
  # forage composite from fixed set
  forage_scaled <- data_selected %>%
    select(any_of(forage_cols_all)) %>%
    mutate(across(everything(), safe_scale))
  
  data_selected <- data_selected %>%
    mutate(F_raw = rowMeans(as.matrix(forage_scaled), na.rm = TRUE))
  
  T_ref <- mean(data_selected$SST_raw, na.rm = TRUE)
  
  data_selected <- data_selected %>%
    mutate(
      T_sst   = safe_scale(SST_raw),
      T_use   = T_sst,
      SSL_z   = safe_scale(SSL_raw),
      Shark_z = safe_scale(Shark_raw),
      F_z     = safe_scale(F_raw)
    )
  
  # SSL switching shape
  p_switch <- plogis(k0 + kT * data_selected$T_use - kF * data_selected$F_z)
  I_SSL <- data_selected$SSL_z * p_switch
  
  # Shark metabolic shape
  O_t <- pmax(0.2, 1 + 0.3 * data_selected$T_use)
  M_t <- Q10^((data_selected$SST_raw - T_ref)/10)
  I_Shark <- data_selected$Shark_z * M_t * O_t
  
  fit_dat <- data_selected %>%
    mutate(
      p_switch_ssl = p_switch,
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
  
  if (nrow(fit_dat) < 12) return(NULL)
  
  mod <- lm(
    reformulate(c(col_salmon_juv, "I_SSL", "I_Shark"), response = col_salmon_adult),
    data = fit_dat
  )
  
  td <- broom::tidy(mod)
  sm <- summary(mod)
  
  get_term <- function(term, field) {
    out <- td %>% filter(.data$term == term) %>% pull(!!sym(field))
    if (length(out)==0) NA_real_ else out[1]
  }
  
  tibble(
    sst_col = col_sst,
    forcing = "T_sst",
    ssl_col = col_ssl,
    forage_label = paste(forage_fixed, collapse = " + "),
    shark_col = col_shark,
    k0 = k0, kT = kT, kF = kF, Q10 = Q10,
    n = nrow(fit_dat),
    aic = AIC(mod),
    r2 = sm$r.squared,
    b_cpue = get_term(col_salmon_juv, "estimate"),
    b_ssl = get_term("I_SSL", "estimate"),
    b_shark = get_term("I_Shark", "estimate"),
    p_cpue = get_term(col_salmon_juv, "p.value"),
    p_ssl = get_term("I_SSL", "p.value"),
    p_shark = get_term("I_Shark", "p.value")
  )
}

# ---------------------------
# 3) Parameter grid
# ---------------------------
param_grid <- tidyr::crossing(
  shark_col = shark_cands,
  k0 = c(-1.5, -1.0, -0.5, 0.0, 0.5, 1.0),
  kT = c(0.5, 1.0, 1.5, 2.0),
  kF = c(0.5, 1.0, 1.5, 2.0),
  Q10 = c(1.5, 2.0, 2.5)
)

res <- purrr::pmap_dfr(
  param_grid,
  function(shark_col, k0, kT, kF, Q10) {
    fit <- fit_shape_combo(col_shark  = shark_col, k0 = k0, kT = kT, kF = kF, Q10 = Q10)
    if (is.null(fit)) tibble() else fit
  }
)

if (nrow(res) == 0) stop("No models fit in function-shape sweep.")

# ---------------------------
# 4) Rank + short/long summaries
# ---------------------------
res <- res %>%
  mutate(
    model_group = case_when(
      n == 19 ~ "short",
      n == 24 ~ "long",
      TRUE ~ NA_character_
    ),
    sig_cpue  = !is.na(p_cpue)  & p_cpue  < 0.05,
    sig_ssl   = !is.na(p_ssl)   & p_ssl   < 0.05,
    sig_shark = !is.na(p_shark) & p_shark < 0.05,
    sig_all3  = sig_cpue & sig_ssl & sig_shark
  )

res_ranked <- res %>% arrange(aic, desc(r2))
write_csv(res, file.path(out_dir, "phase3_function_shape_all.csv"))
write_csv(res_ranked, file.path(out_dir, "phase3_function_shape_ranked.csv"))

# within-group AIC weights
res_w <- res %>%
  filter(!is.na(model_group), !is.na(aic)) %>%
  group_by(model_group) %>%
  mutate(
    delta_aic = aic - min(aic, na.rm = TRUE),
    rel_like = exp(-0.5 * delta_aic),
    aic_weight = rel_like / sum(rel_like, na.rm = TRUE)
  ) %>%
  ungroup()

write_csv(res_w, file.path(out_dir, "phase3_function_shape_with_aicw.csv"))

# parameter importance by group
param_importance <- bind_rows(
  res_w %>% group_by(model_group, k0)  %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="k0", level=as.character(k0)) %>% select(model_group, param, level, importance_aicw),
  res_w %>% group_by(model_group, kT)  %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="kT", level=as.character(kT)) %>% select(model_group, param, level, importance_aicw),
  res_w %>% group_by(model_group, kF)  %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="kF", level=as.character(kF)) %>% select(model_group, param, level, importance_aicw),
  res_w %>% group_by(model_group, Q10) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="Q10", level=as.character(Q10)) %>% select(model_group, param, level, importance_aicw),
  res_w %>% group_by(model_group, shark_col) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_col", level=shark_col) %>% select(model_group, param, level, importance_aicw)
) %>%
  group_by(model_group, param) %>%
  mutate(
    importance_rel_top1 = importance_aicw / max(importance_aicw, na.rm = TRUE),
    rank = rank(-importance_aicw, ties.method = "min")
  ) %>%
  ungroup()

write_csv(param_importance, file.path(out_dir, "phase3_parameter_importance_short_long.csv"))

param_importance %>% filter(importance_rel_top1>0.7)

# ---------------------------
# 5) Graph functional response for best significant model per group
# ---------------------------
best_sig <- res_w %>%
  filter(sig_all3) %>%
  group_by(model_group) %>%
  arrange(aic, desc(r2), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

write_csv(best_sig, file.path(out_dir, "phase3_best_significant_by_group.csv"))

rebuild_best_data <- function(row_one) {
  col_shark <- row_one$shark_col[[1]]
  k0 <- row_one$k0[[1]]
  kT <- row_one$kT[[1]]
  kF <- row_one$kF[[1]]
  Q10 <- row_one$Q10[[1]]
  
  data_selected <- tibble(year = sort(unique(c(ak_yr$year, sem_yr$year, shark_yr$year)))) %>%
    left_join(ak_yr %>% select(year, SST_raw = all_of(col_sst), SSL_raw = all_of(col_ssl), any_of(forage_ak_use)), by = "year") %>%
    left_join(sem_yr %>% select(year, all_of(col_salmon_juv), all_of(col_salmon_adult), any_of(forage_sem_use)), by = "year") %>%
    left_join(shark_yr %>% select(year, Shark_raw = all_of(col_shark)), by = "year")
  
  forage_scaled <- data_selected %>%
    select(any_of(forage_cols_all)) %>%
    mutate(across(everything(), safe_scale))
  data_selected <- data_selected %>% mutate(F_raw = rowMeans(as.matrix(forage_scaled), na.rm = TRUE))
  
  T_ref <- mean(data_selected$SST_raw, na.rm = TRUE)
  
  data_selected <- data_selected %>%
    mutate(
      T_sst   = safe_scale(SST_raw),
      T_use   = T_sst,
      SSL_z   = safe_scale(SSL_raw),
      Shark_z = safe_scale(Shark_raw),
      F_z     = safe_scale(F_raw)
    )
  
  p_switch <- plogis(k0 + kT * data_selected$T_use - kF * data_selected$F_z)
  I_SSL <- data_selected$SSL_z * p_switch
  O_t <- pmax(0.2, 1 + 0.3 * data_selected$T_use)
  M_t <- Q10^((data_selected$SST_raw - T_ref)/10)
  I_Shark <- data_selected$Shark_z * M_t * O_t
  
  data_selected %>%
    mutate(
      p_switch_ssl = p_switch,
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
}

if (nrow(best_sig) > 0) {
  for (i in seq_len(nrow(best_sig))) {
    r <- best_sig[i, ]
    grp <- r$model_group[[1]]
    dat <- rebuild_best_data(r)
    
    p_sw <- ggplot(dat, aes(x = F_z, y = p_switch_ssl, color = T_use)) +
      geom_point(size = 2, alpha = 0.85) +
      scale_color_viridis_c() +
      theme_bw() +
      labs(
        title = paste0("Best ", grp, " model: SSL switching response"),
        subtitle = paste0("k0=", r$k0[[1]], ", kT=", r$kT[[1]], ", kF=", r$kF[[1]], ", Q10=", r$Q10[[1]], ", shark=", r$shark_col[[1]]),
        x = "F_z (forage, scaled)",
        y = "p_switch_ssl",
        color = "T_use"
      )
    
    p_ssl <- ggplot(dat, aes(x = I_SSL, y = .data[[col_salmon_adult]])) +
      geom_point(alpha = 0.8) +
      theme_bw() +
      labs(
        title = paste0("Best ", grp, " model: SAR vs I_SSL"),
        x = "I_SSL", y = col_salmon_adult
      )
    
    p_shk <- ggplot(dat, aes(x = I_Shark, y = .data[[col_salmon_adult]])) +
      geom_point(alpha = 0.8, color = "firebrick") +
      theme_bw() +
      labs(
        title = paste0("Best ", grp, " model: SAR vs I_Shark"),
        x = "I_Shark", y = col_salmon_adult
      )
    
    ggsave(file.path(out_dir, paste0("phase3_best_", grp, "_switch_response.png")), p_sw, width = 8, height = 5, dpi = 150)
    ggsave(file.path(out_dir, paste0("phase3_best_", grp, "_sar_vs_issl.png")), p_ssl, width = 7, height = 5, dpi = 150)
    ggsave(file.path(out_dir, paste0("phase3_best_", grp, "_sar_vs_ishark.png")), p_shk, width = 7, height = 5, dpi = 150)
  }
}

print(p_sw)
print(p_ssl)
print(p_shk)

message("Done.\n",
        "- phase3_function_shape_all.csv\n",
        "- phase3_function_shape_ranked.csv\n",
        "- phase3_function_shape_with_aicw.csv\n",
        "- phase3_parameter_importance_short_long.csv\n",
        "- phase3_best_significant_by_group.csv\n")
