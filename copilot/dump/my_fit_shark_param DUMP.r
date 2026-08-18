#Explore best temperature for sharks and other model decisions

sem_master_data<-read.csv("outputs_4/sem_master_data.csv",row.names=NULL);names(sem_master_data)[1]
ak_dat<-read.csv("copilot/outputs_2/data_all_tested_columns_annual.csv");names(ak_dat)[1]
ak_dat<-left_join(sem_master_data %>% select(year,X09_DFA_HakeAge5Plus),ak_dat, join_by("year"))
head(ak_dat)



suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})


if (!exists("out_dir")) out_dir <- "copilot/outputs_4"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


#Salmon variables
col_salmon_juv <- "x07_dfa_cpue_int_spr_jun_hw"
col_salmon_adult <- "x16_sar"


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
shark_cands <- find_col(ak_dat, c("sleeper_shark","salmon_shark"))
if (length(shark_cands) == 0) shark_cands <- find_col(sem_yr, c("sleeper_shark","salmon_shark"))
if (length(shark_cands) == 0) stop("No shark abundance candidates found.")

# Basic checks
stopifnot(
  col_sst_ssl %in% names(ak_yr),
  col_salmon_juv %in% names(sem_yr),
  col_salmon_adult %in% names(sem_yr)
)

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

  
  #NEW SEM FUNCTION================
  single_model_text <- glue('
             # Structural Model
              X07_DFA_cpue_IntSprJunHW ~  X09_DFA_HakeAge5Plus 
              
              X16_SAR ~ X07_DFA_cpue_IntSprJunHW +  
                        I_Shark  
                           
          
          ')
  
  fit_single <- sem(single_model_text, data = dat, std.lv = TRUE,missing = "ML") 
  summary_fxn(fit_single)
  
  #Original linear model---------
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

write_csv(imp, file.path(out_dir, "shark_param_importance_short_long.csv"))

message("Saved shark sweep outputs in ", out_dir)


