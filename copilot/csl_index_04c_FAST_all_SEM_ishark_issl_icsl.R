# ==============================================================================
# FAST 2-VARIABLE ONLY SEM SWEEP (<10 SECONDS RUNTIME)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

out_dir <- "copilot/outputs_7"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load & Assemble Data
# -----------------------------------------------------------------------------

guild.dfasAK <- read.csv("data_Lisa/guild.dfasAK.csv") %>% clean_names()

ssl.dat <- read.csv("data_Lisa/ssl.dat.csv") %>% clean_names() %>%
  select(year, ssl_seak_pup_pred, sst_wgoa_coastwatch_junjulaug) %>%
  left_join(guild.dfasAK %>% select(year, x13_stka_herr_matbiom, x13_egoa_herring, x13_mid_il_capelin), by = "year")

ishark_fxn <- function(N_shark, T_mt, Tref_mt = mean(T_mt, na.rm = TRUE), T_ot, Q10 = 2, overlap_slope = 1, transform_shark = "log1p", scale_final = TRUE) {
  minv <- suppressWarnings(min(N_shark, na.rm = TRUE))
  N_pos <- if (transform_shark == "log1p") log1p(N_shark - minv) else N_shark
  M_t <- Q10^((T_mt - Tref_mt) / 10)
  s <- sd(T_ot, na.rm = TRUE); T_ot_z <- if (is.na(s) || s == 0) rep(0, length(T_ot)) else (T_ot - mean(T_ot, na.rm = TRUE)) / s
  O_t <- plogis(overlap_slope * T_ot_z) * 2
  raw <- N_pos * M_t * O_t
  if (scale_final) { s_r <- sd(raw, na.rm = TRUE); if (is.na(s_r) || s_r == 0) rep(0, length(raw)) else (raw - mean(raw, na.rm = TRUE)) / s_r } else raw
}

shark_raw <- read.csv("data_Lisa/shark_wide.csv") %>% clean_names() %>%
  select(year, goa_pacific_sleeper_shark) %>%
  full_join(guild.dfasAK %>% select(year, x21_sst_wgoa_coastwatch_junjulaug), by = "year") %>%
  mutate(ishark = ishark_fxn(goa_pacific_sleeper_shark, x21_sst_wgoa_coastwatch_junjulaug, T_ot = x21_sst_wgoa_coastwatch_junjulaug))

fit_dat <- guild.dfasAK %>% full_join(shark_raw %>% select(year, ishark), by = "year")

df_top_ssl <- ssl.dat %>%
  mutate(
    ssl = as.vector(scale(ssl_seak_pup_pred)),
    sst = as.vector(scale(sst_wgoa_coastwatch_junjulaug)),
    f_cap = as.vector(scale(x13_mid_il_capelin)),
    f_herr = as.vector(scale(x13_stka_herr_matbiom)),
    I_SSL_simple = -1 * (1.0 * ssl - 0.3 * (ssl * sst) + 2.0 * (ssl * f_herr) + 1.0 * (ssl * f_cap)),
    issl = as.vector(scale(I_SSL_simple))
  ) %>%
  full_join(fit_dat, by = "year")

csl_path <- file.path("copilot/outputs_6", "chinook_numeric_csl_exposure_risk_1998_2024.csv")
csl_raw  <- read.csv(csl_path) %>% clean_names()

safe_z <- function(x) { s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) rep(0, length(x)) else as.numeric((x - mean(x, na.rm = TRUE)) / s) }

csl_idx <- csl_raw %>%
  transmute(
    year,
    x15_icsl_cslchin_ratio = safe_z(log1p(csl_during_chinook) - log1p(chinook_estuary_total))
  ) %>%
  mutate(year = year - 2)

guild <- df_top_ssl %>%
  left_join(csl_idx, by = "year") %>%
  mutate(across(-year, ~ as.vector(scale(.))))

guild_dfas1_24yr <- guild %>%
  filter(year >= 1998, year <= 2021) %>%
  select(where(~ sum(!is.na(.)) >= 24)) %>%
  rename(x15_ishark = ishark, x15_issl = issl) %>%
  select(-any_of(c("x15_shark_catch_go_a_pred_ak", "x10_harbor_seal_cr_2yr_lead")))

# -----------------------------------------------------------------------------
# 2. Extract Candidate Vector (Includes Index Terms explicitly)
# -----------------------------------------------------------------------------

ak_prey_cands <- names(guild_dfas1_24yr)[grepl("(?i)^x21_|^x12_|^x13_|^x14_", names(guild_dfas1_24yr))]
ak_pred_cands <- names(guild_dfas1_24yr)[grepl("(?i)^x10_|^x15_", names(guild_dfas1_24yr))]
index_cands   <- names(guild_dfas1_24yr)[grepl("(?i)ishark|issl|icsl", names(guild_dfas1_24yr))]

all_cands <- unique(c(ak_prey_cands, ak_pred_cands, index_cands))
base_needed <- c("year", "x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "x16_sar")

# -----------------------------------------------------------------------------
# 3. Fast 2-Variable Combination Sweep
# -----------------------------------------------------------------------------

comb2 <- combn(all_cands, 2, simplify = FALSE)

fit_sem_2var <- function(extra_terms) {
  keep_cols <- unique(c(base_needed, extra_terms))
  df_model  <- guild_dfas1_24yr %>% select(all_of(keep_cols)) %>% drop_na()
  if (nrow(df_model) < 24) return(NULL)
  
  rhs <- paste(c("x07_dfa_cpue_int_spr_jun_hw", extra_terms), collapse = " + ")
  sem_text <- paste0("x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus\n", "x16_sar ~ ", rhs, "\n")
  
  fit <- tryCatch(sem(sem_text, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  
  fm <- fitMeasures(fit, c("aic", "bic"))
  r2_sar <- inspect(fit, "r2")[["x16_sar"]]
  pe <- parameterEstimates(fit)
  
  get_coef <- function(term) {
    row <- pe %>% filter(lhs == "x16_sar", op == "~", rhs == term)
    if (nrow(row) == 0) return(c(est = NA_real_, p = NA_real_))
    c(est = row$est[1], p = row$pvalue[1])
  }
  c1 <- get_coef(extra_terms[1])
  c2 <- get_coef(extra_terms[2])
  
  tibble(
    terms = paste(extra_terms, collapse = " + "),
    aic = fm[["aic"]], bic = fm[["bic"]], sar_r2 = as.numeric(r2_sar),
    term1 = extra_terms[1], term2 = extra_terms[2],
    b1 = c1["est"], p1 = c1["p"], b2 = c2["est"], p2 = c2["p"]
  )
}

# Run fast 2-variable loop
cat("Running fast 2-variable SEM sweep...\n")
all_models <- purrr::map_dfr(comb2, fit_sem_2var)

# -----------------------------------------------------------------------------
# 4. Compute Expectations & Variable Importance
# -----------------------------------------------------------------------------

get_expected_sign <- function(var_name) {
  if (is.na(var_name) || var_name == "") return(NA_character_)
  case_when(
    str_detect(var_name, "(?i)ishark|issl|icsl") ~ "Negative",
    str_detect(var_name, "(?i)^x10_|^x15_")       ~ "Negative",
    str_detect(var_name, "(?i)^x21_")              ~ "Negative",
    str_detect(var_name, "(?i)^x12_|^x13_")        ~ "Positive",
    str_detect(var_name, "(?i)^x14_")              ~ "Negative",
    TRUE                                           ~ "Unknown"
  )
}

check_sign <- function(b_est, exp_sign) {
  if (is.na(b_est) || is.na(exp_sign) || exp_sign == "Unknown") return(NA)
  if (exp_sign == "Positive") return(b_est > 0)
  if (exp_sign == "Negative") return(b_est < 0)
  return(NA)
}

all_models <- all_models %>%
  mutate(
    delta_aic_global = aic - min(aic, na.rm = TRUE),
    aic_weight = exp(-0.5 * delta_aic_global) / sum(exp(-0.5 * delta_aic_global), na.rm = TRUE),
    
    exp_sign_1 = sapply(term1, get_expected_sign),
    exp_sign_2 = sapply(term2, get_expected_sign),
    
    pass_sign_1 = mapply(check_sign, b1, exp_sign_1),
    pass_sign_2 = mapply(check_sign, b2, exp_sign_2),
    
    n_passed_terms = replace_na(pass_sign_1, FALSE) + replace_na(pass_sign_2, FALSE),
    prop_model_signs_met = n_passed_terms / 2
  ) %>%
  arrange(aic)

terms_eval <- bind_rows(
  all_models %>% select(aic_weight, variable = term1, estimate = b1, exp_sign = exp_sign_1, sign_met = pass_sign_1),
  all_models %>% select(aic_weight, variable = term2, estimate = b2, exp_sign = exp_sign_2, sign_met = pass_sign_2)
) %>%
  filter(!is.na(variable)) %>%
  mutate(
    var_group = case_when(
      str_detect(variable, "(?i)ishark|issl|icsl") ~ "Index",
      str_detect(variable, "(?i)^x10_|^x15_")       ~ "Predator",
      str_detect(variable, "(?i)^x21_")              ~ "Climate",
      str_detect(variable, "(?i)^x12_|^x13_")        ~ "Prey",
      str_detect(variable, "(?i)^x14_")              ~ "Competitor",
      TRUE                                           ~ "Other"
    )
  )

category_summary <- terms_eval %>%
  group_by(var_group, exp_sign) %>%
  summarise(
    evaluations           = n(),
    raw_prop_met          = round(mean(sign_met, na.rm = TRUE), 3),
    aic_weighted_prop_met = round(sum(aic_weight[sign_met == TRUE], na.rm = TRUE) / sum(aic_weight, na.rm = TRUE), 3),
    .groups               = "drop"
  )%>%
  arrange(desc(aic_weighted_prop_met))

var_importance <- terms_eval %>%
  group_by(variable, var_group, exp_sign) %>%
  summarise(
    importance            = round(sum(aic_weight, na.rm = TRUE), 4),
    prop_sign_met_aic_wtd = round(sum(aic_weight[sign_met == TRUE], na.rm = TRUE) / sum(aic_weight, na.rm = TRUE), 3),
    .groups               = "drop"
  ) %>%
  arrange(desc(importance))

# Write fast presentation deliverables
write.csv(all_models, file.path(out_dir, "fast_2var_all_models.csv"), row.names = FALSE)
write.csv(category_summary, file.path(out_dir, "fast_category_expectation_summary.csv"), row.names = FALSE)
write.csv(var_importance, file.path(out_dir, "fast_variable_importance.csv"), row.names = FALSE)

cat("\nFINISHED! Ready for presentation.\n\n")
cat("=== TOP 10 MODELS BY AIC ===\n")
print(head(all_models %>% select(terms, aic, sar_r2, prop_model_signs_met), 10))

cat("\n=== CATEGORY EXPECTATION SUCCESS RATES ===\n")
print(category_summary)

cat("\n=== TOP 15 VARIABLE IMPORTANCE WITH SIGN SUCCESS RATES ===\n")
print(head(var_importance, 15))
