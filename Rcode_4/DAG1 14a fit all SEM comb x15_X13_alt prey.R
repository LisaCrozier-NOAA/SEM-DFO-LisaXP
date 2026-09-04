
#NEED TO FIX PREY YEARS TO 2YR LEAD, OR UNLEAD THE PREDATORS

# ==============================================================================
# Script: Systematic Alternate Prey Interaction SEM Sweep
# Tests: Predator + Prey1 + Prey2 + (Predator * Prey1) + (Predator * Prey2)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

out_dir <- "copilot/outputs_8"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load Data & Exclude Bespoke Index Terms
# -----------------------------------------------------------------------------

guild.dfasAK <- read.csv("data_Lisa/guild.dfasAK.csv",row.names=1) %>% clean_names()

ssl.dat <- read.csv("data_Lisa/ssl.dat.csv") %>% clean_names() %>%
  select(year, ssl_seak_pup_pred) %>%
  rename(x15_ssl_seak_pup_pred=ssl_seak_pup_pred)

# Assemble raw guild dataset without precalculated indices
guild <- ssl.dat %>%
  full_join(guild.dfasAK, by = "year") %>%
  mutate(across(-year, ~ as.vector(scale(.))))

guild_dfas1_24yr <- guild %>%
  filter(year >= 1998, year <= 2021) %>%
  select(where(~ sum(!is.na(.)) >= 24)) %>%
  # Strip duplicate column suffixes (.x / .y) if present
  rename_with(~ str_remove(., "\\.[xy]$")) %>%
  # Remove precalculated index terms
  select(-any_of(c("ishark", "issl", "icsl", "x15_ishark", "x15_issl", "x15_icsl_cslchin_ratio", "x15_icsl_risk_eulachon75_shad25")))

# -----------------------------------------------------------------------------
# 2. Extract Candidate Variables
# -----------------------------------------------------------------------------

# Predators (x10_ or x15_)
ak_pred_cands <- names(guild_dfas1_24yr)[grepl("(?i)^x10_|^x15_", names(guild_dfas1_24yr))]

# Prey & Climate/Temperature (x12_, x13_, x14_, x21_)
ak_prey_temp_cands <- names(guild_dfas1_24yr)[grepl("(?i)^x21_|^x12_|^x13_|^x14_", names(guild_dfas1_24yr))]

base_needed <- c("year", "x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "x16_sar")

# -----------------------------------------------------------------------------
# 3. Build Interaction Pairs & Fit SEM Sweep
# -----------------------------------------------------------------------------

# Combine all 2-variable combinations of Prey/Temperature
prey_pairs <- combn(ak_prey_temp_cands, 2, simplify = FALSE)

fit_interaction_sem <- function(pred_var, prey_pair) {
  prey1 <- prey_pair[1]
  prey2 <- prey_pair[2]
  
  # Construct interaction product names
  int1_name <- paste0("int_", pred_var, "_x_", prey1)
  int2_name <- paste0("int_", pred_var, "_x_", prey2)
  
  # Build modeling dataframe with explicit product terms
  df_model <- guild_dfas1_24yr %>%
    select(all_of(c(base_needed, pred_var, prey1, prey2))) %>%
    mutate(
      int1 = .data[[pred_var]] * .data[[prey1]],
      int2 = .data[[pred_var]] * .data[[prey2]]
    ) %>%
    drop_na()
  
  if (nrow(df_model) < 24) return(NULL)
  
  # SEM Equation text including main effects + interaction products
  sem_text <- paste0(
    "x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus\n",
    "x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + ", pred_var, " + ", prey1, " + ", prey2, " + int1 + int2\n"
  )
  
  fit <- tryCatch(
    sem(sem_text, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  
  fm <- fitMeasures(fit, c("aic", "bic"))
  r2_sar <- inspect(fit, "r2")[["x16_sar"]]
  pe <- parameterEstimates(fit)
  
  get_coef <- function(term) {
    row <- pe %>% filter(lhs == "x16_sar", op == "~", rhs == term)
    if (nrow(row) == 0) return(c(est = NA_real_, p = NA_real_))
    c(est = row$est[1], p = row$pvalue[1])
  }
  
  b_pred <- get_coef(pred_var)
  b_p1   <- get_coef(prey1)
  b_p2   <- get_coef(prey2)
  b_i1   <- get_coef("int1")
  b_i2   <- get_coef("int2")
  
  tibble(
    predator = pred_var,
    prey1 = prey1,
    prey2 = prey2,
    model_formula = paste0(pred_var, " * (", prey1, " + ", prey2, ")"),
    aic = fm[["aic"]],
    bic = fm[["bic"]],
    sar_r2 = as.numeric(r2_sar),
    b_pred = b_pred["est"], p_pred = b_pred["p"],
    b_prey1 = b_p1["est"],  p_prey1 = b_p1["p"],
    b_prey2 = b_p2["est"],  p_prey2 = b_p2["p"],
    b_int1 = b_i1["est"],   p_int1 = b_i1["p"],
    b_int2 = b_i2["est"],   p_int2 = b_i2["p"]
  )
}

# Run loop over all Predator x (Prey1, Prey2) combinations
cat("Running Alternate Prey Interaction SEM Sweep...\n")

all_interaction_models <- expand_grid(
  pred = ak_pred_cands,
  pair_idx = seq_along(prey_pairs)
) %>%
  mutate(
    res = purrr::map2(pred, pair_idx, ~ fit_interaction_sem(.x, prey_pairs[[.y]]))
  ) %>%
  pull(res) %>%
  bind_rows()

# Compute AIC weights
all_interaction_models <- all_interaction_models %>%
  mutate(
    delta_aic = aic - min(aic, na.rm = TRUE),
    aic_weight = exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic), na.rm = TRUE)
  ) %>%
  arrange(aic)

# Save deliverables
write.csv(all_interaction_models, file.path(out_dir, "alternate_prey_interaction_models_full.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Console Summary
# -----------------------------------------------------------------------------

cat("\n====================================================================\n")
cat(" TOP 10 ALTERNATE PREY INTERACTION MODELS BY AIC                    \n")
cat("====================================================================\n")

print(
  all_interaction_models %>%
    select(model_formula, aic, delta_aic, sar_r2, b_int1, p_int1, b_int2, p_int2) %>%
    slice_head(n = 10)
)

print(
  all_interaction_models %>%
    filter(b_int1<0, b_int2<0) %>%
    select(model_formula, delta_aic, b_int1, b_int2) %>%
    slice_head(n = 10)
)


#RESULTS-----------
# model_formula                                                                                  delta_aic  b_int1  b_int2
# 1 x15_halibut_biomass_age8plus_2yr_lead_pred_ak * (x13_wgoa_cap_pcod + x13_stka_herr_matbiom)         0     2.32   -1.19  
# 2 x15_ssl_seak_pup_pred * (x13_wgoa_cap_pcod + x13_stka_herr_matbiom)                                 2.59 -0.889   1.59  
# 3 x15_ssl_seak_pup_pred * (x21_sst_wgoa_coastwatch_junjulaug + x13_stka_herr_matbiom)                 3.75  0.177   0.822 
# 4 x15_dfa_sleeper_shark_bsai_pred_ak * (x13_dfa_wgoa_dfa_seabirds + x13_stka_herr_matbiom)            4.34 -0.0805 -1.16  
# 5 x15_ssl_seak_pup_pred * (x13_dfa_wgoa_dfa_seabirds + x13_stka_herr_matbiom)                         4.92  0.0392  1.26  
# 6 x15_shark_catch_go_a_pred_ak * (x13_hexagram_eai + x21_sst_wgoa_coastwatch_junjulaug)               5.18 -0.359   0.343 
# 7 x15_arrowtooth_flounder_biomass_pred_ak * (x13_dfa_wgoa_dfa_seabirds + x13_egoa_herring)            5.90 -0.367  -2.04  
# 8 x15_spiny_dogfish_go_a_pred_ak * (x13_hexagram_eai + x21_sst_wgoa_coastwatch_junjulaug)             6.23 -0.336   0.353 
# 9 x10_dfa_ssl_est_wholerange_2yr_lead * (x21_sst_wgoa_coastwatch_junjulaug + x13_stka_herr_matb…      6.59  0.443   0.0328
#10 x15_halibut_biomass_age8plus_2yr_lead_pred_ak * (x13_wgoa_cap_pcod + x13_egoa_herring)              6.90  2.17   -2.06  

#JUST NEGATIVE COEF
# model_formula                                                                                  delta_aic  b_int1  b_int2
# 1 x15_dfa_sleeper_shark_bsai_pred_ak * (x13_dfa_wgoa_dfa_seabirds + x13_stka_herr_matbiom)            4.34 -0.0805 -1.16  
# 2 x15_arrowtooth_flounder_biomass_pred_ak * (x13_dfa_wgoa_dfa_seabirds + x13_egoa_herring)            5.90 -0.367  -2.04  
# 3 x15_arrowtooth_flounder_biomass_pred_ak * (x21_sst_wgoa_coastwatch_junjulaug + x13_stka_herr_…      7.58 -0.367  -0.0810
# 4 x15_halibut_biomass_age8plus_2yr_lead_pred_ak * (x21_sst_wgoa_coastwatch_junjulaug + x13_stka…      7.59 -0.427  -0.155 
# 5 x15_dfa_sleeper_shark_bsai_pred_ak * (x21_sst_wgoa_coastwatch_junjulaug + x13_stka_herr_matbi…      8.20 -0.131  -0.415 
  

#NEED TO FIX PREY YEARS TO 2YR LEAD, OR UNLEAD THE PREDATORS                                                                                                                                                                                                                                                                                  