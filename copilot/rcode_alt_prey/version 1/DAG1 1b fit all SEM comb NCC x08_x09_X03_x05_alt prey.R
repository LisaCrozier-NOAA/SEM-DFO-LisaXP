# ==============================================================================
# Script: Systematic Alternate Prey Interaction SEM Sweep - NCC Domain
# Predicts Smolt Metric (x07_dfa_cpue) using NCC Predators, Prey, & Climate
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

out_dir <- "copilot/outputs_8"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load Data & Clean
# -----------------------------------------------------------------------------

guild.dfasAK <- read.csv("data_Lisa/guild.dfasAK.csv",row.names=1) %>% clean_names()

ssl.dat <- read.csv("data_Lisa/ssl.dat.csv") %>% clean_names() %>%
  select(year, ssl_seak_pup_pred) %>%
  rename(x15_ssl_seak_pup_pred=ssl_seak_pup_pred)

guild <- ssl.dat %>%
  full_join(guild.dfasAK, by = "year") %>%
  mutate(across(-year, ~ as.vector(scale(.))))

guild_dfas1_24yr <- guild %>%
  filter(year >= 1998, year <= 2021) %>%
  select(where(~ sum(!is.na(.)) >= 24)) %>%
  rename_with(~ str_remove(., "\\.[xy]$"))

# -----------------------------------------------------------------------------
# 2. Extract Candidate Vectors for NCC
# -----------------------------------------------------------------------------

# NCC Predators (x08_, x09_, x10_, x11_)
ncc_pred_cands <- names(guild_dfas1_24yr)[grepl("(?i)^x08_|^x09_|^x10_|^x11_", names(guild_dfas1_24yr))]

# NCC Prey & Temperature (x21_, x03_, x05_)
ncc_prey_temp_cands <- names(guild_dfas1_24yr)[grepl("(?i)^x21_|^x03_|^x05_", names(guild_dfas1_24yr))]

base_needed <- c("year", "x07_dfa_cpue_int_spr_jun_hw", "x16_sar")

# -----------------------------------------------------------------------------
# 3. Fit NCC Interaction SEM Sweep
# -----------------------------------------------------------------------------

prey_pairs_ncc <- combn(ncc_prey_temp_cands, 2, simplify = FALSE)

fit_ncc_interaction_sem <- function(pred_var, prey_pair) {
  prey1 <- prey_pair[1]
  prey2 <- prey_pair[2]
  
  df_model <- guild_dfas1_24yr %>%
    select(all_of(c(base_needed, pred_var, prey1, prey2))) %>%
    mutate(
      int1 = .data[[pred_var]] * .data[[prey1]],
      int2 = .data[[pred_var]] * .data[[prey2]]
    ) %>%
    drop_na()
  
  if (nrow(df_model) < 24) return(NULL)
  
  # SEM Equation: Interaction predicts smolt metric (x07), which then predicts adult return (x16)
  sem_text <- paste0(
    "x07_dfa_cpue_int_spr_jun_hw ~ ", pred_var, " + ", prey1, " + ", prey2, " + int1 + int2\n",
    "x16_sar ~ x07_dfa_cpue_int_spr_jun_hw\n"
  )
  
  fit <- tryCatch(
    sem(sem_text, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  
  fm <- fitMeasures(fit, c("aic", "bic"))
  r2_smolt <- inspect(fit, "r2")[["x07_dfa_cpue_int_spr_jun_hw"]]
  r2_sar   <- inspect(fit, "r2")[["x16_sar"]]
  pe <- parameterEstimates(fit)
  
  get_coef <- function(term) {
    row <- pe %>% filter(lhs == "x07_dfa_cpue_int_spr_jun_hw", op == "~", rhs == term)
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
    model_formula = paste0("x07 ~ ", pred_var, " * (", prey1, " + ", prey2, ")"),
    aic = fm[["aic"]],
    bic = fm[["bic"]],
    smolt_r2 = as.numeric(r2_smolt),
    sar_r2   = as.numeric(r2_sar),
    b_pred = b_pred["est"], p_pred = b_pred["p"],
    b_prey1 = b_p1["est"],  p_prey1 = b_p1["p"],
    b_prey2 = b_p2["est"],  p_prey2 = b_p2["p"],
    b_int1 = b_i1["est"],   p_int1 = b_i1["p"],
    b_int2 = b_i2["est"],   p_int2 = b_i2["p"]
  )
}

cat("Running NCC Alternate Prey Interaction SEM Sweep...\n")

all_ncc_models <- expand_grid(
  pred = ncc_pred_cands,
  pair_idx = seq_along(prey_pairs_ncc)
) %>%
  mutate(
    res = purrr::map2(pred, pair_idx, ~ fit_ncc_interaction_sem(.x, prey_pairs_ncc[[.y]]))
  ) %>%
  pull(res) %>%
  bind_rows()

# Compute AIC weights
all_ncc_models <- all_ncc_models %>%
  mutate(
    delta_aic = aic - min(aic, na.rm = TRUE),
    aic_weight = exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic), na.rm = TRUE)
  ) %>%
  arrange(aic)

# Export deliverables
write.csv(all_ncc_models, file.path(out_dir, "ncc_alternate_prey_interaction_models.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Console Summary
# -----------------------------------------------------------------------------

cat("\n====================================================================\n")
cat(" TOP 10 NCC ALTERNATE PREY INTERACTION MODELS BY AIC               \n")
cat("====================================================================\n")

print(
  all_ncc_models %>%
    select(model_formula, aic, delta_aic, smolt_r2, b_int1, p_int1, b_int2, p_int2) %>%
    slice_head(n = 10)
)

print(
  all_ncc_models %>%
    filter(b_int1<0, b_int2<0) %>%
    select(model_formula, delta_aic,b_pred, b_int1, b_int2) %>%
    slice_head(n = 10)
)

#comment on expected signs:
#1. Main Effects

    #Predator ($\beta_{pred}$): EXPECTED NEGATIVE ($< 0$)
    #Reasoning: Higher abundance of NCC predators (e.g., Hake, Seabirds, Marine Mammals) directly increases predation pressure on outmigrating smolts, reducing smolt survival/CPUE.
    
    #Prey ($\beta_{prey1}, \beta_{prey2}$): EXPECTED POSITIVE ($> 0$)
    #Reasoning: Higher abundance of alternative prey (e.g., Anchovy, Smelt, Krill) provides food web support or dilutes predator focus away from salmon smolts.
    
    #Temperature / SST ($\beta_{temp}$): EXPECTED NEGATIVE ($< 0$)
    #Reasoning: Warmer ocean temperatures in the NCC generally signal poorer marine survival conditions (reduced upwelling, lower prey quality).

#2. Interaction Terms (The Alternate Prey / Buffer Mechanics)
    #The expected sign of the interaction products ($\text{Predator} \times \text{Prey}$) depends on how the interaction is defined mathematically:
    
    #Predator $\times$ Alternative Prey ($\beta_{int}$): EXPECTED POSITIVE ($> 0$)
    #Reasoning (Buffering Effect): While the predator main effect is negative ($< 0$), abundant alternative prey reduces the per-capita rate at which predators consume salmon smolts. 
    #A positive interaction coefficient ($\beta_{int} > 0$) mathematically offsets/mitigates the negative main effect of the predator when alternative prey abundance is high.
    
    #Predator $\times$ Temperature ($\beta_{int}$): EXPECTED NEGATIVE ($< 0$)
    #Reasoning (Synergistic Stress): Warm SST increases metabolic demand in predators (like Hake), causing higher consumption rates. 
    #High predator abundance under warm conditions compounds the risk, yielding a negative interaction coefficient ($\beta_{int} < 0$).
