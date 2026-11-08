library(tidyverse)
library(lavaan)

# ------------------------------------------------------------------------------
# 1. SETUP DATA & EXACT COLUMN IDENTIFIERS
# ------------------------------------------------------------------------------
data_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/sem_altprey_data_1998_2021_10212026.csv"
raw_df    <- read.csv(data_path, stringsAsFactors = FALSE)
raw_df    <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

# Filter 24-year complete dataset
df_24 <- raw_df %>% filter(complete.cases(.))

scale_bounds <- function(x, a = 0.2, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

# Base min-max standardized dataframe (0 to 1 for regressions)
df_base <- df_24 %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_bounds(.x, 0.0, 1.0)))

# Fixed anchors
var_growth    <- "X06_DFA_IGF_mu_smoltyr"
var_prey_ncc  <- "X05_DFA_abundSardine_smoltyr"
var_cpue      <- "X07_cpue_IntSprJunHW_smoltyr"
var_hake      <- "X09_DFA_HakeAge5Plus_smoltyr"
var_sar       <- "X16_DFA_sar_UCsprCh_smoltyr"
var_prey_ak   <- "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"

# ------------------------------------------------------------------------------
# 2. GENERALIZED FIT ENGINE: DAG 1B + MODERATION
# ------------------------------------------------------------------------------
fit_dag1b_moderated <- function(pred_col, 
                                altprey_cols = NULL, 
                                sst_col = NULL, 
                                include_prey_ak = TRUE, 
                                df_in = df_base, 
                                df_raw_in = df_24) {
  
  if (!pred_col %in% names(df_in)) {
    warning("Predator not found in df_base: ", pred_col)
    return(NULL)
  }
  
  # A. Calculate Effective Predator (P_eff)
  p_vec <- df_in[[pred_col]]
  
  # Temperature gating: scale [0.2, 1.0]
  if (!is.null(sst_col) && sst_col != "None" && sst_col %in% names(df_raw_in)) {
    gate_temp <- scale_bounds(df_raw_in[[sst_col]], a = 0.2, b = 1.0)
    p_vec <- p_vec * gate_temp
  }
  
  # Alternate prey buffering gate: 1 - 1.0 * AltPrey_scaled [0.2, 1.0]
  if (!is.null(altprey_cols) && length(altprey_cols) > 0) {
    valid_prey <- intersect(altprey_cols, names(df_in))
    if (length(valid_prey) > 0) {
      if (length(valid_prey) == 1) {
        prey_sum <- df_in[[valid_prey[1]]]
      } else {
        prey_sum <- rowSums(df_in[, valid_prey, drop = FALSE], na.rm = TRUE)
      }
      gate_alt <- 1.0 - 1.0 * scale_bounds(prey_sum, a = 0.2, b = 1.0)
      p_vec <- p_vec * gate_alt
    }
  }
  
  df_in$P_eff <- p_vec
  
  # B. Build DAG 1B Structural Syntax
  if (include_prey_ak && var_prey_ak %in% names(df_in)) {
    sem_syntax <- paste0(
      var_growth, " ~ ", var_prey_ncc, "\n",
      var_cpue,   " ~ ", var_hake, "\n",
      var_sar,    " ~ ", var_growth, " + ", var_cpue, " + ", var_prey_ak, " + P_eff"
    )
  } else {
    sem_syntax <- paste0(
      var_growth, " ~ ", var_prey_ncc, "\n",
      var_cpue,   " ~ ", var_hake, "\n",
      var_sar,    " ~ ", var_growth, " + ", var_cpue, " + P_eff"
    )
  }
  
  # C. Fit with Lavaan
  fit <- tryCatch(
    sem(sem_syntax, data = df_in, missing = "ML", warn = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  
  pe  <- parameterEstimates(fit, standardized = TRUE)
  rsq <- inspect(fit, "r2")
  
  pred_est   <- pe %>% filter(lhs == var_sar, rhs == "P_eff", op == "~")
  growth_est <- pe %>% filter(lhs == var_sar, rhs == var_growth, op == "~")
  preyak_est <- pe %>% filter(lhs == var_sar, rhs == var_prey_ak, op == "~")
  
  tibble(
    predator         = pred_col,
    sst_driver       = ifelse(is.null(sst_col), "None", sst_col),
    alt_prey         = ifelse(is.null(altprey_cols), "None", paste(altprey_cols, collapse = " + ")),
    has_prey_ak      = include_prey_ak,
    aic              = round(as.numeric(fitMeasures(fit, "aic")), 2),
    r2_sar           = round(as.numeric(rsq[[var_sar]]), 3),
    beta_peff_std    = round(pred_est$std.all[1], 3),
    p_peff           = round(pred_est$pvalue[1], 4),
    beta_growth_std  = round(growth_est$std.all[1], 3),
    p_growth         = round(growth_est$pvalue[1], 4),
    beta_preyak_std  = ifelse(include_prey_ak && nrow(preyak_est) > 0, round(preyak_est$std.all[1], 3), NA_real_),
    p_preyak         = ifelse(include_prey_ak && nrow(preyak_est) > 0, round(preyak_est$pvalue[1], 4), NA_real_),
    pred_sign_neg    = (pred_est$std.all[1] < 0),
    growth_sign_pos  = (growth_est$std.all[1] > 0 && growth_est$pvalue[1] < 0.05)
  )
}

# ------------------------------------------------------------------------------
# 3. RUN THE TEST BATTERY (Sablefish & Sleeper Shark)
# ------------------------------------------------------------------------------
test_runs <- list(
  # 1. Sablefish: Baseline vs Moderated (with & without Prey AK)
  list(label = "Sablefish - Raw (With PreyAK)",
       pred = "X15_sablefishAK_biomass3plus_smoltyr", prey = NULL, sst = NULL, ak = TRUE),
  list(label = "Sablefish - Moderated (With PreyAK)",
       pred = "X15_sablefishAK_biomass3plus_smoltyr", 
       prey = c("X13_pollockGOA_age1and2_smoltyr", "sitkaHerring_EGoA_smoltyr"), 
       sst  = "SST_NCC_aprjun_DFObccm_adultyr", ak = TRUE),
  list(label = "Sablefish - Moderated (Dropped PreyAK)",
       pred = "X15_sablefishAK_biomass3plus_smoltyr", 
       prey = c("X13_pollockGOA_age1and2_smoltyr", "sitkaHerring_EGoA_smoltyr"), 
       sst  = "SST_NCC_aprjun_DFObccm_adultyr", ak = FALSE),
  
  # 2. Sleeper Shark: Baseline vs Moderated (with & without Prey AK)
  list(label = "Sleeper Shark - Raw (With PreyAK)",
       pred = "X15_DFA_sleeperSharkBSAI_predAK_smoltyr", prey = NULL, sst = NULL, ak = TRUE),
  list(label = "Sleeper Shark - Moderated (With PreyAK)",
       pred = "X15_DFA_sleeperSharkBSAI_predAK_smoltyr", 
       prey = c("capelin_WGoA_smoltyr"), 
       sst  = "X21_sst_egoa_junjulaug_smoltyr", ak = TRUE),
  list(label = "Sleeper Shark - Moderated (Dropped PreyAK)",
       pred = "X15_DFA_sleeperSharkBSAI_predAK_smoltyr", 
       prey = c("capelin_WGoA_smoltyr"), 
       sst  = "X21_sst_egoa_junjulaug_smoltyr", ak = FALSE)
)

compiled_results <- list()

for (tr in test_runs) {
  res <- fit_dag1b_moderated(
    pred_col        = tr$pred, 
    altprey_cols    = tr$prey, 
    sst_col         = tr$sst, 
    include_prey_ak = tr$ak
  )
  if (!is.null(res)) {
    compiled_results[[length(compiled_results) + 1]] <- res %>% mutate(scenario = tr$label)
  }
}

final_comp_table <- bind_rows(compiled_results) %>%
  select(scenario, aic, r2_sar, beta_peff_std, p_peff, beta_growth_std, p_growth, beta_preyak_std, p_preyak, pred_sign_neg, growth_sign_pos)

cat("\n--- DAG 1B MODERATION COMPARISON TABLE ---\n")
print(as.data.frame(final_comp_table))
