

#3. R Script: Integrated DAG 1B Moderation Sweep
#This script runs a test sweep over your candidate Alaska/NCC predators and their alternate prey/temperature drivers, 
#evaluating:Whether adding the moderating gate turns the predator $\beta$ negative ($\beta < 0$).
#Whether the Growth (IGF) -> SAR link remains robust and positive ($\beta > 0$, $p < 0.05$).
#The comparison between Keeping PreyAK vs. Dropping PreyAK.


library(tidyverse)
library(lavaan)

# ------------------------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------------------------
data_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/sem_altprey_data_1998_2021_10212026.csv"
lookup_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/Routput_for_paper/data/guild_predator_lookup_table_24yr.csv"

raw_df <- read.csv(data_path, stringsAsFactors = FALSE)
raw_df <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

# Filter 24-year complete dataset
df_24 <- raw_df %>% filter(complete.cases(.))

scale_bounds <- function(x, a = 0.2, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

# Base min-max standardized dataframe
df_base <- df_24 %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_bounds(.x, 0.0, 1.0)))

# Fixed anchors
var_growth    <- "X06_DFA_IGF_mu"
var_prey_ncc  <- "X05_DFA_abundSardine"
var_cpue      <- "X07_DFA_cpue_IntSprJunHW"
var_hake      <- "X09_DFA_HakeAge5Plus_smoltyr"
var_sar       <- "X16_DFA_sar_UCsprCh_smoltyr"

# Candidate Prey AK series to test
var_prey_ak   <- "X12_DFA_EGOA_DFA_lowerTrophic" # or your preferred AK prey series

# ------------------------------------------------------------------------------
# 2. Evaluation Engine: Integrated DAG 1B with Moderated PredAK
# ------------------------------------------------------------------------------
fit_dag1b_moderated <- function(pred_col, 
                                altprey_cols = NULL, 
                                sst_col = NULL, 
                                include_prey_ak = TRUE, 
                                df_in = df_base, 
                                df_raw_in = df_24) {
  
  if (!pred_col %in% names(df_in)) return(NULL)
  
  # A. Calculate P_eff
  p_vec <- df_in[[pred_col]]
  
  # Temperature gating (if specified)
  if (!is.null(sst_col) && sst_col %in% names(df_raw_in)) {
    gate_temp <- scale_bounds(df_raw_in[[sst_col]], a = 0.2, b = 1.0)
    p_vec <- p_vec * gate_temp
  }
  
  # Alternate prey buffering gate (if specified)
  if (!is.null(altprey_cols) && length(altprey_cols) > 0) {
    valid_prey <- intersect(altprey_cols, names(df_in))
    if (length(valid_prey) > 0) {
      prey_sum <- rowSums(df_in[, valid_prey, drop = FALSE], na.rm = TRUE)
      gate_alt <- 1.0 - 1.0 * scale_bounds(prey_sum, a = 0.2, b = 1.0)
      p_vec <- p_vec * gate_alt
    }
  }
  
  df_in$P_eff <- p_vec
  
  # B. Build DAG 1B Syntax
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
  
  # C. Fit Model
  fit <- tryCatch(
    sem(sem_syntax, data = df_in, missing = "ML", warn = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  
  pe <- parameterEstimates(fit, standardized = TRUE)
  rsq <- inspect(fit, "r2")
  
  pred_est   <- pe %>% filter(lhs == var_sar, rhs == "P_eff", op == "~")
  growth_est <- pe %>% filter(lhs == var_sar, rhs == var_growth, op == "~")
  
  tibble(
    predator         = pred_col,
    sst_driver       = ifelse(is.null(sst_col), "None", sst_col),
    alt_prey         = ifelse(is.null(altprey_cols), "None", paste(altprey_cols, collapse = " + ")),
    has_prey_ak      = include_prey_ak,
    aic              = fitMeasures(fit, "aic"),
    r2_sar           = as.numeric(rsq[[var_sar]]),
    beta_peff_std    = pred_est$std.all[1],
    p_peff           = pred_est$pvalue[1],
    beta_growth_std  = growth_est$std.all[1],
    p_growth         = growth_est$pvalue[1],
    pred_sign_ok     = (pred_est$std.all[1] < 0 && pred_est$pvalue[1] < 0.10),
    growth_sign_ok   = (growth_est$std.all[1] > 0 && growth_est$pvalue[1] < 0.10)
  )
}

# ------------------------------------------------------------------------------
# 3. Test Example: Sablefish & Sleeper Shark (Raw vs. Moderated, With/Without PreyAK)
# ------------------------------------------------------------------------------
test_cases <- list(
  # Sablefish with pollock/herring alternate prey
  list(pred = "X15_sablefishAK_biomass3plus_smoltyr",
       prey = c("X13_pollockGOA_age1and2_smoltyr", "sitkaHerring_EGoA_smoltyr"),
       temp = "SST_NCC_aprjun_DFObccm_adultyr"),
  
  # Sleeper Shark with Capelin alternate prey
  list(pred = "X15_DFA_sleeperSharkBSAI_predAK_smoltyr",
       prey = c("capelin_WGoA_smoltyr"),
       temp = "X21_sst_egoa_junjulaug_smoltyr")
)

results <- list()

for (tc in test_cases) {
  # 1. Raw Baseline (With PreyAK)
  results[[length(results) + 1]] <- fit_dag1b_moderated(tc$pred, altprey_cols = NULL, sst_col = NULL, include_prey_ak = TRUE) %>%
    mutate(model_desc = "Raw Baseline (With PreyAK)")
  
  # 2. Moderated (With PreyAK)
  results[[length(results) + 1]] <- fit_dag1b_moderated(tc$pred, altprey_cols = tc$prey, sst_col = tc$temp, include_prey_ak = TRUE) %>%
    mutate(model_desc = "Moderated (With PreyAK)")
  
  # 3. Moderated (Dropped PreyAK)
  results[[length(results) + 1]] <- fit_dag1b_moderated(tc$pred, altprey_cols = tc$prey, sst_col = tc$temp, include_prey_ak = FALSE) %>%
    mutate(model_desc = "Moderated (Dropped PreyAK)")
}

comparison_table <- bind_rows(results) %>%
  select(predator, model_desc, aic, r2_sar, beta_peff_std, p_peff, beta_growth_std, p_growth, pred_sign_ok, growth_sign_ok)

print(as.data.frame(comparison_table))