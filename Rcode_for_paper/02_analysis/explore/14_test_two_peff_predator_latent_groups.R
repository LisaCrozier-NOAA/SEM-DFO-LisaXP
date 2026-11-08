# ==============================================================================
# Script: 02_analysis/14_test_two_peff_predator_latent_groups.R
# Purpose: Test latent models for Group 1 (Mammals) and Group 2 (Demersal),
#          prune non-significant indicators, and compare against the observed baseline.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# Combine into your standardized dataset (df_fit) with CPUE and Pink
# ------------------------------------------------------------------------------
# 1. SPECIFY MODELS
# ------------------------------------------------------------------------------

# Model 1: Group 1 Alone (Mammal Latent)
syntax_g1_latent <- "
  Latent_Mammal =~ DFA_Steller_s_l_WS_adultyr + 
                   harbour_seal_DFO_WCVI_smoltyr + 
                   Killer_whales_NR_BC_adultyr + 
                   ssl_seak_pup_pred_smoltyr

  SAR ~ CPUE + Pink + Latent_Mammal
"

# Model 2: Group 2 Alone (Demersal / Shark Latent - Equated Loadings)
syntax_g2_latent <- "
  Latent_Demersal =~ a*DFA_sleeperSharkBSAI_predAK_smoltyr + 
                     a*sablefishAK_biomass3plus_smoltyr

  SAR ~ CPUE + Pink + Latent_Demersal
"

# Model 3: Dual Latent Model (Both Marine Mammal & Demersal Latents Co-varying)
syntax_dual_latent <- "
  Latent_Mammal =~ DFA_Steller_s_l_WS_adultyr + 
                   harbour_seal_DFO_WCVI_smoltyr + 
                   Killer_whales_NR_BC_adultyr + 
                   ssl_seak_pup_pred_smoltyr

  Latent_Demersal =~ a*DFA_sleeperSharkBSAI_predAK_smoltyr + 
                     a*sablefishAK_biomass3plus_smoltyr

  Latent_Mammal ~~ Latent_Demersal

  SAR ~ CPUE + Pink + Latent_Mammal + Latent_Demersal
"

# Model 4: Top Observed Benchmark (Single Sablefish Peff alone)
syntax_obs_sable <- "
  SAR ~ CPUE + Pink + sablefishAK_biomass3plus_smoltyr
"

# ------------------------------------------------------------------------------
# 2. FIT AND EXTRACT DIAGNOSTICS
# ------------------------------------------------------------------------------
run_peff_model <- function(syntax, name) {
  fit <- tryCatch(
    sem(syntax, data = df_preds_scaled, missing = "ML", bounds = "pos.var",
        optim.method = "BFGS", control = list(iter.max = 2000)),
    error = function(e) NULL
  )
  if (is.null(fit) || !lavInspect(fit, "converged")) {
    return(tibble(Model = name, AIC = NA, CFI = NA, SRMR = NA, SAR_R2 = NA, 
                  Heywood = "FAILED", Paths = "Did not converge"))
  }
  pe <- parameterEstimates(fit, standardized = TRUE)
  loadings  <- pe %>% filter(op == "=~")
  variances <- pe %>% filter(op == "~~", lhs == rhs, lhs != "SAR")
  has_heywood <- any(variances$est <= 0.001) || any(abs(loadings$std.all) > 1.0)
  
  paths <- pe %>%
    filter(lhs == "SAR", op == "~") %>%
    mutate(desc = sprintf("%s: b=%.3f (p=%.3f)", rhs, std.all, pvalue)) %>%
    pull(desc) %>%
    paste(collapse = " | ")
  
  tibble(
    Model   = name,
    AIC     = round(as.numeric(fitMeasures(fit, "aic")), 2),
    CFI     = round(as.numeric(fitMeasures(fit, "cfi")), 3),
    SRMR    = round(as.numeric(fitMeasures(fit, "srmr")), 3),
    SAR_R2  = round(as.numeric(inspect(fit, "r2")["SAR"]), 3),
    Heywood = ifelse(has_heywood, "YES", "NO"),
    Paths   = paths
  )
}

# Attach CPUE, Pink, and SAR to df_preds for fitting
df_preds_scaled <- df_preds %>%
  mutate(SAR  = df_fit$SAR,
         CPUE = df_fit$CPUE,
         Pink = df_fit$Pink) %>%
  mutate(across(everything(), ~ as.numeric(scale(.x))))

models_list <- list(
  "M1_Mammal_Latent"    = syntax_g1_latent,
  "M2_Demersal_Latent"  = syntax_g2_latent,
  "M3_Dual_Latent"      = syntax_dual_latent,
  "M4_Observed_Sable"   = syntax_obs_sable
)

benchmark_results <- map2_df(models_list, names(models_list), run_peff_model) %>%
  mutate(delta_AIC = round(AIC - min(AIC, na.rm = TRUE), 2)) %>%
  arrange(AIC)

cat("\n====================================================================\n")
cat("BENCHMARK: CLUSTER-DERIVED PEFF LATENTS VS. OBSERVED BENCHMARK\n")
cat("====================================================================\n")
print(as.data.frame(benchmark_results))


#LISA CODE====================
#manual pruning-------------------