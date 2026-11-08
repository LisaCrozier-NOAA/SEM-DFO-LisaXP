# ==============================================================================
# Script: 02_analysis/06f_test_growth_removal_and_competitor.R
# Purpose: Test if removing Growth from SAR regression elevates Latent_AK,
#          and compare against decoupling Pink Salmon as a distinct competitor.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")

# ------------------------------------------------------------------------------
# 1. SPECIFY MODEL SYNTAXES
# ------------------------------------------------------------------------------

# Model A: Baseline with Latent_AK (4 indicators including Pink) + Growth
syntax_A <- "
  Latent_AK =~ WGOA_lower + WGOA_mid + EGOA_seabirds + X14_pinkSalmonNorthAmerica
  SAR ~ Growth + CPUE + Latent_AK + P_ocean_eff
"

# Model B: Drop Growth from SAR (Testing if AK absorbs early growth variance)
syntax_B <- "
  Latent_AK =~ WGOA_lower + WGOA_mid + EGOA_seabirds + X14_pinkSalmonNorthAmerica
  SAR ~ CPUE + Latent_AK + P_ocean_eff
"

# Model C: Decouple Pink Salmon as a separate direct competitor node
# Latent_PreyAK retains the 3 biological foraging DFAs (equating to prevent negative variance)
syntax_C <- "
  Latent_PreyAK =~ WGOA_lower + WGOA_mid + EGOA_seabirds
  SAR ~ Growth + CPUE + Latent_PreyAK + X14_pinkSalmonNorthAmerica + P_ocean_eff
"

# Model D: Decouple Pink Salmon AND Drop Growth
syntax_D <- "
  Latent_PreyAK =~ WGOA_lower + WGOA_mid + EGOA_seabirds
  SAR ~ CPUE + Latent_PreyAK + X14_pinkSalmonNorthAmerica + P_ocean_eff
"

models <- list(
  "Model_A_Full_with_Growth"   = syntax_A,
  "Model_B_Drop_Growth"        = syntax_B,
  "Model_C_Separate_Pink"      = syntax_C,
  "Model_D_Separate_Pink_NoGr" = syntax_D
)

# ------------------------------------------------------------------------------
# 2. FIT MODELS AND EXTRACT COMPARATIVE METRICS
# ------------------------------------------------------------------------------
comparison_results <- list()

for (m_name in names(models)) {
  fit <- tryCatch(
    sem(
      models[[m_name]], 
      data         = df_fit, 
      std.lv       = TRUE, 
      missing      = "ML", 
      bounds       = "pos.var",
      optim.method = "BFGS",
      control      = list(iter.max = 1000)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !lavInspect(fit, "converged")) {
    cat(sprintf("[-] %s: Failed to converge.\n", m_name))
    next
  }
  
  pe <- parameterEstimates(fit, standardized = TRUE)
  
  # Extract structural paths predicting SAR
  sar_paths <- pe %>% 
    filter(lhs == "SAR", op == "~") %>%
    mutate(path_desc = sprintf("%s: beta=%.3f (p=%.3f)", rhs, std.all, pvalue)) %>%
    pull(path_desc) %>%
    paste(collapse = " | ")
  
  comparison_results[[m_name]] <- tibble(
    Model     = m_name,
    AIC       = round(as.numeric(fitMeasures(fit, "aic")), 2),
    p_Chi2    = round(as.numeric(fitMeasures(fit, "pvalue")), 4),
    CFI       = round(as.numeric(fitMeasures(fit, "cfi")), 3),
    SRMR      = round(as.numeric(fitMeasures(fit, "srmr")), 3),
    SAR_R2    = round(as.numeric(inspect(fit, "r2")["SAR"]), 3),
    SAR_Paths = sar_paths
  )
}

summary_table <- bind_rows(comparison_results)
cat("\n====================================================================\n")
cat("MODEL COMPARISON: TESTING GROWTH REMOVAL AND DECOUPLED COMPETITION\n")
cat("====================================================================\n")
print(as.data.frame(summary_table))

write_csv(summary_table, file.path(data_dir, "obj4_test_growth_and_pink_architectures.csv"))


# ==============================================================================
# Model E: Pruned Observed Multi-Stage Bottleneck Architecture
# ==============================================================================

syntax_E <- "
  SAR ~ CPUE + X14_pinkSalmonNorthAmerica + P_ocean_eff
"

fit_E <- sem(
  syntax_E, 
  data = df_fit, 
  missing = "ML", 
  control = list(iter.max = 1000)
)

summary(fit_E, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)


# ==============================================================================
# Testing Decapods as an Independent Node in the Clean Model E Backbone----------
# ==============================================================================

syntax_decapod_test <- "
  SAR ~ CPUE + X14_pinkSalmonNorthAmerica + P_ocean_eff + X12_biomassDecaShelfSum
"

fit_decapod <- sem(syntax_decapod_test, data = df_fit, missing = "ML")
summary(fit_decapod, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
