

# ==============================================================================
# Script: 02_analysis/06g_plot_latent_models_clust1_predAK_mammals.R
#       see 06a_cluster_diagrams_semnodes_v2.R for cluster diagram
# ==============================================================================

#Takeaway on fit_M1_predAK5: those 5 marine mammal indicators share a single dominant coastwide trajectory, 
#and together they provide a clean, statistically significant top-down effect on Chinook SAR ($\beta = -0.454, p = 0.011$) 
#alongside Early Marine Growth and Smolt CPUE.

# ==============================================================================
# Step 5: Expand predAK into a 5-Indicator Marine Mammal Latent Factor
# ==============================================================================
model_1B_predAK5 <- "
  # Latent Measurement Model (predAK)
  predAK =~ DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z +
            DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z +
            Killer_whales_NR_BC_adultyr__peff_Temp_z +
            harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z +
            ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z

  # Central Structural Backbone
  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
  #                               X12_DFA_EGOA_DFA_lowerTrophic_smoltyr + 
                                 predAK
                                 
# # Allow latent predAK to covary with the other 3 outer drivers (3 parameters)
#   predAK ~~ X02_DFA_SeaNettle_smoltyr + 
#             X09_DFA_HakeAge5Plus_smoltyr                                 
"

fit_1B_predAK5 <- sem(
  model_1B_predAK5,
  data         = sem_data,
  std.lv       = FALSE,
  missing      = "ML",
  fixed.x      = TRUE,
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_1B_predAK5, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
sem_output_summary_fxn(fit_1B_predAK5)

cat("\nFree parameters (npar):", fitMeasures(fit_1B_predAK5, "npar"),
    "| df:", fitMeasures(fit_1B_predAK5, "df"), "\n")

# Check top 5 modification indices to see if any indicator pair has residual correlation
cat("\nTop 5 Modification Indices:\n")
print(head(modindices(fit_1B_predAK5, sort. = TRUE), 5))


#=============
# ==============================================================================
# Step 6: Fix Variance Floor, Check Top Modification Indices, & Log Model Fit--------
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# Initialize model storage bank & comparison log (if not already in workspace)
if (!exists("sem_model_bank")) sem_model_bank <- list()
if (!exists("model_comparison_table")) model_comparison_table <- tibble()

# ------------------------------------------------------------------------------
# 1. FIT MODEL WITH SMALL POSITIVE RESIDUAL VARIANCE FLOOR (0.01)
# ------------------------------------------------------------------------------
model_1B_predAK5_noPreyAK <- "
  # Latent Measurement Model (predAK)
  predAK =~ DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z +
            DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z +
            Killer_whales_NR_BC_adultyr__peff_Temp_z +
            harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z +
            ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z

  # Fix boundary residual variance to a small positive floor (eliminates Heywood warning)
  DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z ~~ 0.01 * DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z

  # Central Structural Backbone (without preyAK)
  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 predAK
"

fit_M1_predAK5 <- sem(
  model_1B_predAK5_noPreyAK,
  data         = sem_data,
  std.lv       = TRUE,
  missing      = "ML",
  fixed.x      = TRUE,
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

# ------------------------------------------------------------------------------
# 2. INSPECT TOP 10 MODIFICATION INDICES (MI)
# ------------------------------------------------------------------------------
cat("====================================================================\n")
cat("TOP 10 MODIFICATION INDICES (M1_predAK5_noPreyAK)\n")
cat("====================================================================\n")
mi_top10 <- modindices(fit_M1_predAK5, sort. = TRUE) %>%
  select(lhs, op, rhs, mi, epc, sepc.all) %>%
  head(10)
print(as.data.frame(mi_top10), row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. HELPER TO LOG & COMPARE MODELS ACROSS ITERATIONS
# ------------------------------------------------------------------------------
log_sem_model <- function(fit_obj, model_name, notes = "") {
  sem_model_bank[[model_name]] <<- fit_obj
  r2 <- inspect(fit_obj, "r2")
  
  new_row <- tibble(
    Model     = model_name,
    npar      = as.integer(fitMeasures(fit_obj, "npar")),
    df        = as.integer(fitMeasures(fit_obj, "df")),
    Chisq     = round(as.numeric(fitMeasures(fit_obj, "chisq")), 2),
    p_Chisq   = round(as.numeric(fitMeasures(fit_obj, "pvalue")), 3),
    AIC       = round(as.numeric(fitMeasures(fit_obj, "aic")), 2),
    BIC       = round(as.numeric(fitMeasures(fit_obj, "bic")), 2),
    CFI       = round(as.numeric(fitMeasures(fit_obj, "cfi")), 3),
    SRMR      = round(as.numeric(fitMeasures(fit_obj, "srmr")), 3),
    Growth_R2 = round(as.numeric(r2["X06_DFA_IGF_mu_smoltyr"]), 3),
    CPUE_R2   = round(as.numeric(r2["X07_cpue_IntSprJunHW_smoltyr"]), 3),
    SAR_R2    = round(as.numeric(r2["X16_DFA_sar_UCsprCh_smoltyr"]), 3),
    Notes     = notes
  )
  
  if (nrow(model_comparison_table) > 0) {
    model_comparison_table <<- model_comparison_table %>%
      filter(Model != model_name) %>%
      bind_rows(new_row)
  } else {
    model_comparison_table <<- new_row
  }
  
  print(as.data.frame(model_comparison_table), row.names = FALSE)
}

log_sem_model(
  fit_M1_predAK5,
  model_name = "M1_predAK5_noPreyAK",
  notes      = "5-mammal predAK (var floor 0.01), no preyAK, fixed.x=TRUE"
)
