library(tidyverse)

# 1. Mandatory AK Predator: Must have a negative main predator effect on SAR (x16_sar)
mandatory_ak_topdown_pairs <- all_params_df %>%
  filter(
    lhs == "x16_sar",
    term_type == "Main Predator Effect",
    est < 0,
    pvalue <= 0.05
  ) %>%
  pull(pair_id) %>%
  unique()

# 2. Blacklist: Identify any pairs with a POSITIVE main predator effect anywhere (CPUE or SAR)
pairs_with_positive_predators <- all_params_df %>%
  filter(
    term_type == "Main Predator Effect",
    est > 0
  ) %>%
  pull(pair_id) %>%
  unique()

# 3. Valid Pairs: Has mandatory negative AK predator AND zero positive predator effects anywhere
strictly_valid_topdown_pairs <- setdiff(mandatory_ak_topdown_pairs, pairs_with_positive_predators)

# 4. Filter model summaries for global SEM fit
top_down_summary_results <- summary_results %>%
  filter(
    pair_id %in% strictly_valid_topdown_pairs,
    pvalue >= 0.05,    # Global SEM chi-square p-value
    cfi >= 0.90        # Global CFI
  ) %>%
  arrange(aic) %>%
  mutate(delta_aic = aic - min(aic, na.rm = TRUE))

# 5. Extract full parameter breakdown
top_down_parameters <- all_params_df %>%
  filter(pair_id %in% top_down_summary_results$pair_id) %>%
  left_join(top_down_summary_results %>% select(pair_id, delta_aic), by = "pair_id") %>%
  arrange(delta_aic, pair_id)

# -----------------------------------------------------------------------------
# Print Outputs
# -----------------------------------------------------------------------------
cat("\n====================================================================\n")
cat(" TOP VALID TWO-STAGE MODELS WITH MANDATORY NEGATIVE AK PREDATOR     \n")
cat(" (pvalue >= 0.05, CFI >= 0.90, AK Pred < 0, No Positive Predators) \n")
cat("====================================================================\n")

print(top_down_summary_results %>% as_tibble(), n = Inf)

cat("\n====================================================================\n")
cat(" PARAMETER BREAKDOWN                                                \n")
cat("====================================================================\n")

print(
  top_down_parameters %>%
    select(pair_id, delta_aic, stage = lhs, term, resolved_species, term_type, est, pvalue, effect_direction) %>%
    as_tibble(),
  n = Inf
)




###########################################
###########################################
# ====================================================================
# PARAMETER BREAKDOWN                                                
# ====================================================================
# pair_id delta_aic stage                    term                     resolved_species      term_type                      est   pvalue effect_direction              
# 1      39         0 x07_dfa_cpue_intsprjunhw x09_dfa_hakeage5plus     x09_dfa_hakeage5plus  Main Predator Effect       -0.603  0.000211 Top-Down Predation (-)        
# 2      39         0 x16_sar                  x07_dfa_cpue_intsprjunhw Early CPUE Link       Stage Link (CPUE -> SAR)    0.390  0.00786  Positive Link                 
# 3      39         0 x16_sar                  x11_ssl_seak_pup_pred    x11_ssl_seak_pup_pred Main Predator Effect       -0.991  0.000899 Top-Down Predation (-)        
# 4      39         0 x16_sar                  ak_int_1                 x13_stka_herr_matbiom Alternate Prey Interaction  0.0625 0.00795  Buffering / Prey-Switching (+)
# 5      39         0 x16_sar                  ak_int_2                 x13_mid_il_capelin    Alternate Prey Interaction  0.0744 0.00786  Buffering / Prey-Switching (+)
# 6      39         0 x16_sar                  ak_int_3                 x12_egoa.krill        Alternate Prey Interaction  0.0593 0.0182   Buffering / Prey-Switching (+)
