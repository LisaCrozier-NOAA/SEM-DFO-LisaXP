# ==============================================================================
# Script: 11_compare_previous_table4_vs_minmax.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Compare previous Goal 2 (Table 4) baseline SEM indicators against 
#          the new 5-Guild Min-Max Gating Attenuation models.
# Output: Rcode_for_paper/Routput_for_paper/tables/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(gt)
  library(stringr)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & LOAD DATASETS
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
tbl_out_dir  <- file.path(proj_dir, "tables")
data_out_dir <- file.path(master_out, "data")

# 1. Load Previous Table 4 Baseline Diagnostics
t4_summary_path <- file.path(tbl_out_dir, "Table4_Indicator_Diagnostics_dAIC3.csv")
t4_eval_path    <- file.path(data_out_dir, "Goal2_Indicator_Evaluations_dAIC3.csv")

if (!file.exists(t4_summary_path) && !file.exists(t4_eval_path)) {
  stop("Missing Table 4 baseline files in 'tables/' or 'data/'. Run script 02 first.")
}

t4_summary <- if (file.exists(t4_summary_path)) read_csv(t4_summary_path, show_col_types = FALSE) else NULL

# 2. Load Min-Max Sweep Results (or Checkpoint if still running)
years_included<-"19yr"
sweep_file <- file.path(data_out_dir, paste0("Goal3_NCC_CPUE_Sweep_Results_minmax_combo_size_2_", years_included, ".csv"))
#sweep_file <- file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_24yr_noNA_akferris.csv")
#sweep_path <- file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_19yr_noNA_dougsmoothed.csv")
#chk_file   <- file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_CHECKPOINT.csv")


# sweep_path <- if (file.exists(sweep_file)) sweep_file else chk_file
# sweep_path <- sweep_file
# if (!file.exists(sweep_path)) {
#   stop("No Min-Max sweep results or checkpoint files found in 'data/'.")
# }

cat(sprintf("Loading Min-Max Attenuation Sweep results from: %s\n", basename(sweep_path)))
minmax_results <- read_csv(sweep_file, show_col_types = FALSE)

minmax_results %>% filter()
#did altprey/sst flip signs?
#did altprey/sst produce a better model?
#what are the final selected models?
nrow(minmax_results)
neg_selected<-minmax_results %>% filter(neg_delta_aic< -2)
nrow(neg_selected) #6/19 (16 modelsdaic<0, but not really justified)
neg_selected <- minmax_results %>% filter(neg_pvalue < 0.05,neg_delta_aic< -2)
print(neg_selected %>% select(predator,neg_prey,base_beta,neg_beta_peff,base_r2,neg_final_r2))
predator                               neg_prey                                                           base_beta neg_beta_peff base_r2 neg_final_r2
<chr>                                  <chr>                                                                  <dbl>         <dbl>   <dbl>        <dbl>
1 X09_DFA_HakeAge5Plus_smoltyr           X12_DFA_biomassEuphShelfSum_smoltyr + X09_DFA_HakeAge5Plus_smoltyr   -0.535         -2.38   0.445         0.683
2 X10_Californian_s_l_2yrLead_WS_smoltyr X01_habCompInd_smoltyr + X04_marketsquid_GAM_smoltyr                 -0.0874        -0.914  0.0128        0.206
> 
# ------------------------------------------------------------------------------
# STEP 1: CLEAN AND HARMONIZE PREDATOR INDICATOR NAMES FOR JOINING
# ------------------------------------------------------------------------------
# Clean Table 4 baseline records
t4_predators <- t4_summary %>%
  filter(SEM_Node %in% c("PredNCC", "PredAK")) %>%
  mutate(
    # Normalize indicator string to match data file column formats
    pred_clean = case_when(
      grepl("X09_DFA_HakeAge5Plus", Indicator) ~ "X09_DFA_HakeAge5Plus",
      grepl("X15_PacificCodBiomass", Indicator) ~ "X15_PacificCodBiomass_predAK",
      grepl("X15_sablefishBiomass", Indicator) ~ "X15_sablefishBiomass_predAK",
      grepl("X15_sablefishRecruitment", Indicator) ~ "X15_sablefishRecruitment_predAK",
      grepl("X15_spinyDogfishBSAI", Indicator) ~ "X15_spinyDogfishBSAI_predAK",
      grepl("X15_DFA_sleeperSharks", Indicator) ~ "X15_DFA_sleeperSharks",
      grepl("X10_Harbour_s", Indicator) ~ "X10_Harbour_s_2yrLead_WS",
      TRUE ~ Indicator
    ),
    # Capture original baseline beta (Supported vs Flipped)
    Baseline_Beta = if_else(!is.na(Beta_Supp), Beta_Supp, Beta_Flip),
    Baseline_Status = if_else(!is.na(Beta_Supp), "Supported (-)", "Flipped (+)")
  )

# Match Min-Max results back to Table 4 baseline indicators
minmax_matched <- minmax_results %>%
  mutate(
    pred_clean = case_when(
      grepl("X09_DFA_HakeAge5Plus", predator) ~ "X09_DFA_HakeAge5Plus",
      grepl("X15_PacificCodBiomass", predator) ~ "X15_PacificCodBiomass_predAK",
      grepl("X15_sablefishBiomass", predator) ~ "X15_sablefishBiomass_predAK",
      grepl("X15_sablefishRecruitment", predator) ~ "X15_sablefishRecruitment_predAK",
      grepl("X15_spinyDogfish", predator) ~ "X15_spinyDogfishBSAI_predAK",
      grepl("X15_DFA_sleeperSharks", predator) ~ "X15_DFA_sleeperSharks",
      grepl("X10_Harbour_s", predator) ~ "X10_Harbour_s_2yrLead_WS",
      TRUE ~ predator
    )
  )

# ------------------------------------------------------------------------------
# STEP 2: BUILD SIDE-BY-SIDE COMPARISON TABLE
# ------------------------------------------------------------------------------
comparison_table <- t4_predators %>%
  inner_join(minmax_matched, by = "pred_clean", relationship = "many-to-many") %>%
  select(
    SEM_Node,
    Indicator,
    predator_full = predator,
    N_Models,
    Baseline_Status,
    Baseline_Beta,
    Selected_Gates = selected_prey,
    Fitted_Weights = fitted_weights,
    MinMax_Beta_Peff = beta_peff,
    Pvalue = pvalue,
    MinMax_AIC = aic,
    Base_AIC = base_aic,
    Delta_AIC = delta_aic
  ) %>%
  arrange(SEM_Node, MinMax_AIC)

write_csv(comparison_table, file.path(tbl_out_dir, "Table4_vs_MinMax_Predator_Comparison_19yr.csv"))

comparison_table
# A tibble: 19 × 13
# SEM_Node Indicator                       predator_full                     N_Models Baseline_Status Baseline_Beta Selected_Gates Fitted_Weights MinMax_Beta_Peff  Pvalue MinMax_AIC Base_AIC Delta_AIC
# 1 PredAK   X15_sablefishBiomass_predAK     X15_sablefishBiomass_predAK_smol…        5 Flipped (+)             0.717 X13_pollock_a… 1; 1                     -4.56  2.34e-8     -7.10    10.9     -18.0  
# 2 PredAK   X10_Harbour_s_2yrLead_WS        X10_Harbour_s_2yrLead_WS_smoltyr         3 Supported (-)          -0.336 X01_habCompIn… 1; 0.3                   -1.12  1.67e-5     -1.75     8.40    -10.1  
# 3 PredAK   X15_sablefishBiomass_predAK     X15_sablefishBiomass_predAK_adul…        5 Flipped (+)             0.717 X13_pollock_a… 1; 1                     -2.27  4.60e-5     -0.771   10.7     -11.4  
# 4 PredAK   X10_Harbour_s_2yrLead_WS        X10_Harbour_s_2yrLead_WS_adultyr         3 Supported (-)          -0.336 X01_habCompIn… 0.6; 1                   -0.874 6.15e-5     -0.483    5.88     -6.36 
# 5 PredAK   X15_DFA_sleeperSharks           X15_DFA_sleeperSharks_adultyr           13 Supported (-)          -0.918 X14_pinkSalmo… 0.4; 1                   -0.715 7.01e-5     -0.352    4.44     -4.79 
# 6 PredAK   X15_DFA_sleeperSharks           X15_DFA_sleeperSharks_smoltyr           13 Supported (-)          -0.918 X13_mid_il_ca… 1                        -0.623 2.38e-4      0.918    1.04     -0.122
# 7 PredAK   X15_sablefishRecruitment_predAK X15_sablefishRecruitment_predAK_…        2 Supported (-)          -0.466 X13_pollock_a… 1; 1                     -4.00  1.95e-3      3.30    11.0      -7.65 
# 8 PredAK   X15_spinyDogfishBSAI_predAK     X15_spinyDogfishBSAI_predAK_smol…        2 Flipped (+)             0.405 X13_mid_il_ca… 1; 1                     -0.933 2.39e-3      3.55     9.25     -5.71 
# 9 PredAK   X15_PacificCodBiomass_predAK    X15_PacificCodBiomass_predAK_adu…        4 Flipped (+)             0.429 X13_pollock_a… 1; 1                     -0.794 4.67e-3      4.37    10.9      -6.50 
# 10 PredAK   X15_spinyDogfishBSAI_predAK     X15_spinyDogfishBSAI_predAK_adul…        2 Flipped (+)             0.405 X13_pollock_a… 1; 1                     -0.868 1.98e-2      6.25    10.2      -3.96 
# 11 PredAK   X15_spinyDogfishBSAI_predAK     X15_spinyDogfishGoA_predAK_smolt…        2 Flipped (+)             0.405 X13_pollock_a… 1; 1                     -1.71  2.74e-2      6.69    11.0      -4.30 
# 12 PredAK   X15_PacificCodBiomass_predAK    X15_PacificCodBiomass_predAK_smo…        4 Flipped (+)             0.429 None (Unatten… None                      0.347 5.41e-2      7.62     7.62      0    
# -----------------------------------------------------------------------------
gt_comp_df <- comparison_table %>%
  select(
    Indicator,
    Baseline_Status,
    Baseline_Beta,
    Selected_Gates,
    MinMax_Beta_Peff,
    MinMax_AIC,
    Delta_AIC 
  )

gt_comp <- gt(gt_comp_df) %>%
  tab_header(
    title = md("**19yr Comparative Diagnostic: Additive Table 4 SEMs vs. Min-Max Gated Models**"),
    subtitle = "Evaluating whether spatial overlap / alternate prey gating resolves baseline flipped predator signs."
  ) %>%
  tab_spanner(
    label = md("**Previous Baseline (Table 4 Additive)**"),
    columns = c(Baseline_Status, Baseline_Beta)
  ) %>%
  tab_spanner(
    label = md("**New Min-Max Attenuation Model**"),
    columns = c(Selected_Gates, MinMax_Beta_Peff, MinMax_AIC, Delta_AIC)
  ) %>%
  cols_label(
    Indicator        = md("**Predator Indicator**"),
    Baseline_Status  = md("**Status**"),
    Baseline_Beta    = md("**Base β**"),
    Selected_Gates   = md("**Selected Gate(s)**"),
    MinMax_Beta_Peff = md("**Gated β**"),
    MinMax_AIC       = md("**AIC**"),
    Delta_AIC        = md("**ΔAIC**")
  ) %>%
  # --- ROUND ALL NUMERIC COLUMNS TO 2 DECIMAL PLACES ---
  fmt_number(
    columns = where(is.numeric),
    decimals = 2
  ) %>%
  cols_align(align = "left", columns = c(Indicator, Selected_Gates)) %>%
  cols_align(align = "center", columns = c(Baseline_Status, Baseline_Beta, MinMax_Beta_Peff, MinMax_AIC, Delta_AIC)) %>%
  sub_missing(columns = everything(), missing_text = "-") %>%
  tab_style(
    style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#2E7D32", weight = "bold")),
    locations = cells_body(
      columns = c(MinMax_Beta_Peff, Selected_Gates),
      rows = Baseline_Status == "Flipped (+)" & MinMax_Beta_Peff < 0
    )
  ) %>%
  tab_options(
    table.font.size = px(10),
    heading.title.font.size = px(12),
    column_labels.font.weight = "bold",
    data_row.padding = px(3)
  )

print(gt_comp)

#gtsave(gt_comp, file.path(tbl_out_dir, "Table4_vs_MinMax_Predator_Comparison_24yr.html"))
gtsave(gt_comp, file.path(tbl_out_dir, "Table4_vs_MinMax_Predator_Comparison_19yr.html"))
message("Comparison complete! Results written to: ", tbl_out_dir)

library(webshot2) # Required by gt to render PNGs

# Save as high-res PNG image (adjust vwidth/vheight as needed)
gtsave(
  data = gt_comp,
  filename = file.path(tbl_out_dir, "Table4_vs_MinMax_Predator_Comparison_19yr.png"),
#  filename = file.path(tbl_out_dir, "Table4_vs_MinMax_Predator_Comparison_24yr.png"),
  vwidth = 1000,
  vheight = 800,
  zoom = 2 # Increases image crispness for slides
)
