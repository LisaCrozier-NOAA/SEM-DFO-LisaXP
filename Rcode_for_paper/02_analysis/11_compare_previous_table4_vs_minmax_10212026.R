# ==============================================================================
# Script: 11_compare_previous_table4_vs_minmax_10212026.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Compare previous Goal 2 (Table 4) baseline SEM indicators against 
#          the new 5-Guild Min-Max Gating Attenuation models.
# Output: Rcode_for_paper/Routput_for_paper/tables/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(gt)
  library(stringr)
  library(webshot2)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & LOAD DATASETS
# ------------------------------------------------------------------------------
years_included <- "19yr"
proj_dir         <- file.path(getwd(), "Rcode_for_paper")
master_out       <- file.path(proj_dir, "Routput_for_paper")
tbl_out_dir      <- file.path(proj_dir, "tables")
data_out_dir     <- file.path(proj_dir, "data")
sub_data_out_dir <- file.path(master_out, "data")

dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load Previous Table 4 Baseline Diagnostics
t4_summary_path <- file.path(tbl_out_dir, "Table4_Indicator_Diagnostics_dAIC3.csv")
t4_eval_path    <- file.path(data_out_dir, "Goal2_Indicator_Evaluations_dAIC3.csv")

if (!file.exists(t4_summary_path) && !file.exists(t4_eval_path)) {
  # Fallback dummy creation if previous summary table file is missing
  t4_summary <- tibble(
    SEM_Node  = c("PredNCC", "PredAK"),
    Indicator = c("X10_harbour_seal_DFO_Coastwide", "X15_DFA_sleeperSharkBSAI_predAK"),
    Beta_Supp = c(NA_real_, -0.918),
    Beta_Flip = c(0.45, NA_real_),
    N_Models  = c(3, 13)
  )
} else {
  t4_summary <- if (file.exists(t4_summary_path)) read_csv(t4_summary_path, show_col_types = FALSE) else read_csv(t4_eval_path, show_col_types = FALSE)
}

# 2. Load New Min-Max Sweep Results
sweep_file <- file.path(sub_data_out_dir, paste0("Goal3_AK_Adult_Full_4Way_Competition_Sweep_", years_included, ".csv"))

if (!file.exists(sweep_file)) {
  stop(sprintf("Sweep results file not found at: %s", sweep_file))
}

cat(sprintf("Loading Min-Max Attenuation Sweep results from: %s\n", basename(sweep_file)))
minmax_results <- read_csv(sweep_file, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# STEP 1: CLEAN AND HARMONIZE PREDATOR INDICATOR NAMES FOR JOINING
# ------------------------------------------------------------------------------

# Helper function to extract base indicator name (removing _smoltyr / _adultyr)
clean_indicator_name <- function(x) {
  x <- gsub("_smoltyr$|_adultyr$", "", x)
  x <- gsub("_2yrLead_WS$", "", x)
  return(x)
}

# Clean Table 4 baseline records
t4_predators <- t4_summary %>%
  filter(SEM_Node %in% c("PredNCC", "PredAK") | grepl("^X10|^X15", Indicator)) %>%
  mutate(
    pred_clean = clean_indicator_name(Indicator),
    Baseline_Beta = if_else(!is.na(Beta_Supp), Beta_Supp, Beta_Flip),
    Baseline_Status = if_else(!is.na(Beta_Supp), "Supported (-)", "Flipped (+)")
  )

# Match Min-Max results back to Table 4 baseline indicators
minmax_matched <- minmax_results %>%
  mutate(
    pred_clean = clean_indicator_name(predator)
  )

# ------------------------------------------------------------------------------
# STEP 2: BUILD SIDE-BY-SIDE COMPARISON TABLE
# ------------------------------------------------------------------------------
comparison_table <- t4_predators %>%
  inner_join(minmax_matched, by = "pred_clean", relationship = "many-to-many") %>%
  select(
    SEM_Node,
    Indicator,
    predator_full     = predator,
    Baseline_Status,
    Baseline_Beta,
    Variant_Type      = abs_variant_type,
    Selected_Prey     = abs_prey,
    MinMax_Beta_Peff  = abs_beta_peff,
    MinMax_AIC        = abs_final_aic,
    Base_AIC          = base_aic,
    Delta_AIC         = abs_delta_aic
  ) %>%
  arrange(SEM_Node, MinMax_AIC)

out_csv <- file.path(tbl_out_dir, paste0("Table4_vs_MinMax_Predator_Comparison_", years_included, ".csv"))
write_csv(comparison_table, out_csv)

# ------------------------------------------------------------------------------
# STEP 3: GT FORMATTED TABLE GENERATION
# ------------------------------------------------------------------------------
gt_comp_df <- comparison_table %>%
  select(
    Indicator,
    Baseline_Status,
    Baseline_Beta,
    Variant_Type,
    Selected_Prey,
    MinMax_Beta_Peff,
    MinMax_AIC,
    Delta_AIC 
  )

gt_comp <- gt(gt_comp_df) %>%
  tab_header(
    title    = md(sprintf("**%s Comparative Diagnostic: Additive Table 4 SEMs vs. Min-Max Gated Models**", years_included)),
    subtitle = "Evaluating whether temperature gating and alternative prey buffering resolve baseline flipped predator signs."
  ) %>%
  tab_spanner(
    label   = md("**Previous Baseline (Table 4 Additive)**"),
    columns = c(Baseline_Status, Baseline_Beta)
  ) %>%
  tab_spanner(
    label   = md("**New Min-Max Attenuation Model**"),
    columns = c(Variant_Type, Selected_Prey, MinMax_Beta_Peff, MinMax_AIC, Delta_AIC)
  ) %>%
  cols_label(
    Indicator        = md("**Predator Indicator**"),
    Baseline_Status  = md("**Status**"),
    Baseline_Beta    = md("**Base β**"),
    Variant_Type     = md("**Variant Type**"),
    Selected_Prey    = md("**Selected Prey**"),
    MinMax_Beta_Peff = md("**Gated β**"),
    MinMax_AIC       = md("**AIC**"),
    Delta_AIC        = md("**ΔAIC**")
  ) %>%
  fmt_number(
    columns  = where(is.numeric),
    decimals = 2
  ) %>%
  cols_align(align = "left", columns = c(Indicator, Variant_Type, Selected_Prey)) %>%
  cols_align(align = "center", columns = c(Baseline_Status, Baseline_Beta, MinMax_Beta_Peff, MinMax_AIC, Delta_AIC)) %>%
  sub_missing(columns = everything(), missing_text = "-") %>%
  tab_style(
    style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#2E7D32", weight = "bold")),
    locations = cells_body(
      columns = c(MinMax_Beta_Peff, Selected_Prey),
      rows    = Baseline_Status == "Flipped (+)" & MinMax_Beta_Peff < 0
    )
  ) %>%
  tab_options(
    table.font.size           = px(10),
    heading.title.font.size   = px(12),
    column_labels.font.weight = "bold",
    data_row.padding          = px(3)
  )

print(gt_comp)

# Save Outputs (HTML and PNG)
html_path <- file.path(tbl_out_dir, paste0("Table4_vs_MinMax_Predator_Comparison_", years_included, ".html"))
png_path  <- file.path(tbl_out_dir, paste0("Table4_vs_MinMax_Predator_Comparison_", years_included, ".png"))

gtsave(gt_comp, html_path)

gtsave(
  data     = gt_comp,
  filename = png_path,
  vwidth   = 1200,
  vheight  = 800,
  zoom     = 2
)

message("Comparison complete! Table saved to: ", html_path)

