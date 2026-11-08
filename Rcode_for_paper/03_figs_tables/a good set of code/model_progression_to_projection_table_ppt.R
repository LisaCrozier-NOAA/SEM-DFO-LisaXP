# ==============================================================================
# Script: 03_figs_tables/model_progression_to_projection_table_ppt.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Empirically fit all progression models on the complete 24-year dataset
#          to verify exact R2, AIC, and parameter estimates, including the
#          direct Hake projection models (dropping intermediate CPUE).
# Output: Rcode_for_paper/Routput_for_paper/tables/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(gt)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & PREPARE 24-YEAR COMPLETE DATA
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
data_out_dir <- file.path(proj_dir, "Routput_for_paper", "data")
tbl_out_dir  <- file.path(proj_dir, "Routput_for_paper", "tables")

raw_csv_path  <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
raw_df_master <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df_master <- raw_df_master[, names(raw_df_master) != "" & !is.na(names(raw_df_master))]

# Enforce 24-year complete cases
valid_cols <- names(raw_df_master)[colSums(!is.na(raw_df_master)) == 24]
df_raw <- raw_df_master %>% 
  select(all_of(valid_cols)) %>% 
  filter(complete.cases(.))

scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

# Base min-max [0, 1] dataframe
df_base <- df_raw %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

target_cpue      <- "X07_cpue_IntSprJunHW_smoltyr"
target_sar       <- "X16_DFA_sar_UCsprCh_smoltyr"
hake_col         <- "X09_DFA_HakeAge5Plus_smoltyr"
inriver_pred_col <- "X10_DFA_AllSeaLionsEMB_adultyr"
sst_col          <- "X21_sst_egoa_junjulaug_smoltyr"

# ------------------------------------------------------------------------------
# STEP 1: CONSTRUCT THE TWO EFFECTIVE OCEAN PREDICTORS (a = 0.2, w = 1.0)
# ------------------------------------------------------------------------------
# Thermal gate
sst_gate <- scale_custom_bounds(df_raw[[sst_col]], a = 0.2, b = 1.0)

# A) Sleeper Shark: Capelin AltPrey
capelin_gate <- 1 - scale_custom_bounds(df_base[["capelin_WGoA_smoltyr"]], a = 0.2, b = 1.0)
df_base$P_eff_sleeper <- df_base[["X15_DFA_sleeperSharkBSAI_predAK_smoltyr"]] * sst_gate * capelin_gate

# B) Sablefish 3+: Pollock + Sitka Herring AltPrey
raw_forage_sable <- rowSums(df_base[, c("X13_pollockGOA_age1and2_smoltyr", "sitkaHerring_EGoA_smoltyr")], na.rm = TRUE)
sable_prey_gate  <- 1 - scale_custom_bounds(raw_forage_sable, a = 0.2, b = 1.0)
df_base$P_eff_sable <- df_base[["X15_sablefishAK_biomass3plus_smoltyr"]] * sst_gate * sable_prey_gate

# Raw unmodulated versions
df_base$P_raw_sleeper <- df_base[["X15_DFA_sleeperSharkBSAI_predAK_smoltyr"]]
df_base$P_raw_sable   <- df_base[["X15_sablefishAK_biomass3plus_smoltyr"]]

# ------------------------------------------------------------------------------
# STEP 2: DEFINE ALL EMPIRICAL PROGRESSION MODELS TO TEST
# ------------------------------------------------------------------------------
model_suite <- tribble(
  ~Step_ID, ~Step_Name, ~Model_Type, ~Model_Syntax,
  
  # --- Standard Chain (CPUE mediated by Hake) ---
  "1a", "Baseline (No AK Preds)", "2-Stage (CPUE)",
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue),
  
  "2a", "Raw Sleeper (No Covariates)", "2-Stage (CPUE)",
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_raw_sleeper"),
  
  "2b", "Raw Sablefish (No Covariates)", "2-Stage (CPUE)",
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_raw_sable"),
  
  "3a", "Sleeper + SST + Capelin", "2-Stage (CPUE)",
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_eff_sleeper"),
  
  "3b", "Sablefish + SST + Pollock/Herring", "2-Stage (CPUE)",
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_eff_sable"),
  
  "4a", "Final A: Sleeper + Covs + embSSL", "2-Stage (CPUE)",
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_eff_sleeper + ", inriver_pred_col),
  
  "4b", "Final B: Sablefish + Covs + embSSL", "2-Stage (CPUE)",
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_eff_sable + ", inriver_pred_col),
  
  # --- Projection Formulations (CPUE Bypassed: SAR ~ Hake + ...) ---
  "5a", "Direct Projection: Hake Only", "Projection (No CPUE)",
  paste0(target_sar, " ~ ", hake_col),
  
  "5b", "Direct Projection: Final A (Hake + Sleeper_eff + embSSL)", "Projection (No CPUE)",
  paste0(target_sar, " ~ ", hake_col, " + P_eff_sleeper + ", inriver_pred_col),
  
  "5c", "Direct Projection: Final B (Hake + Sable_eff + embSSL)", "Projection (No CPUE)",
  paste0(target_sar, " ~ ", hake_col, " + P_eff_sable + ", inriver_pred_col)
)

# ------------------------------------------------------------------------------
# STEP 3: FIT MODELS AND EXTRACT PRECISE R2 AND AIC
# ------------------------------------------------------------------------------
fitted_results <- list()

for (i in seq_len(nrow(model_suite))) {
  row <- model_suite[i, ]
  
  fit <- tryCatch(
    sem(row$Model_Syntax, data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(fit) && lavInspect(fit, "converged")) {
    r2_vals  <- inspect(fit, "r2")
    sar_r2   <- if (target_sar %in% names(r2_vals)) as.numeric(r2_vals[[target_sar]]) else NA_real_
    cpue_r2  <- if (target_cpue %in% names(r2_vals)) as.numeric(r2_vals[[target_cpue]]) else NA_real_
    curr_aic <- as.numeric(fitMeasures(fit, "aic"))
    
    # Extract standardized coefficients on SAR equation
    pe <- parameterEstimates(fit, standardized = TRUE) %>% 
      filter(op == "~" & lhs == target_sar)
    
    pe_summary <- pe %>%
      mutate(desc = sprintf("%s (β_std=%.2f, p=%.3f)", rhs, std.all, pvalue)) %>%
      pull(desc) %>%
      paste(collapse = "; ")
    
    fitted_results[[i]] <- tibble(
      Step_ID     = row$Step_ID,
      Step_Name   = row$Step_Name,
      Model_Type  = row$Model_Type,
      SAR_R2      = round(sar_r2, 3),
      CPUE_R2     = round(cpue_r2, 3),
      AIC         = round(curr_aic, 2),
      Estimates   = pe_summary
    )
  }
}

verified_table <- bind_rows(fitted_results)

# ------------------------------------------------------------------------------
# STEP 4: PRINT AUDIT OUTPUT AND BUILD SLIDE TABLE
# ------------------------------------------------------------------------------
cat("\n==============================================================================\n")
cat("      EMPIRICAL VERIFICATION OF PROGRESSION AND PROJECTION MODELS (24-YEAR)   \n")
cat("==============================================================================\n\n")

print(as.data.frame(verified_table %>% select(Step_ID, Step_Name, Model_Type, SAR_R2, AIC, Estimates)))

# Format clean presentation GT table
progression_verified_gt <- verified_table %>%
  select(Step_Name, Model_Type, SAR_R2, AIC, Estimates) %>%
  gt() %>%
  tab_header(
    title = md("**Empirical Verification of Chinook Survival Models (24-Year Series)**"),
    subtitle = "Evaluating variance explained (SAR R²) across structural progression and direct projection models"
  ) %>%
  cols_label(
    Step_Name  = md("**Model Specification**"),
    Model_Type = md("**Framework**"),
    SAR_R2     = md("**SAR R²**"),
    AIC        = md("**AIC**"),
    Estimates  = md("**Fitted Paths on SAR (Std. Estimate & p-value)**")
  ) %>%
  cols_align(align = "left", columns = c(Step_Name, Model_Type, Estimates)) %>%
  cols_align(align = "center", columns = c(SAR_R2, AIC)) %>%
  tab_options(
    table.width                = px(950),
    table.font.size            = px(12),
    heading.title.font.size    = px(15),
    heading.subtitle.font.size = px(11),
    column_labels.font.weight  = "bold",
    data_row.padding           = px(5)
  )

print(progression_verified_gt)

# Save audit files
out_csv <- file.path(data_out_dir, "Goal5_Model_Progression_to_Projection.csv")
write_csv(verified_table, out_csv)
gtsave(progression_verified_gt, file.path(tbl_out_dir, "Table_Verified_Progression_Models_24yr.html"))

cat(sprintf("\nVerification data saved to: %s\n", out_csv))


gtsave(
  data     = progression_verified_gt,
  filename = file.path(tbl_out_dir, "Goal5_Table_Verified_Progression_Models.png"),
  vwidth   = 1000,
  vheight  = 600
)
