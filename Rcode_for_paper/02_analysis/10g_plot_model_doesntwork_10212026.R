# ==============================================================================
# Script: 11_plot_top_sem_models.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Dynamically generate and export path diagrams for top SEM models 
#          featuring dynamic SST gating and alternative prey predictors.
# Output: Rcode_for_paper/Routput_for_paper/figures/sem_top_models/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(semPlot)
})

# ------------------------------------------------------------------------------
# STEP 0: PATHS & CONFIGURATION
# ------------------------------------------------------------------------------
years_included <- "24yr" # "19yr" or "24yr"

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
fig_out_dir  <- file.path(master_out, "figures", "sem_top_models")

if (!dir.exists(fig_out_dir)) dir.create(fig_out_dir, recursive = TRUE)

# Load MinMax Complete Data
data_path <- file.path(data_out_dir, paste0("minmax_complete_data_10212026_", years_included, ".csv"))
df_base   <- read_csv(data_path, show_col_types = FALSE)

# Target CPUE response variable
target_cpue <- names(df_base)[grepl("^X07.*cpue_IntSprJunHW", names(df_base))][1]

# Define Top Models Input Table
top_models <- tribble(
  ~predator, ~abs_variant_type, ~abs_prey,
  "X10_harbour_seal_DFO_Coastwide_smoltyr", "Joint SST + 1 AltPrey", "X05_DFA_abundSardine_smoltyr",
  "X15_DFA_sleeperSharkBSAI_predAK_smoltyr", "Joint SST + 1 AltPrey", "capelin_WGoA_smoltyr",
  "X10_Harbour_s_WS_adultyr", "Raw Baseline", "None",
  "X15_sablefishAK_biomass4plus_smoltyr", "Joint SST + 1 AltPrey", "sitkaHerring_EGoA_smoltyr",
  "X15_sablefishAK_biomass3plus_smoltyr", "Joint SST + 2 AltPrey", "X13_pollockGOA_age1and2_smoltyr + sitkaHerring_EGoA_smoltyr",
  "X15_DFA_halibutBiomassAge8plus_predAK_adultyr", "Joint SST + 1 AltPrey", "capelin_WGoA_adultyr",
  "X11_ssl_seak_pup_pred_smoltyr", "Joint SST + 2 AltPrey", "sitkaHerring_EGoA_smoltyr + capelin_WGoA_smoltyr",
  "X15_salmonSharkGoA_predAK_smoltyr", "Joint SST + 2 AltPrey", "sitkaHerring_EGoA_smoltyr + capelin_WGoA_smoltyr"
)

# Custom bound scaling helper
scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

# ------------------------------------------------------------------------------
# STEP 1: DYNAMIC MODEL PARSING & PLOTTING LOOP
# ------------------------------------------------------------------------------
cat(sprintf("Generating path diagrams for %d top SEM models...\n", nrow(top_models)))

for (i in seq_len(nrow(top_models))) {
  row          <- top_models[i, ]
  pred_col     <- row$predator
  variant_type <- row$abs_variant_type
  raw_prey_str <- row$abs_prey
  
  # Determine life stage suffix (adultyr vs smoltyr)
  stage_sfx <- if (grepl("_adultyr$", pred_col)) "adultyr" else "smoltyr"
  sst_col   <- paste0("X21_sst_egoa_junjulaug_", stage_sfx)
  
  # Verify required columns exist in dataset
  if (!pred_col %in% names(df_base)) {
    warning(sprintf("Predator column '%s' not found in dataset. Skipping.", pred_col))
    next
  }
  
  df_model <- df_base
  
  # Construct Predictor Terms
  rhs_terms <- c()
  
  # 1. Handle SST Gating / SST Main Effects
  if (grepl("SST", variant_type, ignore.case = TRUE)) {
    if (sst_col %in% names(df_model)) {
      # Calculate effective gated predator (P_eff) using standard default parameters
      P_norm <- df_model[[pred_col]]
      S_norm <- scale_custom_bounds(df_model[[sst_col]], a = 0.0, b = 1.0)
      
      df_model$P_eff <- P_norm * (1 - 0.5 * (1 - S_norm))
      rhs_terms <- c(rhs_terms, "P_eff")
    } else {
      # Fallback to base predator if SST column missing
      rhs_terms <- c(rhs_terms, pred_col)
    }
  } else {
    rhs_terms <- c(rhs_terms, pred_col)
  }
  
  # 2. Parse AltPrey Terms
  if (raw_prey_str != "None" && !is.na(raw_prey_str)) {
    prey_vars <- unlist(strsplit(raw_prey_str, "\\s*\\+\\s*"))
    valid_prey <- intersect(prey_vars, names(df_model))
    if (length(valid_prey) > 0) {
      rhs_terms <- c(rhs_terms, valid_prey)
    }
  }
  
  # Build Lavaan Syntax
  model_syntax <- paste0(target_cpue, " ~ ", paste(rhs_terms, collapse = " + "))
  
  # Fit SEM Model
  fit <- tryCatch(
    sem(model_syntax, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !lavInspect(fit, "converged")) {
    warning(sprintf("Model failed to converge for predator '%s'. Skipping plot.", pred_col))
    next
  }
  
  # ----------------------------------------------------------------------------
  # STEP 2: GENERATE AND EXPORT PATH DIAGRAM
  # ----------------------------------------------------------------------------
  clean_pred_name <- gsub("^(X[0-9]{2}_|DFA_)", "", pred_col)
  fig_filename    <- file.path(fig_out_dir, paste0("rank_", i, "_", clean_pred_name, ".png"))
  
  png(filename = fig_filename, width = 2400, height = 1800, res = 300)
  
  semPaths(
    fit,
    what             = "std",
    layout           = "tree",
    edge.label.cex   = 1.1,
    curve            = 2,
    style            = "lisrel",
    nCharNodes       = 0,
    mar              = c(6, 8, 6, 8),
    intercepts       = FALSE,
    residuals        = FALSE,
    edge.color       = "black",
    fade             = FALSE,
    title            = FALSE # Set to FALSE so semPaths doesn't reserve default title space
  )
  
  # 2. Add custom title using base R's title() function
  title(
    main = paste0("Rank ", i, ": ", clean_pred_name, "\n[", variant_type, "]"),
    line = 2,
    cex.main = 1.2
  ) 
  dev.off()
  cat(sprintf("[%d/%d] Saved path diagram: %s\n", i, nrow(top_models), fig_filename))
}

cat(sprintf("\nAll top SEM model path diagrams generated successfully in:\n%s\n", fig_out_dir))
