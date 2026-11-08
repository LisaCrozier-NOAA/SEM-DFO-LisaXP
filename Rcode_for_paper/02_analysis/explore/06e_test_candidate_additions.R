

#RESULTS-------
#pink salmon resolved the heywood cases and improved the model
# --- Baseline_3DFAs ---
#   AIC: 243.89 | CFI: 0.521 | SRMR: 0.284 | SAR R2: 0.658 | SAR ~ AKPrey: beta = -0.315 (p = 0.2792)
#   Candidate Loading: N/A | Heywood Case: YES
#   
#   --- + Decapods ---
#     AIC: 305.46 | CFI: 0.484 | SRMR: 0.308 | SAR R2: 0.644 | SAR ~ AKPrey: beta = -0.358 (p = 0.0420)
#   Candidate Loading: -0.542 (p=0.002) | Heywood Case: YES
#   
#   --- + Mysids ---
#     AIC: 316.93 | CFI: 0.436 | SRMR: 0.272 | SAR R2: 0.647 | SAR ~ AKPrey: beta = -0.301 (p = 0.3226)
#   Candidate Loading: 0.045 (p=0.826) | Heywood Case: YES
#   
#   --- + PinkSalmon ---
#     AIC: 312.10 | CFI: 0.495 | SRMR: 0.278 | SAR R2: 0.732 | SAR ~ AKPrey: beta = -0.396 (p = 0.0796)
#   Candidate Loading: 0.454 (p=0.028) | Heywood Case: NO

# Quick check of Pink Salmon model loadings and covariance
syntax_pink <- "
  Latent_AKPrey =~ WGOA_lower + WGOA_mid + EGOA_seabirds + X14_pinkSalmonNorthAmerica
  SAR ~ Growth + CPUE + Latent_AKPrey + P_ocean_eff
"
fit_pink <- sem(syntax_pink, data = df_fit, std.lv = TRUE, missing = "ML", bounds = "pos.var")
parameterEstimates(fit_pink, standardized = TRUE) %>%
  filter(op == "=~") %>%
  select(rhs, std.all, z, pvalue)

# ==============================================================================
# Script: 02_analysis/06e_test_candidate_additions.R
# Purpose: Sequentially test individual candidate indicators:
#          - X12_biomassDecaShelfSum
#          - X12_biomassMysidShelfSum
#          - X14_pinkSalmonNorthAmerica
#          added to the AK latent block within DAG 1B.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# 0. SETUP & LOAD DATA
# ------------------------------------------------------------------------------
proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")

alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")

df_raw  <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)
peff_df <- read_csv(peff_ts_path, show_col_types = FALSE)

yr_col <- grep("^year$", names(df_raw), value = TRUE, ignore.case = TRUE)
df_raw$year <- if (length(yr_col) > 0) as.numeric(df_raw[[yr_col[1]]]) else as.numeric(rownames(df_raw))

df_24 <- df_raw %>% filter(year >= 1998 & year <= 2021)

find_col <- function(pattern, df) {
  matches <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
  if (length(matches) == 0) stop("Could not find column matching pattern: ", pattern)
  matches[1]
}

# Assemble fixed anchors
df_sem <- tibble(
  year        = df_24$year,
  Growth      = df_24[[find_col("IGF_mu", df_24)]],
  CPUE        = df_24[[find_col("cpue_IntSprJunHW", df_24)]],
  SAR         = df_24[[find_col("sar_UCsprCh", df_24)]],
  P_ocean_eff = peff_df$Peff_sablefishAK_biomass3plus_smoltyr,
  # Baseline 3 retained DFAs
  WGOA_lower    = df_24[[find_col("WGOA.*lower", df_24)]],
  WGOA_mid      = df_24[[find_col("WGOA.*mid", df_24)]],
  EGOA_seabirds = df_24[[find_col("EGOA.*seabird|EGOA.*sb", df_24)]]
)

# Attach the 3 new candidates
new_candidates <- c(
  "X12_biomassDecaShelfSum",
  "X12_biomassMysidShelfSum",
  "X14_pinkSalmonNorthAmerica"
)

for (cand in new_candidates) {
  matched_name <- find_col(cand, df_24)
  df_sem[[cand]] <- df_24[[matched_name]]
}

# Standardize (z-score) all numeric variables
df_fit <- df_sem %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

# ------------------------------------------------------------------------------
# 1. EVALUATION ENGINE
# ------------------------------------------------------------------------------
test_results <- list()

# Baseline 3-indicator block
base_indicators <- c("WGOA_lower", "WGOA_mid", "EGOA_seabirds")

# Run combinations: baseline alone, then baseline + each single candidate
models_to_test <- list(
  "Baseline_3DFAs" = base_indicators,
  "+ Decapods"     = c(base_indicators, "X12_biomassDecaShelfSum"),
  "+ Mysids"       = c(base_indicators, "X12_biomassMysidShelfSum"),
  "+ PinkSalmon"   = c(base_indicators, "X14_pinkSalmonNorthAmerica")
)

cat("\n====================================================================\n")
cat("TESTING CANDIDATE ADDITIONS TO LATENT AKPREY\n")
cat("====================================================================\n")

for (mod_name in names(models_to_test)) {
  indicators <- models_to_test[[mod_name]]
  
  meas_block   <- paste0("Latent_AKPrey =~ ", paste(indicators, collapse = " + "))
  struct_block <- "SAR ~ Growth + CPUE + Latent_AKPrey + P_ocean_eff"
  syntax       <- paste(meas_block, struct_block, sep = "\n")
  
  fit <- tryCatch(
    sem(
      syntax,
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
    cat(sprintf("\n[-] Model '%s': Failed to converge.\n", mod_name))
    next
  }
  
  # Check for Heywood cases (residual variances <= 0 or std.all loadings > 1.0)
  pe <- parameterEstimates(fit, standardized = TRUE)
  loadings <- pe %>% filter(op == "=~", lhs == "Latent_AKPrey")
  var_est  <- pe %>% filter(op == "~~", lhs == rhs, lhs %in% indicators)
  
  has_heywood <- any(var_est$est <= 0.001) || any(abs(loadings$std.all) > 1.0)
  
  # Structural path
  sar_path <- pe %>% filter(lhs == "SAR", rhs == "Latent_AKPrey")
  
  cand_loading <- if (mod_name != "Baseline_3DFAs") {
    cand_name <- tail(indicators, 1)
    cand_row  <- loadings %>% filter(rhs == cand_name)
    sprintf("%.3f (p=%.3f)", cand_row$std.all[1], cand_row$pvalue[1])
  } else {
    "N/A"
  }
  
  test_results[[mod_name]] <- tibble(
    Model          = mod_name,
    AIC            = round(as.numeric(fitMeasures(fit, "aic")), 2),
    p_Chi2         = round(as.numeric(fitMeasures(fit, "pvalue")), 4),
    CFI            = round(as.numeric(fitMeasures(fit, "cfi")), 3),
    SRMR           = round(as.numeric(fitMeasures(fit, "srmr")), 3),
    SAR_R2         = round(as.numeric(inspect(fit, "r2")["SAR"]), 3),
    Path_Beta      = round(sar_path$std.all[1], 3),
    Path_p         = round(sar_path$pvalue[1], 4),
    Candidate_Load = cand_loading,
    Heywood_Case   = ifelse(has_heywood, "YES", "NO")
  )
  
  cat(sprintf("\n--- %s ---\n", mod_name))
  cat(sprintf("AIC: %.2f | CFI: %.3f | SRMR: %.3f | SAR R2: %.3f | SAR ~ AKPrey: beta = %.3f (p = %.4f)\n",
              fitMeasures(fit, "aic"), fitMeasures(fit, "cfi"), fitMeasures(fit, "srmr"),
              inspect(fit, "r2")["SAR"], sar_path$std.all[1], sar_path$pvalue[1]))
  cat(sprintf("Candidate Loading: %s | Heywood Case: %s\n", cand_loading, ifelse(has_heywood, "YES", "NO")))
}

# ------------------------------------------------------------------------------
# 2. COMPARISON SUMMARY TABLE
# ------------------------------------------------------------------------------
summary_df <- bind_rows(test_results)
cat("\n====================================================================\n")
cat("SUMMARY OF CANDIDATE ADDITIONS\n")
cat("====================================================================\n")
print(as.data.frame(summary_df))

write_csv(summary_df, file.path(data_dir, "obj4_candidate_additions_summary.csv"))