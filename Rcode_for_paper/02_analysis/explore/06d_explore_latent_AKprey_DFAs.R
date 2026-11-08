

# ==============================================================================
# Script: 02_analysis/06d_explore_latent_AKprey_DFAs.R
# Purpose: Exploratory single-node latent SEM for Alaska Prey using the intact
#          DFA trends that met the Figure 1 threshold (wi >= 0.05).
#          Backward prunes non-significant factor loadings (=~) one at a time.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# 0. DIRECTORIES & PATHS
# ------------------------------------------------------------------------------
proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")

alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")

if (!file.exists(alt_data_path)) stop("Missing: ", alt_data_path)
if (!file.exists(peff_ts_path))  stop("Missing: ", peff_ts_path)

# ------------------------------------------------------------------------------
# 1. ASSEMBLE 24-YEAR DATASET WITH INTACT DFA TRENDS
# ------------------------------------------------------------------------------
df_raw  <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)
peff_df <- read_csv(peff_ts_path, show_col_types = FALSE)

yr_col <- grep("^year$", names(df_raw), value = TRUE, ignore.case = TRUE)
df_raw$year <- if (length(yr_col) > 0) as.numeric(df_raw[[yr_col[1]]]) else as.numeric(rownames(df_raw))

# Filter to 1998-2021
df_24 <- df_raw %>% filter(year >= 1998 & year <= 2021)

# Helper: locate column flexibly
find_col <- function(pattern, df) {
  matches <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
  if (length(matches) == 0) stop("Could not find column matching: ", pattern)
  matches[1]
}

# Core DAG 1B anchors
df_sem <- tibble(
  year        = df_24$year,
  Growth      = df_24[[find_col("IGF_mu", df_24)]],
  CPUE        = df_24[[find_col("cpue_IntSprJunHW", df_24)]],
  SAR         = df_24[[find_col("sar_UCsprCh", df_24)]],
  P_ocean_eff = peff_df$Peff_sablefishAK_biomass3plus_smoltyr
)

# Attach the 4 Figure 1 intact DFA trends
# Adjust regex patterns if column headers in sem_altprey_data differ slightly:
df_sem$EGOA_lower    <- df_24[[find_col("EGOA.*lower", df_24)]]
df_sem$WGOA_lower    <- df_24[[find_col("WGOA.*lower", df_24)]]
df_sem$WGOA_mid      <- df_24[[find_col("WGOA.*mid", df_24)]]
df_sem$EGOA_seabirds <- df_24[[find_col("EGOA.*seabird|EGOA.*sb", df_24)]]

# Z-score all modeling variables
df_fit <- df_sem %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

cat("Loaded 24-year complete dataset with 4 intact AK prey DFAs:\n")
cat("Indicators: EGOA_lower, WGOA_lower, WGOA_mid, EGOA_seabirds\n")

# Quick pairwise correlation check
cat("\nCorrelations with SAR:\n")
print(round(cor(df_fit %>% select(SAR, EGOA_lower, WGOA_lower, WGOA_mid, EGOA_seabirds)), 3))

# ------------------------------------------------------------------------------
# 2. PRUNING ENGINE FOR LATENT AK PREY
# ------------------------------------------------------------------------------
build_syntax <- function(indicators) {
  if (length(indicators) == 2) {
    # Equating loadings identifies a 2-indicator block and prevents negative variance
    meas_block <- sprintf("Latent_AKPrey =~ a*%s + a*%s", indicators[1], indicators[2])
  } else {
    meas_block <- paste0("Latent_AKPrey =~ ", paste(indicators, collapse = " + "))
  }
  
  # Streamlined DAG 1B structural path focusing directly on late-marine survival
  struct_block <- "SAR ~ Growth + CPUE + Latent_AKPrey + P_ocean_eff"
  
  paste(meas_block, struct_block, sep = "\n")
}

current_indicators <- c("EGOA_lower", "WGOA_lower", "WGOA_mid", "EGOA_seabirds")
step <- 1
step_history <- list()

cat("\n====================================================================\n")
cat("STARTING STEPWISE PRUNING: Latent_AKPrey (Figure 1 DFAs)\n")
cat("====================================================================\n")

repeat {
  cat(sprintf("\n>>> STEP %d: Testing %d Indicators <<<\n", step, length(current_indicators)))
  print(current_indicators)
  
  syntax <- build_syntax(current_indicators)
  
  fit_err <- NULL
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
    error = function(e) {
      fit_err <<- e$message
      return(NULL)
    }
  )
  
  if (is.null(fit)) {
    cat("Model failed to fit. Error message:\n", fit_err, "\n")
    break
  }
  
  if (!lavInspect(fit, "converged")) {
    cat("Optimizer failed to converge at step", step, "\n")
    break
  }
  
  # Fit measures
  aic_val  <- round(as.numeric(fitMeasures(fit, "aic")), 2)
  cfi_val  <- round(as.numeric(fitMeasures(fit, "cfi")), 3)
  srmr_val <- round(as.numeric(fitMeasures(fit, "srmr")), 3)
  p_chi2   <- round(as.numeric(fitMeasures(fit, "pvalue")), 4)
  r2_sar   <- round(as.numeric(inspect(fit, "r2")["SAR"]), 3)
  
  pe <- parameterEstimates(fit, standardized = TRUE)
  
  # Focal structural path: Latent_AKPrey -> SAR
  sar_path  <- pe %>% filter(lhs == "SAR", rhs == "Latent_AKPrey")
  path_beta <- round(sar_path$std.all[1], 3)
  path_se   <- round(sar_path$se[1], 3)
  path_p    <- round(sar_path$pvalue[1], 4)
  
  cat(sprintf("\nModel Fit : AIC = %.2f | p(Chi2) = %.4f | CFI = %.3f | SRMR = %.3f | SAR R2 = %.3f\n",
              aic_val, p_chi2, cfi_val, srmr_val, r2_sar))
  cat(sprintf("Path Beta : SAR ~ Latent_AKPrey: beta = %.3f (SE = %.3f, p = %.4f)\n", 
              path_beta, path_se, path_p))
  
  # Factor loadings
  loadings <- pe %>%
    filter(op == "=~", lhs == "Latent_AKPrey") %>%
    select(indicator = rhs, std_loading = std.all, z_val = z, pvalue) %>%
    arrange(desc(pvalue))
  
  cat("\nFactor Loadings (=~):\n")
  print(as.data.frame(loadings))
  
  # Log step
  step_history[[step]] <- tibble(
    step         = step,
    n_indicators = length(current_indicators),
    AIC          = aic_val,
    CFI          = cfi_val,
    SRMR         = srmr_val,
    p_chi2       = p_chi2,
    SAR_R2       = r2_sar,
    Latent_Beta  = path_beta,
    Latent_p     = path_p,
    worst_ind    = loadings$indicator[1],
    worst_p      = round(loadings$pvalue[1], 4)
  )
  
  if (length(current_indicators) <= 2) {
    cat("\nReached minimum 2 indicators for latent variable. Pruning stops.\n")
    break
  }
  
  worst_ind <- loadings$indicator[1]
  worst_p   <- loadings$pvalue[1]
  
  if (is.na(worst_p) || worst_p < 0.05) {
    cat(sprintf("\n[SUCCESS] All remaining factor loadings are statistically significant (p < 0.05)!\n"))
    break
  }
  
  cat(sprintf("\n>>> Backward Pruning: Dropping '%s' (loading p = %.4f)\n", worst_ind, worst_p))
  current_indicators <- setdiff(current_indicators, worst_ind)
  step <- step + 1
}

# ------------------------------------------------------------------------------
# 3. SAVE AND EXAMINE SUMMARY
# ------------------------------------------------------------------------------
prune_summary <- bind_rows(step_history)
cat("\n====================================================================\n")
cat("FINAL PRUNING TRAJECTORY SUMMARY (DFA LATENT NODE)\n")
cat("====================================================================\n")
print(as.data.frame(prune_summary))

out_path <- file.path(data_dir, "obj4_explore_AKprey_DFAs_pruning_history.csv")
write_csv(prune_summary, out_path)
cat(sprintf("\nSaved summary to: %s\n", out_path))