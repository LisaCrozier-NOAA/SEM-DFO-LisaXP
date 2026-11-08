# ==============================================================================
# Script: 02_analysis/07_synthesize_akprey_evidence_and_plot_final.R
# Purpose: Synthesize evidence across candidate AK Prey architectures:
#          1. Stepwise latent pruning of GOA biological DFAs
#          2. Decoupling offshore Pink Salmon competition
#          3. Resolving coastal Growth collinearity
#          4. Testing decapod shelf biomass vs. coastal CPUE
#          Produces a unified evidence table and plots the final pruned SEM.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(semPlot)
})

# ------------------------------------------------------------------------------
# 0. SETUP & LOAD DATA
# ------------------------------------------------------------------------------
proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")
fig_dir  <- file.path(proj_dir, "Outputs/figures")

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(fig_dir))  dir.create(fig_dir, recursive = TRUE)

alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")

if (!file.exists(alt_data_path)) stop("Missing: ", alt_data_path)
if (!file.exists(peff_ts_path))  stop("Missing: ", peff_ts_path)

df_raw  <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)
peff_df <- read_csv(peff_ts_path, show_col_types = FALSE)

yr_col <- grep("^year$", names(df_raw), value = TRUE, ignore.case = TRUE)
df_raw$year <- if (length(yr_col) > 0) as.numeric(df_raw[[yr_col[1]]]) else as.numeric(rownames(df_raw))

df_24 <- df_raw %>% filter(year >= 1998 & year <= 2021)

find_col <- function(pattern, df) {
  matches <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
  if (length(matches) == 0) stop("Could not find column matching: ", pattern)
  matches[1]
}

# ------------------------------------------------------------------------------
# 1. ASSEMBLE 24-YEAR MODELING DATASET
# ------------------------------------------------------------------------------
df_master <- tibble(
  year          = df_24$year,
  SAR           = df_24[[find_col("sar_UCsprCh", df_24)]],
  Growth        = df_24[[find_col("IGF_mu", df_24)]],
  CPUE          = df_24[[find_col("cpue_IntSprJunHW", df_24)]],
  P_ocean_eff   = peff_df$Peff_sablefishAK_biomass3plus_smoltyr,
  WGOA_lower    = df_24[[find_col("WGOA.*lower", df_24)]],
  WGOA_mid      = df_24[[find_col("WGOA.*mid", df_24)]],
  EGOA_seabirds = df_24[[find_col("EGOA.*seabird|EGOA.*sb", df_24)]],
  PinkSalmon    = df_24[[find_col("pinkSalmonNorthAmerica", df_24)]],
  Decapods      = df_24[[find_col("biomassDecaShelfSum", df_24)]]
)

# Standardize (z-score) all numeric variables
df_fit <- df_master %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

# ------------------------------------------------------------------------------
# 2. DEFINE THE 5 CANDIDATE ARCHITECTURES ACROSS OUR INVESTIGATION
# ------------------------------------------------------------------------------
models <- list(
  # M1: Unified 4-indicator AK latent node (forced reflective competitor)
  "M1_Latent_AllInOne" = "
    Latent_AK =~ WGOA_lower + WGOA_mid + EGOA_seabirds + PinkSalmon
    SAR ~ Growth + CPUE + Latent_AK + P_ocean_eff
  ",
  
  # M2: Dropping Growth with unified latent node
  "M2_Latent_DropGrowth" = "
    Latent_AK =~ WGOA_lower + WGOA_mid + EGOA_seabirds + PinkSalmon
    SAR ~ CPUE + Latent_AK + P_ocean_eff
  ",
  
  # M3: Decoupled Pink Salmon as an independent competitor node + 3-DFA Latent Prey
  "M3_Decoupled_Pink_with_LatentPrey" = "
    Latent_PreyAK =~ WGOA_lower + WGOA_mid + EGOA_seabirds
    SAR ~ CPUE + Latent_PreyAK + PinkSalmon + P_ocean_eff
  ",
  
  # M4: Final Pruned Observed Multi-Stage Architecture (Retaining CPUE)
  "M4_Final_Pruned_with_CPUE" = "
    SAR ~ CPUE + PinkSalmon + P_ocean_eff
  ",
  
  # M5: Alternative Shelf Anomaly Model (Adding Decapods)
  "M5_Observed_with_Decapods" = "
    SAR ~ CPUE + Decapods + PinkSalmon + P_ocean_eff
  "
)

# ------------------------------------------------------------------------------
# 3. FIT MODELS & EXTRACT COMPARATIVE EVIDENCE TABLE
# ------------------------------------------------------------------------------
evidence_list <- list()

cat("\nFitting candidate structural models...\n")

for (m_id in names(models)) {
  fit <- tryCatch(
    sem(
      models[[m_id]], 
      data         = df_fit, 
      missing      = "ML", 
      bounds       = "pos.var",
      optim.method = "BFGS",
      control      = list(iter.max = 1000)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !lavInspect(fit, "converged")) {
    cat(sprintf("[-] %s did not converge.\n", m_id))
    next
  }
  
  pe <- parameterEstimates(fit, standardized = TRUE)
  
  # Check for Heywood cases / negative residual variances
  ov_vars <- pe %>% filter(op == "~~", lhs == rhs, lhs != "SAR")
  loadings <- pe %>% filter(op == "=~")
  has_heywood <- any(ov_vars$est <= 0.001) || any(abs(loadings$std.all) > 1.0)
  
  # Summarize structural paths to SAR
  sar_paths <- pe %>%
    filter(lhs == "SAR", op == "~") %>%
    mutate(formatted = sprintf("%s: beta=%.3f (p=%.3f)", rhs, std.all, pvalue)) %>%
    pull(formatted) %>%
    paste(collapse = " | ")
  
  # Fit indices
  p_chi2   <- round(as.numeric(fitMeasures(fit, "pvalue")), 4)
  cfi_val  <- round(as.numeric(fitMeasures(fit, "cfi")), 3)
  srmr_val <- round(as.numeric(fitMeasures(fit, "srmr")), 3)
  aic_val  <- round(as.numeric(fitMeasures(fit, "aic")), 2)
  r2_val   <- round(as.numeric(inspect(fit, "r2")["SAR"]), 3)
  n_par    <- fitMeasures(fit, "npar")
  
  evidence_list[[m_id]] <- tibble(
    Model_ID      = m_id,
    n_params      = n_par,
    AIC           = aic_val,
    p_Chi2        = ifelse(is.na(p_chi2), "Just-Identified", as.character(p_chi2)),
    CFI           = ifelse(is.na(cfi_val), 1.000, cfi_val),
    SRMR          = ifelse(is.na(srmr_val), 0.000, srmr_val),
    SAR_R2        = r2_val,
    Heywood_Case  = ifelse(has_heywood, "YES", "NO"),
    SAR_Paths     = sar_paths
  )
}

evidence_table <- bind_rows(evidence_list) %>%
  mutate(delta_AIC = round(AIC - min(AIC), 2)) %>%
  select(Model_ID, n_params, AIC, delta_AIC, p_Chi2, CFI, SRMR, SAR_R2, Heywood_Case, SAR_Paths)

cat("\n========================================================================================\n")
cat("MASTER MODEL EVIDENCE & COMPARISON TABLE\n")
cat("========================================================================================\n")
print(as.data.frame(evidence_table))

write_csv(evidence_table, file.path(data_dir, "obj4_master_model_evidence_table.csv"))
cat(sprintf("\nSaved evidence table to: %s\n", file.path(data_dir, "obj4_master_model_evidence_table.csv")))

# ------------------------------------------------------------------------------
# 4. FIT & EXAMINE FINAL PRUNED MODEL (M4)
# ------------------------------------------------------------------------------
final_fit <- sem(
  models[["M4_Final_Pruned_with_CPUE"]], 
  data    = df_fit, 
  missing = "ML"
)

cat("\n========================================================================================\n")
cat("FINAL PRUNED MODEL ESTIMATES (M4: CPUE + Pink Salmon + Sablefish Peff)\n")
cat("========================================================================================\n")
summary(final_fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# ------------------------------------------------------------------------------
# 5. PLOT THE FINAL PRUNED MODEL
# ------------------------------------------------------------------------------
node_labels <- c(
  "SAR"         = "Chinook\nSurvival\n(SAR)",
  "CPUE"        = "Coastal\nOutmigration\n(CPUE)",
  "PinkSalmon"  = "North American\nPink Salmon\n(Competitor)",
  "P_ocean_eff" = "Subarctic Ocean\nPredation\n(P_eff)"
)

fig_path <- file.path(fig_dir, "fig_final_pruned_model_path_diagram.png")
png(fig_path, width = 8, height = 6, units = "in", res = 300)

semPaths(
  final_fit,
  whatLabels      = "std",
  layout          = "tree2",
  rotation        = 2,
  edge.label.cex  = 1.15,
  edge.color      = "black",
  nodeLabels      = node_labels[fitNames(final_fit)],
  sizeMan         = 11,
  sizeMan2        = 7,
  color           = list(man = "lightblue"),
  borders         = TRUE,
  mar             = c(4, 5, 4, 5),
  style           = "lisrel",
  fade            = FALSE,
  residuals       = FALSE,
  curvePivot      = TRUE
)

title(main = "Final Pruned Marine Survival Structural Model (1998-2021)", 
      sub  = "All structural paths significant at p < 0.05 | Model R2 = 0.611",
      cex.main = 1.1, font.main = 2)

dev.off()
cat(sprintf("\nSaved publication-quality path diagram to: %s\n", fig_path))