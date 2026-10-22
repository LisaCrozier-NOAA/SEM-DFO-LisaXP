# ==============================================================================
# Script: 10f_plot_top_ak_adult_models.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Generate visual diagnostics and time-series plots for top-ranking
#          AK & Adult-Stage SAR models from the competition sweep.
# Output: Rcode_for_paper/Routput_for_paper/plots/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP & DATA LOADING
# ------------------------------------------------------------------------------
years_included <- "24yr"
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
plot_out_dir <- file.path(master_out, "plots")
dir.create(plot_out_dir, showWarnings = FALSE, recursive = TRUE)

raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df       <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

valid_cols  <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
filtered_df <- raw_df %>% select(all_of(valid_cols)) %>% filter(complete.cases(.))

scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

minmax_base <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

target_cpue <- "X07_cpue_IntSprJunHW_smoltyr"
target_sar  <- "X16_DFA_sar_UCsprCh_smoltyr"
hake_col    <- "X09_DFA_HakeAge5Plus_smoltyr"

sst_smolt   <- "X21_sst_egoa_junjulaug_smoltyr"
sst_adult   <- "X21_sst_egoa_junjulaug_adultyr"

# ------------------------------------------------------------------------------
# STEP 1: DEFINE TOP 6 MODELS FOR GRAPHING
# ------------------------------------------------------------------------------
top_models_config <- list(
  list(
    pred = "X10_harbour_seal_DFO_Coastwide_smoltyr", 
    type = "Joint SST + 1 AltPrey", 
    env  = sst_smolt, 
    prey = c("X05_DFA_abundSardine_smoltyr")
  ),
  list(
    pred = "X15_DFA_sleeperSharkBSAI_predAK_smoltyr", 
    type = "Joint SST + 1 AltPrey", 
    env  = sst_smolt, 
    prey = c("capelin_WGoA_smoltyr")
  ),
  list(
    pred = "X10_Harbour_s_WS_adultyr", 
    type = "Raw Baseline", 
    env  = NULL, 
    prey = NULL
  ),
  list(
    pred = "X15_sablefishAK_biomass4plus_smoltyr", 
    type = "Joint SST + 1 AltPrey", 
    env  = sst_smolt, 
    prey = c("sitkaHerring_EGoA_smoltyr")
  ),
  list(
    pred = "X15_sablefishAK_biomass3plus_smoltyr", 
    type = "Joint SST + 2 AltPrey", 
    env  = sst_smolt, 
    prey = c("X13_pollockGOA_age1and2_smoltyr", "sitkaHerring_EGoA_smoltyr")
  ),
  list(
    pred = "X15_DFA_halibutBiomassAge8plus_predAK_adultyr", 
    type = "Joint SST + 1 AltPrey", 
    env  = sst_adult, 
    prey = c("capelin_WGoA_adultyr")
  )
)

# ------------------------------------------------------------------------------
# STEP 2: GENERATE TIME-SERIES FIT PLOTS
# ------------------------------------------------------------------------------
plot_list <- list()

for (i in seq_along(top_models_config)) {
  cfg  <- top_models_config[[i]]
  pred <- cfg$pred
  
  if (!pred %in% names(minmax_base)) {
    warning(sprintf("Predator '%s' not found in dataset. Skipping.", pred))
    next
  }
  
  P_norm <- minmax_base[[pred]]
  
  # Compute Gated P_eff (Fixing w = 1.0)
  if (grepl("Joint", cfg$type) && !is.null(cfg$env) && cfg$env %in% names(filtered_df)) {
    env_gate   <- scale_custom_bounds(filtered_df[[cfg$env]], a = 0.2, b = 1.0)
    
    if (!is.null(cfg$prey) && length(cfg$prey) > 0) {
      forage_raw  <- rowSums(minmax_base[, cfg$prey, drop = FALSE], na.rm = TRUE)
      forage_pool <- scale_custom_bounds(forage_raw, 0.0, 1.0)
      P_eff       <- P_norm * env_gate * (1 - 1.0 * forage_pool)
    } else {
      P_eff       <- P_norm * env_gate
    }
  } else if (!is.null(cfg$prey) && length(cfg$prey) > 0) {
    forage_raw  <- rowSums(minmax_base[, cfg$prey, drop = FALSE], na.rm = TRUE)
    forage_pool <- scale_custom_bounds(forage_raw, 0.0, 1.0)
    P_eff       <- P_norm * (1 - 1.0 * forage_pool)
  } else {
    P_eff       <- P_norm
  }
  
  df_plot       <- minmax_base
  df_plot$P_eff <- P_eff
  df_plot$year  <- filtered_df$year
  
  # 1. Fit Raw Baseline SEM
  base_sem_str <- paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + ", pred)
  base_fit     <- sem(base_sem_str, data = df_plot, std.lv = TRUE, missing = "ML", warn = FALSE)
  base_aic     <- round(as.numeric(fitMeasures(base_fit, "aic")), 2)
  base_r2      <- round(as.numeric(inspect(base_fit, "r2")[[target_sar]]), 3)
  
  # 2. Fit Gated SEM
  gated_sem_str <- paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_eff")
  gated_fit     <- sem(gated_sem_str, data = df_plot, std.lv = TRUE, missing = "ML", warn = FALSE)
  gated_aic     <- round(as.numeric(fitMeasures(gated_fit, "aic")), 2)
  gated_r2      <- round(as.numeric(inspect(gated_fit, "r2")[[target_sar]]), 3)
  
  # Calculate Fitted SAR manually from lavaan parameter estimates
  pe     <- parameterEstimates(gated_fit)
  b_cpue <- pe %>% filter(lhs == target_sar & op == "~" & rhs == target_cpue) %>% pull(est)
  b_peff <- pe %>% filter(lhs == target_sar & op == "~" & rhs == "P_eff") %>% pull(est)
  
  df_plot$SAR_fitted <- b_cpue * df_plot[[target_cpue]] + b_peff * df_plot$P_eff
  
  # Clean Title
  clean_title <- gsub("^X[0-9]{2}_|^DFA_", "", pred)
  clean_title <- gsub("_predAK|_smoltyr|_adultyr", "", clean_title)
  
  p <- ggplot(df_plot, aes(x = year)) +
    geom_line(aes(y = .data[[target_sar]], color = "Observed SAR"), size = 1.1) +
    geom_line(aes(y = SAR_fitted, color = "Model Fitted SAR"), linetype = "dashed", size = 1.0) +
    geom_col(aes(y = P_eff * max(.data[[target_sar]], na.rm = TRUE), fill = "Gated Pred Load (Peff)"), alpha = 0.25) +
    scale_color_manual(values = c("Observed SAR" = "black", "Model Fitted SAR" = "firebrick")) +
    scale_fill_manual(values = c("Gated Pred Load (Peff)" = "blue")) +
    labs(
      title    = paste0(i, ". ", clean_title),
      subtitle = paste0(cfg$type, " | AIC: ", base_aic, " -> ", gated_aic, " | R2: ", base_r2, " -> ", gated_r2),
      x = "Year", y = "Standardized SAR"
    ) +
    theme_minimal() +
    theme(
      plot.title      = element_text(size = 11, face = "bold"),
      plot.subtitle   = element_text(size = 9, color = "grey30"),
      legend.position = "bottom",
      legend.title    = element_blank()
    )
  
  plot_list[[i]] <- p
}

# Combine into 3x2 grid
combined_plot <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(
    title    = "Top 6 AK & Adult-Stage Predator SEM Models (SAR Dynamics)",
    subtitle = "Comparing Observed vs. Fitted SAR under Temperature Gating & Alternative Prey Buffering",
    theme    = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(combined_plot)

plot_file <- file.path(plot_out_dir, paste0("Top6_AK_Adult_SEM_Diagnostic_Plots_", years_included, ".png"))
ggsave(plot_file, combined_plot, width = 14, height = 12, dpi = 300)

cat(sprintf("\nDiagnostic time-series plots successfully written to: %s\n", plot_file))