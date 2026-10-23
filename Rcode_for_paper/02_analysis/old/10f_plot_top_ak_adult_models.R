# ==============================================================================
# Script: 10f_plot_top_ak_adult_models.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Generate visual diagnostics and time-series plots for top-ranking
#          AK & Adult-Stage SAR models from the 4-way competition sweep.
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

raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df       <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

if ("X01_habCompInd_smoltyr" %in% names(raw_df)) raw_df$X01_habCompInd_smoltyr <- -1 * raw_df$X01_habCompInd_smoltyr
if ("X01_habCompInd_adultyr" %in% names(raw_df)) raw_df$X01_habCompInd_adultyr <- -1 * raw_df$X01_habCompInd_adultyr

valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
filtered_df <- raw_df %>% select(all_of(valid_cols)) %>% filter(complete.cases(.))

scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

minmax_base <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"
hake_col    <- "X09_DFA_HakeAge5Plus_smoltyr"

# ------------------------------------------------------------------------------
# STEP 1: DEFINE TOP 6 MODELS FOR GRAPHING
# ------------------------------------------------------------------------------
top_models_config <- list(
  list(pred = "X10_Harbour_s_2yrLead_WS_adultyr", type = "Joint SST + 2 AltPrey", env = "X01_habCompInd_adultyr", prey = c("X12_DFA_biomassEuphShelfSum_adultyr", "X04_marketsquid_GAM_adultyr")),
  list(pred = "X10_Californian_s_l_2yrLead_WS_adultyr", type = "AltPrey Only (1 Prey)", env = NULL, prey = c("X09_DFA_HakeAge5Plus_adultyr")),
  list(pred = "X15_ArrowtoothFlounderBiomass_predAK_adultyr", type = "AltPrey Only (1 Prey)", env = NULL, prey = c("X13_capelin_WGoA_adultyr")),
  list(pred = "X15_sablefishBiomass_predAK_smoltyr", type = "Joint SST + 2 AltPrey", env = "X21_sst_egoa_junjulaug_smoltyr", prey = c("X13_pollock_age1plus_smoltyr", "X13_capelin_WGoA_smoltyr")),
  list(pred = "X15_DFA_sleeperSharks_adultyr", type = "Joint SST + 2 AltPrey", env = "X21_sst_egoa_junjulaug_adultyr", prey = c("X13_pollock_age1plus_adultyr", "X13_capelin_WGoA_adultyr")),
  list(pred = "X11_ssl_seak_pup_pred", type = "Joint SST + 2 AltPrey", env = "X21_sst_egoa_junjulaug_smoltyr", 
       prey = c("X13_sitkaHerring_EGoA_smoltyr","X13_capelin_WGoA_smoltyr"))
)


# ------------------------------------------------------------------------------
# STEP 2: GENERATE TIME-SERIES FIT PLOTS
# ------------------------------------------------------------------------------
plot_list <- list()

for (i in seq_along(top_models_config)) {
  cfg  <- top_models_config[[i]]
  pred <- cfg$pred
  
  P_norm <- minmax_base[[pred]]
  
  # Compute Gated P_eff (Fixing w = 1.0)
  if (grepl("Joint", cfg$type)) {
    env_gate <- scale_custom_bounds(filtered_df[[cfg$env]], a = 0.2, b = 1.0)
    forage_raw <- rowSums(minmax_base[, cfg$prey, drop = FALSE], na.rm = TRUE)
    forage_pool <- scale_custom_bounds(forage_raw, 0.0, 1.0)
    P_eff <- P_norm * env_gate * (1 - 1.0 * forage_pool)
  } else {
    forage_raw <- rowSums(minmax_base[, cfg$prey, drop = FALSE], na.rm = TRUE)
    forage_pool <- scale_custom_bounds(forage_raw, 0.0, 1.0)
    P_eff <- P_norm * (1 - 1.0 * forage_pool)
  }
  
  df_plot <- minmax_base
  df_plot$P_eff <- P_eff
  df_plot$Year  <- filtered_df$Year
  
  # 1. Fit Raw Baseline SEM
  base_sem_str <- paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + ", pred)
  base_fit <- sem(base_sem_str, data = df_plot, std.lv = TRUE, missing = "ML", warn = FALSE)
  base_aic <- round(as.numeric(fitMeasures(base_fit, "aic")), 2)
  base_r2  <- round(as.numeric(inspect(base_fit, "r2")[[target_sar]]), 3)
  
  # 2. Fit Gated SEM
  gated_sem_str <- paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + P_eff")
  gated_fit <- sem(gated_sem_str, data = df_plot, std.lv = TRUE, missing = "ML", warn = FALSE)
  gated_aic <- round(as.numeric(fitMeasures(gated_fit, "aic")), 2)
  gated_r2  <- round(as.numeric(inspect(gated_fit, "r2")[[target_sar]]), 3)
  
  # Calculate Fitted SAR manually from lavaan parameter estimates
  pe <- parameterEstimates(gated_fit)
  b_cpue <- pe %>% filter(lhs == target_sar & op == "~" & rhs == target_cpue) %>% pull(est)
  b_peff <- pe %>% filter(lhs == target_sar & op == "~" & rhs == "P_eff") %>% pull(est)
  
  df_plot$SAR_fitted <- b_cpue * df_plot[[target_cpue]] + b_peff * df_plot$P_eff
  
  # Clean Title
  clean_title <- gsub("^X10_|^X15_", "", pred)
  clean_title <- gsub("_predAK|_2yrLead_WS", "", clean_title)
  
  p <- ggplot(df_plot, aes(x = Year)) +
    geom_line(aes(y = X16_SAR, color = "Observed SAR"), size = 1.1) +
    geom_line(aes(y = SAR_fitted, color = "Model Fitted SAR"), linetype = "dashed", size = 1.0) +
    geom_col(aes(y = P_eff * max(X16_SAR, na.rm = TRUE), fill = "Gated Pred Load (Peff)"), alpha = 0.25) +
    scale_color_manual(values = c("Observed SAR" = "black", "Model Fitted SAR" = "firebrick")) +
    scale_fill_manual(values = c("Gated Pred Load (Peff)" = "blue")) +
    labs(
      title = paste0(i, ". ", clean_title),
      subtitle = paste0(cfg$type, " | AIC: ", base_aic, " -> ", gated_aic, " | R2: ", base_r2, " -> ", gated_r2),
      x = "Year", y = "Standardized SAR"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  plot_list[[i]] <- p
}

# Combine into 3x2 grid
combined_plot <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(
    title = "Top 6 AK & Adult-Stage Predator SEM Models (SAR Dynamics)",
    subtitle = "Comparing Observed vs. Fitted SAR under Temperature Gating & Alternative Prey Buffering",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(combined_plot)

plot_file <- file.path(plot_out_dir, paste0("Top6_AK_Adult_SEM_Diagnostic_Plots_", years_included, ".png"))
ggsave(plot_file, combined_plot, width = 14, height = 12, dpi = 300)

cat(sprintf("\nDiagnostic time-series plots successfully written to: %s\n", plot_file))