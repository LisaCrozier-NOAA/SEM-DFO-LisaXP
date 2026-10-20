# ==============================================================================
# Script: plot_top4_joint_sem_and_ingredients.R
# Purpose: Hardcode SEM models for top 4 Joint (SST + 2 AltPrey) predators,
#          fit lavaan models, and plot raw ingredient time series with Peff on top.
# Output: Rcode_for_paper/Routput_for_paper/plots/Top4_Joint_Models_Ingredients.png
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP & DATA PREPARATION
# ------------------------------------------------------------------------------
years_included <- "24yr"
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
plot_out_dir <- file.path(master_out, "plots")
dir.create(plot_out_dir, showWarnings = FALSE, recursive = TRUE)

raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df       <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

# Invert HCI so High = Warm / Compressed
if ("X01_habCompInd_smoltyr" %in% names(raw_df)) raw_df$X01_habCompInd_smoltyr <- -1 * raw_df$X01_habCompInd_smoltyr
if ("X01_habCompInd_adultyr" %in% names(raw_df)) raw_df$X01_habCompInd_adultyr <- -1 * raw_df$X01_habCompInd_adultyr

valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
filtered_df <- raw_df %>% select(all_of(valid_cols)) %>% filter(complete.cases(.))

scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

minmax_df <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

target_cpue <- "X07_DFA_cpue_IntSprJunHW"
target_sar  <- "X16_SAR"
hake_col    <- "X09_DFA_HakeAge5Plus_smoltyr"

# ------------------------------------------------------------------------------
# STEP 1: HARDCODE THE 4 JOINT MODELS
# ------------------------------------------------------------------------------

# --- MODEL 1: Harbor Seals Adultyr ---
# P_norm = X10_Harbour_s_2yrLead_WS_adultyr
# Env    = X01_habCompInd_adultyr
# Prey   = X12_DFA_biomassEuphShelfSum_adultyr + X04_marketsquid_GAM_adultyr
env_gate_1   <- scale_custom_bounds(filtered_df$X01_habCompInd_adultyr, a = 0.2, b = 1.0)
forage_raw_1 <- rowSums(minmax_df[, c("X12_DFA_biomassEuphShelfSum_adultyr", "X04_marketsquid_GAM_adultyr")], na.rm = TRUE)
forage_pool_1<- scale_custom_bounds(forage_raw_1, 0.0, 1.0)

minmax_df$Peff_M1 <- minmax_df$X10_Harbour_s_2yrLead_WS_adultyr * env_gate_1 * (1 - 1.0 * forage_pool_1)

m1_sem <- sem(
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + Peff_M1"),
  data = minmax_df, std.lv = TRUE, missing = "ML", warn = FALSE
)

# --- MODEL 2: Sablefish Biomass Smoltyr ---
# P_norm = X15_sablefishBiomass_predAK_smoltyr
# Env    = X21_sst_egoa_junjulaug_smoltyr
# Prey   = X13_pollock_age1plus_smoltyr + X13_capelin_WGoA_smoltyr
env_gate_2   <- scale_custom_bounds(filtered_df$X21_sst_egoa_junjulaug_smoltyr, a = 0.2, b = 1.0)
forage_raw_2 <- rowSums(minmax_df[, c("X13_pollock_age1plus_smoltyr", "X13_capelin_WGoA_smoltyr")], na.rm = TRUE)
forage_pool_2<- scale_custom_bounds(forage_raw_2, 0.0, 1.0)

minmax_df$Peff_M2 <- minmax_df$X15_sablefishBiomass_predAK_smoltyr * env_gate_2 * (1 - 1.0 * forage_pool_2)

m2_sem <- sem(
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + Peff_M2"),
  data = minmax_df, std.lv = TRUE, missing = "ML", warn = FALSE
)

# --- MODEL 3: Sleeper Sharks Adultyr ---
# P_norm = X15_DFA_sleeperSharks_adultyr
# Env    = X21_sst_egoa_junjulaug_adultyr
# Prey   = X13_pollock_age1plus_adultyr + X13_capelin_WGoA_adultyr
env_gate_3   <- scale_custom_bounds(filtered_df$X21_sst_egoa_junjulaug_adultyr, a = 0.2, b = 1.0)
forage_raw_3 <- rowSums(minmax_df[, c("X13_pollock_age1plus_adultyr", "X13_capelin_WGoA_adultyr")], na.rm = TRUE)
forage_pool_3<- scale_custom_bounds(forage_raw_3, 0.0, 1.0)

minmax_df$Peff_M3 <- minmax_df$X15_DFA_sleeperSharks_adultyr * env_gate_3 * (1 - 1.0 * forage_pool_3)

m3_sem <- sem(
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + Peff_M3"),
  data = minmax_df, std.lv = TRUE, missing = "ML", warn = FALSE
)

# --- MODEL 4: Steller Sea Lion Pups (SEAK) ---
# P_norm = X11_ssl_seak_pup_pred
# Env    = X21_sst_egoa_junjulaug_smoltyr
# Prey   = X13_sitkaHerring_EGoA_smoltyr + X13_capelin_WGoA_smoltyr
env_gate_4   <- scale_custom_bounds(filtered_df$X21_sst_egoa_junjulaug_smoltyr, a = 0.2, b = 1.0)
forage_raw_4 <- rowSums(minmax_df[, c("X13_sitkaHerring_EGoA_smoltyr", "X13_capelin_WGoA_smoltyr")], na.rm = TRUE)
forage_pool_4<- scale_custom_bounds(forage_raw_4, 0.0, 1.0)

minmax_df$Peff_M4 <- minmax_df$X11_ssl_seak_pup_pred * env_gate_4 * (1 - 1.0 * forage_pool_4)

m4_sem <- sem(
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + Peff_M4"),
  data = minmax_df, std.lv = TRUE, missing = "ML", warn = FALSE
)

# ------------------------------------------------------------------------------
# STEP 2: PLOT INGREDIENTS VS PEFF FOR ALL 4 MODELS
# ------------------------------------------------------------------------------

make_ingredient_plot <- function(df_scaled, df_raw_env, title_str, pred_var, env_var, prey_vars, peff_var) {
  
  plot_df <- tibble(
    Year              = filtered_df$Year,
    Raw_Predator      = df_scaled[[pred_var]],
    Env_Gate          = scale_custom_bounds(df_raw_env[[env_var]], a = 0.2, b = 1.0),
    AltPrey_Pool      = scale_custom_bounds(rowSums(df_scaled[, prey_vars, drop = FALSE], na.rm = TRUE), 0.0, 1.0),
    Effective_Peff    = df_scaled[[peff_var]]
  ) %>%
    pivot_longer(-Year, names_to = "Component", values_to = "Scaled_Value") %>%
    mutate(
      Component = factor(Component, levels = c("Raw_Predator", "Env_Gate", "AltPrey_Pool", "Effective_Peff")),
      Is_Peff   = ifelse(Component == "Effective_Peff", "Top (Peff)", "Ingredient")
    )
  
  ggplot(plot_df, aes(x = Year, y = Scaled_Value, color = Component, size = Is_Peff, linetype = Is_Peff)) +
    geom_line() +
    geom_point(data = filter(plot_df, Component == "Effective_Peff"), size = 2) +
    scale_color_manual(values = c(
      "Raw_Predator"   = "#E69F00", # Orange
      "Env_Gate"       = "#56B4E9", # Sky Blue
      "AltPrey_Pool"   = "#009E73", # Green
      "Effective_Peff" = "#D55E00"  # Thick Red/Dark Orange
    )) +
    scale_size_manual(values = c("Ingredient" = 0.7, "Top (Peff)" = 1.6)) +
    scale_linetype_manual(values = c("Ingredient" = "dashed", "Top (Peff)" = "solid")) +
    labs(
      title = title_str,
      x = "Year", y = "Scaled Value [0, 1]"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

p1 <- make_ingredient_plot(
  minmax_df, filtered_df, 
  "1. Harbour Seals Adultyr", 
  "X10_Harbour_s_2yrLead_WS_adultyr", "X01_habCompInd_adultyr", 
  c("X12_DFA_biomassEuphShelfSum_adultyr", "X04_marketsquid_GAM_adultyr"), "Peff_M1"
)

p2 <- make_ingredient_plot(
  minmax_df, filtered_df, 
  "2. Sablefish Biomass Smoltyr", 
  "X15_sablefishBiomass_predAK_smoltyr", "X21_sst_egoa_junjulaug_smoltyr", 
  c("X13_pollock_age1plus_smoltyr", "X13_capelin_WGoA_smoltyr"), "Peff_M2"
)

p3 <- make_ingredient_plot(
  minmax_df, filtered_df, 
  "3. Sleeper Sharks Adultyr", 
  "X15_DFA_sleeperSharks_adultyr", "X21_sst_egoa_junjulaug_adultyr", 
  c("X13_pollock_age1plus_adultyr", "X13_capelin_WGoA_adultyr"), "Peff_M3"
)

p4 <- make_ingredient_plot(
  minmax_df, filtered_df, 
  "4. SSL SEAK Pup Pred", 
  "X11_ssl_seak_pup_pred", "X21_sst_egoa_junjulaug_smoltyr", 
  c("X13_sitkaHerring_EGoA_smoltyr", "X13_capelin_WGoA_smoltyr"), "Peff_M4"
)

combined_ingredients <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Ingredient Breakdown vs. Effective Predator Load (Peff)",
    subtitle = "Dashed lines = Raw Ingredients | Thick Solid Line = Final Gated Peff Load",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(combined_ingredients)
plot_file <- file.path(plot_out_dir, "Top4_Joint_Models_Ingredients.png")
ggsave(plot_file, combined_ingredients, width = 14, height = 10, dpi = 300)

cat(sprintf("\nIngredient plots successfully saved to: %s\n", plot_file))