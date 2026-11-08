# ==============================================================================
# Script: plot_top4_joint_sem_and_ingredients.R
# Purpose: Fit lavaan models for top Joint (SST + AltPrey) predators,
#          including model variant type and prey lists in plot labels.
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

raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df       <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]

# Filter 24-year complete case columns
valid_cols  <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
filtered_df <- raw_df %>% select(all_of(valid_cols)) %>% filter(complete.cases(.))

scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

minmax_df <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

# Standardized Target Column Names
target_cpue <- "X07_cpue_IntSprJunHW_smoltyr"
target_sar  <- "X16_DFA_sar_UCsprCh_smoltyr"
hake_col    <- "X09_DFA_HakeAge5Plus_smoltyr"

# Standardized SST Predictor Variables
sst_smolt   <- "X21_sst_egoa_junjulaug_smoltyr"
sst_adult   <- "X21_sst_egoa_junjulaug_adultyr"

# ------------------------------------------------------------------------------
# STEP 1: CALCULATE PEFF AND FIT THE 4 JOINT MODELS
# ------------------------------------------------------------------------------

# --- MODEL 1: Harbour Seal (DFO Coastwide) Smoltyr ---
env_gate_1    <- scale_custom_bounds(filtered_df[[sst_smolt]], a = 0.2, b = 1.0)
forage_pool_1 <- minmax_df[["X05_DFA_abundSardine_smoltyr"]]

minmax_df$Peff_M1 <- minmax_df$X10_harbour_seal_DFO_Coastwide_smoltyr * env_gate_1 * (1 - 1.0 * forage_pool_1)

m1_sem <- sem(
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + Peff_M1"),
  data = minmax_df, std.lv = TRUE, missing = "ML", warn = FALSE
)

# --- MODEL 2: Sablefish Biomass (3+) Smoltyr ---
env_gate_2    <- scale_custom_bounds(filtered_df[[sst_smolt]], a = 0.2, b = 1.0)
forage_raw_2  <- rowSums(minmax_df[, c("X13_pollockGOA_age1and2_smoltyr", "sitkaHerring_EGoA_smoltyr")], na.rm = TRUE)
forage_pool_2 <- scale_custom_bounds(forage_raw_2, 0.0, 1.0)

minmax_df$Peff_M2 <- minmax_df$X15_sablefishAK_biomass3plus_smoltyr * env_gate_2 * (1 - 1.0 * forage_pool_2)

m2_sem <- sem(
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + Peff_M2"),
  data = minmax_df, std.lv = TRUE, missing = "ML", warn = FALSE
)

# --- MODEL 3: Sleeper Shark BSAI Smoltyr ---
env_gate_3    <- scale_custom_bounds(filtered_df[[sst_smolt]], a = 0.2, b = 1.0)
forage_pool_3 <- minmax_df[["capelin_WGoA_smoltyr"]]

minmax_df$Peff_M3 <- minmax_df$X15_DFA_sleeperSharkBSAI_predAK_smoltyr * env_gate_3 * (1 - 1.0 * forage_pool_3)

m3_sem <- sem(
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + Peff_M3"),
  data = minmax_df, std.lv = TRUE, missing = "ML", warn = FALSE
)

# --- MODEL 4: Halibut Biomass (Age 8+) Adultyr ---
env_gate_4    <- scale_custom_bounds(filtered_df[[sst_adult]], a = 0.2, b = 1.0)
forage_pool_4 <- minmax_df[["capelin_WGoA_adultyr"]]

minmax_df$Peff_M4 <- minmax_df$X15_DFA_halibutBiomassAge8plus_predAK_adultyr * env_gate_4 * (1 - 1.0 * forage_pool_4)

m4_sem <- sem(
  paste0(target_cpue, " ~ ", hake_col, "\n", target_sar, " ~ ", target_cpue, " + Peff_M4"),
  data = minmax_df, std.lv = TRUE, missing = "ML", warn = FALSE
)

# ------------------------------------------------------------------------------
# STEP 2: PLOT INGREDIENTS WITH VARIANT & PREY LABELS
# ------------------------------------------------------------------------------

make_ingredient_plot <- function(df_scaled, df_raw_env, title_str, variant_str, pred_var, env_var, prey_vars, peff_var) {
  
  # Format prey list for subtitle
  clean_prey_names <- gsub("^(X[0-9]{2}_|DFA_)", "", prey_vars)
  prey_label <- paste(clean_prey_names, collapse = " + ")
  subtitle_str <- paste0("Variant: ", variant_str, " | Prey: ", prey_label)
  
  plot_df <- tibble(
    year              = filtered_df$year,
    Raw_Predator      = df_scaled[[pred_var]],
    Env_Gate          = scale_custom_bounds(df_raw_env[[env_var]], a = 0.2, b = 1.0),
    AltPrey_Pool      = scale_custom_bounds(rowSums(df_scaled[, prey_vars, drop = FALSE], na.rm = TRUE), 0.0, 1.0),
    Effective_Peff    = df_scaled[[peff_var]]
  ) %>%
    pivot_longer(-year, names_to = "Component", values_to = "Scaled_Value") %>%
    mutate(
      Component = factor(Component, levels = c("Raw_Predator", "Env_Gate", "AltPrey_Pool", "Effective_Peff")),
      Is_Peff   = ifelse(Component == "Effective_Peff", "Top (Peff)", "Ingredient")
    )
  
  ggplot(plot_df, aes(x = year, y = Scaled_Value, color = Component, size = Is_Peff, linetype = Is_Peff)) +
    geom_line() +
    geom_point(data = filter(plot_df, Component == "Effective_Peff"), size = 2) +
    scale_color_manual(values = c(
      "Raw_Predator"   = "#E69F00", # Orange
      "Env_Gate"       = "#56B4E9", # Sky Blue
      "AltPrey_Pool"   = "#009E73", # Green
      "Effective_Peff" = "#D55E00"  # Red/Dark Orange
    )) +
    scale_size_manual(values = c("Ingredient" = 0.7, "Top (Peff)" = 1.6)) +
    scale_linetype_manual(values = c("Ingredient" = "dashed", "Top (Peff)" = "solid")) +
    labs(
      title    = title_str,
      subtitle = subtitle_str,
      x = "Year", y = "Scaled Value [0, 1]"
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, face = "italic", color = "gray20"),
      legend.position = "bottom",
      legend.title  = element_blank()
    )
}

p1 <- make_ingredient_plot(
  minmax_df, filtered_df, 
  title_str   = "1. Harbour Seal (DFO Coastwide - Smoltyr)", 
  variant_str = "Joint SST + 1 AltPrey",
  pred_var    = "X10_harbour_seal_DFO_Coastwide_smoltyr", 
  env_var     = sst_smolt, 
  prey_vars   = c("X05_DFA_abundSardine_smoltyr"), 
  peff_var    = "Peff_M1"
)

p2 <- make_ingredient_plot(
  minmax_df, filtered_df, 
  title_str   = "2. Sablefish Biomass 3+ (Smoltyr)", 
  variant_str = "Joint SST + 2 AltPrey",
  pred_var    = "X15_sablefishAK_biomass3plus_smoltyr", 
  env_var     = sst_smolt, 
  prey_vars   = c("X13_pollockGOA_age1and2_smoltyr", "sitkaHerring_EGoA_smoltyr"), 
  peff_var    = "Peff_M2"
)

p3 <- make_ingredient_plot(
  minmax_df, filtered_df, 
  title_str   = "3. Sleeper Shark BSAI (Smoltyr)", 
  variant_str = "Joint SST + 1 AltPrey",
  pred_var    = "X15_DFA_sleeperSharkBSAI_predAK_smoltyr", 
  env_var     = sst_smolt, 
  prey_vars   = c("capelin_WGoA_smoltyr"), 
  peff_var    = "Peff_M3"
)

p4 <- make_ingredient_plot(
  minmax_df, filtered_df, 
  title_str   = "4. Halibut Biomass Age 8+ (Adultyr)", 
  variant_str = "Joint SST + 1 AltPrey",
  pred_var    = "X15_DFA_halibutBiomassAge8plus_predAK_adultyr", 
  env_var     = sst_adult, 
  prey_vars   = c("capelin_WGoA_adultyr"), 
  peff_var    = "Peff_M4"
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