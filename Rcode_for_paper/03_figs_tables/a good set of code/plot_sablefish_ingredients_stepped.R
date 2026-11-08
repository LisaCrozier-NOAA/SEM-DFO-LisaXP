# ==============================================================================
# Script: plot_sablefish_ingredients_stepped.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Generate 4 sequential slide figures building up the Sablefish (3+)
#          + Pollock/Herring + embSSL dual-arena model using rphylopic.
#          Step 4 de-emphasizes ingredients to highlight Peff bars.
# Output: Rcode_for_paper/Routput_for_paper/plots/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(rphylopic)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & DATA
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
data_out_dir <- file.path(proj_dir, "Routput_for_paper", "data")
plot_out_dir <- file.path(proj_dir, "Routput_for_paper", "plots")
dir.create(plot_out_dir, showWarnings = FALSE, recursive = TRUE)

raw_csv_path  <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
raw_df_master <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
raw_df_master <- raw_df_master[, names(raw_df_master) != "" & !is.na(names(raw_df_master))]

# 24-year complete cases
valid_cols <- names(raw_df_master)[colSums(!is.na(raw_df_master)) == 24]
df_raw <- raw_df_master %>% 
  select(all_of(valid_cols)) %>% 
  filter(complete.cases(.))

scale_custom_bounds <- function(x, a = 0.0, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

df_base <- df_raw %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

# Variable definitions
target_cpue      <- "X07_cpue_IntSprJunHW_smoltyr"
target_sar       <- "X16_DFA_sar_UCsprCh_smoltyr"
hake_col         <- "X09_DFA_HakeAge5Plus_smoltyr"
inriver_pred_col <- "X10_DFA_AllSeaLionsEMB_adultyr"
sst_col          <- "X21_sst_egoa_junjulaug_smoltyr"
pred_col         <- "X15_sablefishAK_biomass3plus_smoltyr"
prey_cols        <- c("X13_pollockGOA_age1and2_smoltyr", "sitkaHerring_EGoA_smoltyr")

# ------------------------------------------------------------------------------
# STEP 1: CONSTRUCT PEFF & FIT THE DUAL-ARENA SEM
# ------------------------------------------------------------------------------
# SST Gate: a = 0.2, w = 1.0
temp_gate <- scale_custom_bounds(df_raw[[sst_col]], a = 0.2, b = 1.0)

# Alternative Prey Gate (Sum of Pollock + Sitka Herring, a = 0.2, w = 1.0)
raw_forage_sable <- rowSums(df_base[, prey_cols], na.rm = TRUE)
prey_gate        <- 1 - scale_custom_bounds(raw_forage_sable, a = 0.2, b = 1.0)

# Effective Predation
df_base$P_eff <- df_base[[pred_col]] * temp_gate * prey_gate

# Fit Dual-Arena Model
sem_formula <- paste0(
  target_cpue, " ~ ", hake_col, "\n",
  target_sar,  " ~ ", target_cpue, " + P_eff + ", inriver_pred_col
)
fit_dual <- sem(sem_formula, data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE)

pe <- parameterEstimates(fit_dual, standardized = TRUE)
r2_sar  <- round(as.numeric(inspect(fit_dual, "r2")[[target_sar]]), 3)
aic_val <- round(as.numeric(fitMeasures(fit_dual, "aic")), 1)

pe_ocean   <- pe %>% filter(op == "~" & lhs == target_sar & rhs == "P_eff")
pe_inriver <- pe %>% filter(op == "~" & lhs == target_sar & rhs == inriver_pred_col)

stat_subtitle <- sprintf(
  "SAR R² = %.3f | AIC = %.1f | P_eff (β_std = %.2f, p = %.4f) | Bonneville Sea Lions (β_std = %.2f, p = %.3f)",
  r2_sar, aic_val, pe_ocean$std.all[1], pe_ocean$pvalue[1], pe_inriver$std.all[1], pe_inriver$pvalue[1]
)

# ------------------------------------------------------------------------------
# STEP 2: PREPARE DATA & PALETTES
# ------------------------------------------------------------------------------
component_levels <- c(
  "1. Raw Sablefish (3+)",
  "2. Summer SST (EGoA)",
  "3. Pollock + Herring Buffer",
  "4. Effective Predation (Peff)"
)

ts_df <- tibble(
  year                              = df_raw$year,
  `1. Raw Sablefish (3+)`           = df_base[[pred_col]],
  `2. Summer SST (EGoA)`            = temp_gate,
  `3. Pollock + Herring Buffer`     = scale_custom_bounds(raw_forage_sable, a = 0.2, b = 1.0),
  `4. Effective Predation (Peff)`   = df_base$P_eff
) %>%
  pivot_longer(-year, names_to = "Component", values_to = "Value") %>%
  mutate(Component = factor(Component, levels = component_levels))

# Unified palette: Deep Steel Blue (#2B547E) for predator, Midnight Navy (#08306B) for Peff
palette_colors <- c(
  "1. Raw Sablefish (3+)"           = "#2B547E", # Deep Steel Blue
  "2. Summer SST (EGoA)"            = "#C00000", # Crimson Red
  "3. Pollock + Herring Buffer"     = "#1B7837", # Deep Forest Green
  "4. Effective Predation (Peff)"   = "#08306B"  # Midnight Navy Blue
)

# Silhouette UUIDs via rphylopic
sable_uuid  <- get_uuid(name = "Anoplopoma fimbria") # Sablefish / Blackcod
herring_uuid <- get_uuid(name = "Clupea")             # Herring for forage buffer

# ------------------------------------------------------------------------------
# STEP 3: FUNCTION TO BUILD AND SAVE SEQUENTIAL PLOTS
# ------------------------------------------------------------------------------
make_stepped_plot <- function(active_components, step_num, step_label) {
  
  df_sub <- ts_df %>% filter(Component %in% active_components)
  
  # Step 4 knocks down background lines to spotlight Peff bars
  line_alpha <- if (step_num == 4) 0.28 else 0.95
  line_width <- if (step_num == 4) 0.95 else 1.50
  
  p <- ggplot() +
    geom_hline(yintercept = seq(0, 1, by = 0.2), color = "gray85", linewidth = 0.3)
  
  # Step 4: Peff bars (High opacity & prominence)
  if ("4. Effective Predation (Peff)" %in% active_components) {
    df_peff <- df_sub %>% filter(Component == "4. Effective Predation (Peff)")
    p <- p +
      geom_col(
        data = df_peff,
        aes(x = year, y = Value, fill = Component, color = Component),
        width = 0.62,
        alpha = 0.85,
        linewidth = 0.6
      )
  }
  
  # Steps 1 to 3: Ingredient lines
  df_lines <- df_sub %>% filter(Component != "4. Effective Predation (Peff)")
  if (nrow(df_lines) > 0) {
    p <- p +
      geom_line(
        data = df_lines,
        aes(x = year, y = Value, color = Component),
        linewidth = line_width,
        alpha = line_alpha
      )
  }
  
  # Sablefish silhouette & label ONLY on Step 1
  if (step_num %in% c(1,4)) {
    p <- p +
      add_phylopic(
        uuid   = sable_uuid,
        x      = 2018.5,
        y      = 0.98,
        height = 0.14,
        fill   = "#2B547E",
        alpha  = 0.85
      ) +
      annotate("text", x = 2018.8, y = 0.88, label = "Sablefish (Age 3+)", 
               size = 3.4, fontface = "bold.italic", color = "#2B547E")
  }
  
  # Herring/Forage silhouette on Step 3 (silent icon, disappears on Step 4)
  if (step_num == 3) {
    p <- p +
      add_phylopic(
        uuid   = herring_uuid,
        x      = 2018.5,
        y      = 0.35,
        height = 0.08,
        fill   = "#1B7837",
        alpha  = 0.75
      )
  }
  
  # Scales, Unified Legend, and Axis Limits
  p <- p +
    scale_x_continuous(breaks = seq(1998, 2021, by = 2), expand = expansion(mult = c(0.02, 0.04))) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.2)) +
    scale_color_manual(values = palette_colors, drop = FALSE) +
    scale_fill_manual(values = palette_colors, drop = FALSE) +
    guides(
      color = guide_legend(
        override.aes = list(
          linewidth = if (step_num == 4) c(0.95, 0.95, 0.95, 0) else c(1.5, 1.5, 1.5, 0),
          alpha     = if (step_num == 4) c(0.35, 0.35, 0.35, 1.0) else c(1.0, 1.0, 1.0, 1.0),
          fill      = c(NA, NA, NA, alpha("#08306B", 0.85)),
          shape     = c(NA, NA, NA, 22),
          size      = c(NA, NA, NA, 6)
        )
      ),
      fill = "none"
    ) +
    labs(
      title    = paste0("Model B Progression: Alaska Sablefish (Step ", step_num, "/4)"),
      subtitle = if (step_num == 4) stat_subtitle else paste0("Adding component: ", step_label),
      x        = "Ocean Migration Year (Smolt Year)",
      y        = "Scaled Environmental / Predation Index [0, 1]",
      color    = "Ingredient Sequence"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title         = element_text(face = "bold", size = 16, color = "black"),
      plot.subtitle      = element_text(size = 11, face = "italic", color = ifelse(step_num == 4, "#08306B", "gray30"), margin = margin(b = 10)),
      axis.title         = element_text(face = "bold", size = 13),
      axis.text          = element_text(size = 11, color = "black"),
      axis.line          = element_line(color = "black", linewidth = 0.5),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      
      # Side Legend Layout
      legend.position    = "right",
      legend.title       = element_text(face = "bold", size = 12),
      legend.text        = element_text(size = 11, face = "bold"),
      legend.key.width   = unit(1.2, "cm"),
      legend.key.height  = unit(0.9, "cm"),
      legend.background  = element_rect(fill = "gray98", color = "gray80", linewidth = 0.5),
      plot.margin        = margin(15, 15, 15, 15)
    )
  
  # 16:9 Widescreen slide dimensions
  file_png <- file.path(plot_out_dir, sprintf("Sablefish_Step%d_%s.png", step_num, gsub(" ", "", step_label)))
  ggsave(
    filename = file_png,
    plot     = p,
    width    = 13.33,
    height   = 7.5,
    units    = "in",
    dpi      = 300,
    bg       = "white"
  )
  cat(sprintf("[Saved for PPT] %s\n", basename(file_png)))
  
  return(p)
}

# ------------------------------------------------------------------------------
# STEP 4: GENERATE ALL 4 SEQUENTIAL SLIDES
# ------------------------------------------------------------------------------
p_step1 <- make_stepped_plot(component_levels[1],   1, "Raw Predator")
p_step2 <- make_stepped_plot(component_levels[1:2], 2, "Temperature Gate")
p_step3 <- make_stepped_plot(component_levels[1:3], 3, "Alternative Prey")
p_step4 <- make_stepped_plot(component_levels[1:4], 4, "Effective Load")

print(p_step1)
print(p_step2)
print(p_step3)
print(p_step4)
