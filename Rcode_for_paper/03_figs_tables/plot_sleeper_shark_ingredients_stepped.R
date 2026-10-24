# ==============================================================================
# Script: plot_sleeper_shark_ingredients_stepped.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Generate 4 sequential slide figures building up Sleeper Shark model.
#          Step 4 de-emphasizes ingredients to highlight Peff bars.
#          Predator line, silhouette, and Peff bars share a unified blue theme.
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
pred_col         <- "X15_DFA_sleeperSharkBSAI_predAK_smoltyr"
prey_col         <- "capelin_WGoA_smoltyr"

# ------------------------------------------------------------------------------
# STEP 1: CONSTRUCT PEFF & FIT THE DUAL-ARENA SEM
# ------------------------------------------------------------------------------
# SST Gate: a = 0.2, w = 1.0
temp_gate <- scale_custom_bounds(df_raw[[sst_col]], a = 0.2, b = 1.0)

# Alternative Prey Gate: a = 0.2, w = 1.0
prey_gate <- 1 - scale_custom_bounds(df_base[[prey_col]], a = 0.2, b = 1.0)

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
# STEP 2: PREPARE DATA & BLUE THEMED PALETTES
# ------------------------------------------------------------------------------
component_levels <- c(
  "1. Raw Sleeper Shark",
  "2. Summer SST (EGoA)",
  "3. Capelin Buffer (WGoA)",
  "4. Effective Predation (Peff)"
)

ts_df <- tibble(
  year                            = df_raw$year,
  `1. Raw Sleeper Shark`          = df_base[[pred_col]],
  `2. Summer SST (EGoA)`          = temp_gate,
  `3. Capelin Buffer (WGoA)`      = scale_custom_bounds(df_base[[prey_col]], a = 0.2, b = 1.0),
  `4. Effective Predation (Peff)` = df_base$P_eff
) %>%
  pivot_longer(-year, names_to = "Component", values_to = "Value") %>%
  mutate(Component = factor(Component, levels = component_levels))

# Unified palette: Steel blue for raw predator, deep navy for Peff
palette_colors <- c(
  "1. Raw Sleeper Shark"          = "#2B547E", # Steel Blue (lighter predator intensity)
  "2. Summer SST (EGoA)"          = "#C00000", # Crimson Red
  "3. Capelin Buffer (WGoA)"      = "#1B7837", # Deep Forest Green
  "4. Effective Predation (Peff)" = "#08306B"  # Midnight Navy Blue (deepest intensity)
)

# Fetch silhouettes
shark_uuid   <- get_uuid(name = "sharks")
capelin_uuid <- get_uuid(name = "capelin")

# ------------------------------------------------------------------------------
# STEP 3: FUNCTION TO BUILD AND SAVE SEQUENTIAL PLOTS
# ------------------------------------------------------------------------------
make_stepped_plot <- function(active_components, step_num, step_label) {
  
  df_sub <- ts_df %>% filter(Component %in% active_components)
  
  # When Step 4 is active, knock down the alpha and width of previous lines
  line_alpha <- if (step_num == 4) 0.28 else 0.95
  line_width <- if (step_num == 4) 0.95 else 1.50
  
  p <- ggplot() +
    geom_hline(yintercept = seq(0, 1, by = 0.2), color = "gray85", linewidth = 0.3)
  
  # Step 4: Draw Peff bars (High opacity & prominence)
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
  
  # Steps 1 to 3: Draw ingredient lines
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
  
  # Sleeper Shark silhouette & label ONLY on Step 1 (Matching Steel Blue hue)
  if (step_num %in% c(1,4)) {
    p <- p +
      add_phylopic(
        uuid   = shark_uuid,
        x      = 2018.5,
        y      = 1.00,
        height = 0.16,
        fill   = "#08306B",
        alpha  = 0.80
      ) +
      annotate("text", x = 2018.8, y = 0.89, label = "Pacific Sleeper Shark", 
               size = 3.4, fontface = "bold.italic", color = "#08306B")
  }
  
  # Capelin silhouette on Step 3 (omitted on Step 4 to keep focus on Peff)
  if (step_num == 3) {
    p <- p +
      add_phylopic(
        uuid   = capelin_uuid,
        x      = 2018.5,
        y      = 0.35,
        height = 0.09,
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
      title    = paste0("Model A Progression: Pacific Sleeper Shark (Step ", step_num, "/4)"),
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
  
  file_png <- file.path(plot_out_dir, sprintf("SleeperShark_Step%d_%s.png", step_num, gsub(" ", "", step_label)))
  ggsave(file_png, p, width = 12, height = 6.2, dpi = 300)
  cat(sprintf("[Saved] Step %d figure: %s\n", step_num, basename(file_png)))
  
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
print(p_step4 +
        add_phylopic(
          uuid   = shark_uuid,
          x      = 2018.5,
          y      = 1.0,
          height  = 0.16,
          fill   = "#08306B",
          alpha  = 0.75
        ))


#Muted & Slate Undertones (Closest in Character)
# Dark Slate Blue (#483D8B): Slightly purplish-gray, retaining the dusty, subdued feel of steel blue while being significantly deeper.
#   
#   Deep Steel Blue (#2B547E): A direct, darker step down in value from steel blue without shifting hue toward purple or green.
#     
#     Air Force Blue / Slate Navy (#27408B): A desaturated, classic military blue that sits right between slate and navy.
#       
#       Rich & Classic Undertones
#       Prussian Blue (#003153): A deep, rich, slightly teal-leaning dark blue historically used in marine blueprints and drafting.
#         
#         Midnight Blue (#191970): An extremely dark, near-black blue with subtle royal undertones.
#           
#           Navy Blue (#000080): The standard dark, pure navy with no gray dampening.
#             
#             For figures or data visualizations where steel blue is already in use, deep steel blue (#2B547E) or dark slate blue (#483D8B) will look the most cohesive without clashing.