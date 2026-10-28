# ==============================================================================
# Script: plot_sleeper_shark_stack.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Stacked barchart decomposing Sleeper Shark raw index into:
#          1) Base: Final Peff (Midnight Navy)
#          2) Middle: Reduction from Capelin buffer (Deep Green)
#          3) Top: Reduction from Thermal gate (Dark Orange)
#          Total bar height equals the original Scaled Raw Predator Index.
# Output: Rcode_for_paper/Routput_for_paper/plots/SleeperShark_Waterfall_Stack.png
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
# STEP 1: COMPUTE GATES AND DECOMPOSE REDUCTIONS
# ------------------------------------------------------------------------------
# Raw predator index
P_raw <- df_base[[pred_col]]

# Gates (bounded between 0.2 and 1.0)
G_temp <- scale_custom_bounds(df_raw[[sst_col]], a = 0.2, b = 1.0)
G_prey <- 1 - scale_custom_bounds(df_base[[prey_col]], a = 0.2, b = 1.0)

# Multiplicative steps
P_after_temp <- P_raw * G_temp
P_eff        <- P_after_temp * G_prey

# Additive slices:
# 1. Base: Effective load
# 2. Middle: Portion mitigated by Capelin
# 3. Top: Portion mitigated by Cold Temperature
delta_temp <- P_raw - P_after_temp
delta_prey <- P_after_temp - P_eff

# Verify exact reconstruction
stopifnot(all(abs((P_eff + delta_prey + delta_temp) - P_raw) < 1e-9))

# ==============================================================================
# Script: plot_sleeper_shark_waterfall_stack.R
# Purpose: Stacked barchart with Peff at the base (Midnight Navy),
#          Capelin buffer in the middle (Forest Green), and
#          Temperature reduction at the top (Dark Orange).
# ==============================================================================
# ------------------------------------------------------------------------------
# STEP 2 & 3: CORRECTED STACK ORDER & MAPPING
# ------------------------------------------------------------------------------
# 1. Order factor levels from Base to Top
stack_levels <- c(
  "Effective Predation (Peff)",              # 1. Physical Bottom
  "Reduced by Capelin (Alternative Prey)",   # 2. Middle
  "Reduced by Temperature (Thermal Refuge)"  # 3. Physical Top
)

stack_df <- tibble(
  year  = df_raw$year,
  P_raw = P_raw,
  `Effective Predation (Peff)`              = P_eff,
  `Reduced by Capelin (Alternative Prey)`   = delta_prey,
  `Reduced by Temperature (Thermal Refuge)` = delta_temp
) %>%
  pivot_longer(
    cols      = all_of(stack_levels),
    names_to  = "Component",
    values_to = "Amount"
  ) %>%
  mutate(Component = factor(Component, levels = stack_levels))

stack_colors <- c(
  "Effective Predation (Peff)"              = "#08306B", # Midnight Navy
  "Reduced by Capelin (Alternative Prey)"   = "#1B7837", # Forest Green
  "Reduced by Temperature (Thermal Refuge)" = "#D95F02"  # Warm Orange
)

# 2. Build plot with position_stack(reverse = TRUE)
p_stack <- ggplot() +
  geom_hline(yintercept = seq(0, 1, by = 0.2), color = "gray85", linewidth = 0.3) +
  geom_col(
    data        = stack_df,
    aes(x = year, y = Amount, fill = Component),
    position    = position_stack(reverse = TRUE), # <- Forces Level 1 (Peff) to the baseline
    width       = 0.68,
    color       = "white",
    linewidth   = 0.25,
    alpha       = 0.92
  ) +
  geom_step(
    data        = tibble(year = df_raw$year, P_raw = P_raw),
    aes(x = year, y = P_raw),
    direction   = "mid",
    color       = "gray30",
    linetype    = "dotted",
    linewidth   = 0.6
  ) +
  add_phylopic(
    uuid   = shark_uuid,
    x      = 2018.5,
    y      = 0.98,
    height = 0.14,
    fill   = "#08306B",
    alpha  = 0.85
  ) +
  annotate(
    "text", x = 2018.8, y = 0.89,
    label = "Pacific Sleeper Shark Model",
    size  = 3.4, fontface = "bold.italic", color = "#08306B"
  ) +
  scale_x_continuous(breaks = seq(1998, 2021, by = 2), expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = stack_colors) +
  guides(
    fill = guide_legend(
      reverse   = TRUE, # Matches visual top-to-bottom order (Orange on top, Navy on bottom)
      title     = "Component Partition",
      keywidth  = unit(1.2, "cm"),
      keyheight = unit(0.9, "cm")
    )
  ) +
  labs(
    title    = "Decomposition of Sleeper Shark Predation Index",
    subtitle = "Total stack height = Raw Shark Index | Base bar = Final Effective Predation Load (Peff)",
    x        = "Ocean Migration Year (Smolt Year)",
    y        = "Predation Index Partition [0, 1]"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title         = element_text(face = "bold", size = 16, color = "black"),
    plot.subtitle      = element_text(size = 11, face = "italic", color = "gray25", margin = margin(b = 12)),
    axis.title         = element_text(face = "bold", size = 13),
    axis.text          = element_text(size = 11, color = "black"),
    axis.line          = element_line(color = "black", linewidth = 0.5),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    legend.position    = "right",
    legend.title       = element_text(face = "bold", size = 12),
    legend.text        = element_text(size = 10.5, face = "bold"),
    legend.background  = element_rect(fill = "gray98", color = "gray80", linewidth = 0.5),
    plot.margin        = margin(15, 15, 15, 15)
  )

file_out <- file.path(plot_out_dir, "SleeperShark_Waterfall_Stack.png")
ggsave(file_out, p_stack, width = 13.33, height = 7.5, units = "in", dpi = 300, bg = "white")


# # ------------------------------------------------------------------------------
# # STEP 4: SAVE OUTPUT (16:9 Slide Proportions)
# # ------------------------------------------------------------------------------
# print(p_stack)
# 
# file_out <- file.path(plot_out_dir, "SleeperShark_Waterfall_Stack.png")
# ggsave(
#   filename = file_out,
#   plot     = p_stack,
#   width    = 13.33,
#   height   = 7.5,
#   units    = "in",
#   dpi      = 300,
#   bg       = "white"
# )
# 
# cat(sprintf("\n[Saved] Stacked waterfall figure saved to: %s\n", file_out))
