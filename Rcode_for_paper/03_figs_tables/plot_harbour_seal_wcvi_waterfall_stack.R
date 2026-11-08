# ==============================================================================
# Script: plot_harbour_seal_wcvi_waterfall_stack.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Stacked barchart decomposing WCVI Harbour Seal raw index into:
#          1) Base: Final Peff (Midnight Navy #08306B)
#          2) Top: Reduction from Euphausiids + Sardine buffer (Deep Forest Green #1B7837)
#          (Alt-prey only modulation; no thermal gate).
#          Total bar height equals the original Scaled Raw Predator Index.
# Output: Rcode_for_paper/Routput_for_paper/plots/HarbourSeal_WCVI_Waterfall_Stack.png
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
pred_col  <- "X10_harbour_seal_DFO_WCVI_smoltyr"
prey_cols <- c("biomassEuphShelfSum_smoltyr", "X05_DFA_abundSardine_smoltyr")

# ------------------------------------------------------------------------------
# STEP 1: COMPUTE GATES AND DECOMPOSE REDUCTIONS
# ------------------------------------------------------------------------------
# Raw predator index
P_raw <- df_base[[pred_col]]

# Alternative Prey Gate (Sum of Euphausiids + Sardine, a = 0.2, b = 1.0)
raw_forage_seal <- rowSums(df_base[, prey_cols], na.rm = TRUE)
G_prey          <- 1 - scale_custom_bounds(raw_forage_seal, a = 0.2, b = 1.0)

# Multiplicative step (Alt-prey only)
P_eff <- P_raw * G_prey

# Additive slices:
# 1. Base: Effective load
# 2. Top: Portion mitigated by Euphausiids + Sardine buffer
delta_prey <- P_raw - P_eff

# Verify exact reconstruction
stopifnot(all(abs((P_eff + delta_prey) - P_raw) < 1e-9))

# ------------------------------------------------------------------------------
# STEP 2: PREPARE STACKED DATA TABLE
# ------------------------------------------------------------------------------
stack_levels <- c(
  "Effective Predation (Peff)",                  # 1. Physical Bottom
  "Reduced by Euphausiids + Sardine (Alt Prey)"  # 2. Physical Top
)

stack_df <- tibble(
  year  = df_raw$year,
  P_raw = P_raw,
  `Effective Predation (Peff)`                 = P_eff,
  `Reduced by Euphausiids + Sardine (Alt Prey)` = delta_prey
) %>%
  pivot_longer(
    cols      = all_of(stack_levels),
    names_to  = "Component",
    values_to = "Amount"
  ) %>%
  mutate(Component = factor(Component, levels = stack_levels))

stack_colors <- c(
  "Effective Predation (Peff)"                 = "#08306B", # Deep Midnight Navy (Base)
  "Reduced by Euphausiids + Sardine (Alt Prey)" = "#1B7837"  # Deep Forest Green (Top)
)

seal_uuid <- get_uuid(name = "Phoca vitulina")

# ------------------------------------------------------------------------------
# STEP 3: ASSEMBLE PLOT
# ------------------------------------------------------------------------------
p_stack <- ggplot() +
  geom_hline(yintercept = seq(0, 1, by = 0.2), color = "gray85", linewidth = 0.3) +
  
  # Stacked columns: position_stack(reverse = TRUE) locks Peff to the bottom
  geom_col(
    data        = stack_df,
    aes(x = year, y = Amount, fill = Component),
    position    = position_stack(reverse = TRUE),
    width       = 0.68,
    color       = "white",
    linewidth   = 0.25,
    alpha       = 0.92
  ) +
  
  # Reference line tracing the total Raw Harbour Seal Index over the bars
  geom_step(
    data        = tibble(year = df_raw$year, P_raw = P_raw),
    aes(x = year, y = P_raw),
    direction   = "mid",
    color       = "gray30",
    linetype    = "dotted",
    linewidth   = 0.6
  ) +
  
  # Harbour Seal silhouette in upper margin
  add_phylopic(
    uuid       = seal_uuid,
    x          = 2018.5,
    y          = 0.96,
    height     = 0.14,
    fill       = "#08306B",
    alpha      = 0.85,
    horizontal = TRUE
  ) +
  annotate(
    "text", x = 2018.8, y = 0.86,
    label = "DFO WCVI Harbour Seal",
    size  = 3.4, fontface = "bold.italic", color = "#08306B"
  ) +
  
  # Scales & Legend
  scale_x_continuous(breaks = seq(1998, 2021, by = 2), expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = stack_colors) +
  
  guides(
    fill = guide_legend(
      reverse   = TRUE, # Keeps visual top-to-bottom legend matching the bars (Green top, Navy bottom)
      title     = "Component Partition",
      keywidth  = unit(1.2, "cm"),
      keyheight = unit(0.9, "cm")
    )
  ) +
  
  labs(
    title    = "Decomposition of WCVI Harbour Seal Predation Index",
    subtitle = "Total stack height = Raw DFO WCVI Seal Index | Base bar = Final Effective Predation Load (Peff)",
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

# ------------------------------------------------------------------------------
# STEP 4: SAVE OUTPUT (16:9 Slide Proportions)
# ------------------------------------------------------------------------------
print(p_stack)

file_out <- file.path(plot_out_dir, "HarbourSeal_WCVI_Waterfall_Stack.png")
ggsave(
  filename = file_out,
  plot     = p_stack,
  width    = 13.33,
  height   = 7.5,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

cat(sprintf("\n[Saved] Stacked waterfall figure saved to: %s\n", file_out))