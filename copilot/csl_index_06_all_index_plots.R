
#Load packages---------
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})



out_dir <- "copilot/outputs_7"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


guild.dfasAK_ishark_issl<-read.csv(   file.path("copilot/outputs_7", "guild.dfasAK_ishark_issl.csv"), row.names = NULL)


sort(names(guild.dfasAK_ishark_issl))
#Add csl index------

csl_path   <- file.path("copilot/outputs_6", "chinook_numeric_csl_exposure_risk_1998_2024.csv")

csl_raw <- read.csv(csl_path) %>% clean_names()

# ---------------------------
# 1) Build lag-2 CSL index (idxE10)
# For SAR year t, use CSL/chinook at t+2 adult return year
# ---------------------------
safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - m) / s)
}

csl_idx <- csl_raw %>%
  transmute(
    year,
    eul_z       = safe_z(eulachon_during_chinook),
    shad_z      = safe_z(shad_during_chinook),
    x15_icsl_risk_eulachon75_shad25 = -(0.75 * eul_z + 0.25 * shad_z),
    x15_icsl_cslchin_ratio = safe_z(log1p(csl_during_chinook) - log1p(chinook_estuary_total))
  ) %>%
  mutate(year = year - 2)  # align t+2 CSL to SAR year t


# ---------------------------
# 2) Assemble guild table as in original workflow------
# ---------------------------
guild_raw<-guild.dfasAK_ishark_issl

guild <- guild_raw %>%
  left_join(csl_idx, by = "year") %>%
  mutate(across(-year, ~ as.vector(scale(.))))

guild_dfas1_24yr <- guild %>%
  filter(year >= 1998, year <= 2021) %>%
  select(where(~ sum(!is.na(.)) >= 24))%>%
  rename(X15_ishark=ishark, X15_issl=issl) %>%
  select(-x15_shark_catch_go_a_pred_ak,-x10_harbor_seal_cr_2yr_lead)

sort(names(guild_dfas1_24yr))
write.csv(guild_dfas1_24yr,      file.path(out_dir, "guild_dfas1_24yr.csv"), row.names = FALSE)


#plot-----------
# ==============================================================================
# Script: Time Series Comparison - Transformed/Composite Indices vs. Raw Counts
# Categories:
# 1) Shark metrics
# 2) Steller Sea Lion (SSL), Herring (Herr), & Capelin (Cap) metrics
# 3) California Sea Lion (CSL) metrics
# ==============================================================================

library(tidyverse)
library(patchwork)
library(scales)

out_dir <- "copilot/outputs_6"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load Data & Helper Functions
# -----------------------------------------------------------------------------

# Load updated guild dataset
path_guild <- file.path(out_dir, "guild_dfas1_24yr.csv")
if (!file.exists(path_guild)) {
  path_guild <- file.path(out_dir, "guild_dfas1_24yr.csv")
}

guild_raw <- read.csv(path_guild)

# Safe z-score standardization helper (so counts & indices overlay smoothly)
safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - m) / s)
}

# -----------------------------------------------------------------------------
# 2. Category 1: Shark Metrics (Indices vs Raw Counts)
# -----------------------------------------------------------------------------

shark_cols <- names(guild_raw)[grepl("shark", names(guild_raw))]

df_shark <- guild_raw %>%
  select(year, all_of(shark_cols)) %>%
  pivot_longer(-year, names_to = "metric", values_to = "raw_val") %>%
  group_by(metric) %>%
  mutate(z_score = safe_z(raw_val)) %>%
  ungroup()

p_shark <- ggplot(df_shark, aes(x = year, y = z_score, color = metric, linetype = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  labs(
    title = "Category 1: Shark Metrics (Standardized Z-Score Comparison)",
    subtitle = "Comparing transformed/composite shark DFA indices against underlying raw count series",
    x = "Year",
    y = "Standardized Value (Z-Score)",
    color = "Metric",
    linetype = "Metric"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", color = "#0f172a"),
    panel.grid.minor = element_blank()
  )

print(p_shark)
ggsave(file.path(out_dir, "comparison_category1_sharks.png"), p_shark, width = 11, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 3. Category 2: SSL, Herring, & Capelin Metrics
# -----------------------------------------------------------------------------

cat2_pattern <- "ssl"
cat2_cols <- names(guild_raw)[grepl(cat2_pattern, names(guild_raw))]

df_cat2 <- guild_raw %>%
  select(year, all_of(cat2_cols)) %>%
  pivot_longer(-year, names_to = "metric", values_to = "raw_val") %>%
  group_by(metric) %>%
  mutate(z_score = safe_z(raw_val)) %>%
  ungroup()

p_cat2 <- ggplot(df_cat2, aes(x = year, y = z_score, color = metric, linetype = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  labs(
    title = "Category 2: SSL",
    subtitle = "Comparing Steller Sea Lion indices vs raw abundances",
    x = "Year",
    y = "Standardized Value (Z-Score)",
    color = "Metric",
    linetype = "Metric"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", color = "#0f172a"),
    panel.grid.minor = element_blank()
  )

print(p_cat2)
ggsave(file.path(out_dir, "comparison_category2_ssl.png"), p_cat2, width = 11, height = 6, dpi = 300)

# 3. Category 2: SSL, Herring, & Capelin Metrics
# -----------------------------------------------------------------------------
cat2_pattern <- "x13_egoa_herring|x13_stka_herr|x13_mid_il_capelin|x13_wgoa_cap_pcod"

cat2_cols <- names(guild_raw)[grepl(cat2_pattern, names(guild_raw))]

df_cat2 <- guild_raw %>%
  select(year, all_of(cat2_cols)) %>%
  pivot_longer(-year, names_to = "metric", values_to = "raw_val") %>%
  group_by(metric) %>%
  mutate(z_score = safe_z(raw_val)) %>%
  ungroup()

p_cat2 <- ggplot(df_cat2, aes(x = year, y = z_score, color = metric, linetype = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  labs(
    title = "Category 2: SSL, Herring, & Capelin Metrics (Standardized Z-Score Comparison)",
    subtitle = "Herring, and Capelin indices",
    x = "Year",
    y = "Standardized Value (Z-Score)",
    color = "Metric",
    linetype = "Metric"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", color = "#0f172a"),
    panel.grid.minor = element_blank()
  )

print(p_cat2)
ggsave(file.path(out_dir, "comparison_category2_herr_cap.png"), p_cat2, width = 11, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 4. Category 3: California Sea Lion (CSL) Metrics
# -----------------------------------------------------------------------------

csl_cols <- names(guild_raw)[grepl("csl|x10_dfa_ssl", names(guild_raw))]

df_csl <- guild_raw %>%
  select(year, all_of(csl_cols)) %>%
  pivot_longer(-year, names_to = "metric", values_to = "raw_val") %>%
  group_by(metric) %>%
  mutate(z_score = safe_z(raw_val)) %>%
  ungroup()

p_csl <- ggplot(df_csl, aes(x = year, y = z_score, color = metric, linetype = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  labs(
    title = "Category 3: California Sea Lion (CSL) Metrics (Standardized Z-Score Comparison)",
    subtitle = "Comparing log ratios (e.g., idxE10 CSL:Chinook) and risk indices vs raw CSL counts",
    x = "Year",
    y = "Standardized Value (Z-Score)",
    color = "Metric",
    linetype = "Metric"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", color = "#0f172a"),
    panel.grid.minor = element_blank()
  )

print(p_csl)
ggsave(file.path(out_dir, "comparison_category3_csl.png"), p_csl, width = 11, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 5. Combined Multi-Panel View for Presentation
# -----------------------------------------------------------------------------

p_combined_all <- (p_shark / p_cat2 / p_csl) + plot_layout(heights = c(1, 1, 1))

ggsave(
  file.path(out_dir, "comparison_all_3_categories_combined.png"),
  p_combined_all,
  width = 12, height = 14, dpi = 300
)
