# ==============================================================================
# Script: 00b_diagnostic_daic_cutoffs_VALID_ONLY.R
# Purpose: Inspect N valid models per dAIC threshold and plot valid distributions
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# 1. Setup Paths & Load Crosswalk
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
doug_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26"
fig_out_dir  <- file.path(proj_dir, "Outputs", "figures")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

crosswalk_path <- file.path(proj_dir, "metadata", "master_name_crosswalk.csv")
if (!file.exists(crosswalk_path)) {
  stop("Missing master crosswalk file at: ", crosswalk_path)
}
master_crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)
lookup_dict <- c(
  setNames(master_crosswalk$LisaName, master_crosswalk$dfa_cols),
  setNames(master_crosswalk$LisaName, master_crosswalk$DFAname)
)

model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short")

# ------------------------------------------------------------------------------
# 2. Load Raw Results and Apply Valid Model Filter
# ------------------------------------------------------------------------------
rankings_all <- map_dfr(model_names, function(name) {
  f_path <- file.path(doug_dir, name, "SEMresultsByClus.csv")
  if (!file.exists(f_path)) return(tibble())
  read_csv(f_path, show_col_types = FALSE) %>% mutate(model_id = name)
})

rankings_long <- rankings_all %>%
  select(model_id, modNum, AIC, ends_with("indNames")) %>%
  pivot_longer(
    cols = ends_with("indNames"),
    names_to = "generic_node",
    values_to = "var_name"
  ) %>%
  mutate(
    LisaName = case_when(
      var_name %in% names(lookup_dict) ~ lookup_dict[var_name],
      paste0("X", var_name) %in% names(lookup_dict) ~ lookup_dict[paste0("X", var_name)],
      TRUE ~ var_name
    )
  )

# Purge specific exclusions
cr_seal_keys <- rankings_long %>%
  filter(LisaName %in% c("X15_sablefishAK_biomass4plus", "X15_pollockGOA_SSB_lag1")) %>%
  select(model_id, modNum) %>%
  distinct()

# Filter for strictly VALID models first
valid_rankings <- rankings_all %>%
  anti_join(cr_seal_keys, by = c("model_id", "modNum")) %>%
  filter(
    !is.na(AIC),
    p_Chi2 >= 0.05,
    CFI    >= 0.80,
    SRMR   <= 0.10,
    AGFI   >= 0.80
  ) %>%
  group_by(model_id) %>%
  mutate(
    min_valid_AIC = min(AIC, na.rm = TRUE),
    dAIC_valid    = AIC - min_valid_AIC
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 3. Model Counts Across Various dAIC Cutoffs (Valid Models Only)
# ------------------------------------------------------------------------------
thresholds <- c(1, 2, 3, 4, 5, 7, 10)

cutoff_table_valid <- map_dfr(thresholds, function(cut) {
  valid_rankings %>%
    group_by(model_id) %>%
    summarise(
      Threshold = paste0("dAIC <= ", cut),
      N_Valid   = sum(dAIC_valid <= cut),
      .groups   = "drop"
    )
}) %>%
  pivot_wider(names_from = model_id, values_from = N_Valid) %>%
  arrange(as.numeric(str_extract(Threshold, "\\d+")))

cat("\n================ VALID MODEL COUNTS BY dAIC THRESHOLD ================\n")
print(as.data.frame(cutoff_table_valid))

# ------------------------------------------------------------------------------
# 4. Multi-Panel Histogram Plot (Valid Models)
# ------------------------------------------------------------------------------
p_hist_valid <- ggplot(valid_rankings , aes(x = dAIC_valid)) + #%>% filter(dAIC_valid <= 15)
  geom_histogram(binwidth = 0.5, fill = "#2E7D32", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 3, color = "#1976D2", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = 5, color = "#D32F2F", linetype = "dotted", linewidth = 0.9) +
  facet_wrap(~ model_id, scales = "free_y") +
  theme_minimal(base_size = 11) +
  labs(
    title = expression(paste("Distribution of VALID Candidate SEM Models across ", Delta, "AIC")),
    subtitle = "Blue dashed line: dAIC = 3 | Red dotted line: dAIC = 5 (Zoomed to dAIC <= 15)",
    x = expression(Delta * "AIC (relative to top valid model)"),
    y = "Number of Valid Models"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_out_dir, "obj1_fig_daic_histograms.png"), plot = p_hist_valid, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_out_dir, "obj1_fig_daic_histograms.pdf"), plot = p_hist_valid, width = 8, height = 6)

cat("\nValid model histograms saved to: ", fig_out_dir, "\n")