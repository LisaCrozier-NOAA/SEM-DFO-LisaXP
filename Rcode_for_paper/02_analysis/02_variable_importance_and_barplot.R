# ==============================================================================
# Script: 02_analysis/02_variable_importance_and_barplot.R
# Purpose: Objective 2 - Variable Importance (Akaike Weights) & Faceted Barchart
# Outputs:
#   - Outputs/figures/obj2_fig1_variable_importance_barchart.[pdf|png]
#   - Outputs/data/obj2_data_variable_importance_long.csv
#   - Outputs/data/obj2_data_variable_importance_wide.csv
# ==============================================================================

library(tidyverse)
library(ggh4x)
library(forcats)

# ------------------------------------------------------------------------------
# STEP 0: DEFINE PATHS & LOAD CROSSWALK
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
doug_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26"

master_out   <- file.path(proj_dir, "Outputs")
data_out_dir <- file.path(master_out, "data")
tbl_out_dir  <- file.path(master_out, "tables")
fig_out_dir  <- file.path(master_out, "figures")

dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)
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

cat("Running Objective 2 Variable Importance Analysis across Valid Models...\n")

# ------------------------------------------------------------------------------
# STEP 1: LOAD COMPLETENESS LOOKUP
# ------------------------------------------------------------------------------
completeness_file <- file.path(doug_dir, model_names[1], "completeness.csv")

if (file.exists(completeness_file)) {
  completeness_raw <- read_csv(completeness_file, show_col_types = FALSE)
  
  completeness_df <- completeness_raw %>%
    select(-1) %>%
    summarise(across(everything(), ~ !any(is.na(.)))) %>%
    pivot_longer(
      cols = everything(),
      names_to = "var",
      values_to = "complete"
    ) %>%
    left_join(master_crosswalk %>% select(DFAname, LisaName), by = c("var" = "DFAname")) %>%
    filter(!is.na(LisaName))
} else {
  warning("completeness.csv not found at: ", completeness_file)
}

# ------------------------------------------------------------------------------
# STEP 2: LOAD RANKINGS & FILTER VALID MODELS (2 BASELINE EXCLUSIONS)
# ------------------------------------------------------------------------------
rankings_all <- map_dfr(model_names, function(name) {
  f_path <- file.path(doug_dir, name, "SEMresultsByClus.csv")
  if (!file.exists(f_path)) return(tibble())
  read_csv(f_path, show_col_types = FALSE) %>% mutate(model_id = name)
})

rankings_long <- rankings_all %>%
  select(
    model_id, modNum, AIC, p_Chi2, CFI, SRMR, AGFI,
    PreyNCCindNames, PredNCCindNames, GrowthIndNames, 
    AbundanceIndNames, PreyAKindNames, PredAKindNames
  ) %>%
  pivot_longer(
    cols = ends_with("indNames"),
    names_to = "generic_node",
    values_to = "var_name"
  ) %>%
  mutate(
    generic_node = str_remove(generic_node, "(?i)indnames"),
    LisaName = case_when(
      var_name %in% names(lookup_dict) ~ lookup_dict[var_name],
      paste0("X", var_name) %in% names(lookup_dict) ~ lookup_dict[paste0("X", var_name)],
      TRUE ~ var_name
    )
  )

cr_seal_keys <- rankings_long %>%
  filter(LisaName %in% c("X15_sablefishAK_biomass4plus", "X15_pollockGOA_SSB_lag1")) %>%
  select(model_id, modNum) %>%
  distinct()

valid_rankings_long <- rankings_long %>%
  anti_join(cr_seal_keys, by = c("model_id", "modNum")) %>%
  filter(
    !is.na(AIC),
    p_Chi2 >= 0.05,
    CFI    >= 0.80,
    SRMR   <= 0.10,
    AGFI   >= 0.80
  )

write_csv(valid_rankings_long, file.path(data_out_dir, "obj2_data_valid_rankings_long.csv"))

# ------------------------------------------------------------------------------
# STEP 3: CALCULATE AKAIKE WEIGHTS ACROSS VALID MODELS
# ------------------------------------------------------------------------------
valid_models_with_weights <- valid_rankings_long %>%
  select(model_id, modNum, AIC) %>%
  distinct() %>%
  group_by(model_id) %>%
  mutate(
    delta_AIC     = AIC - min(AIC, na.rm = TRUE),
    rel_lik       = exp(-0.5 * delta_AIC),
    akaike_weight = rel_lik / sum(rel_lik, na.rm = TRUE)
  ) %>%
  ungroup()

# Long dataset (per model run and node)
variable_importance_long <- valid_models_with_weights %>%
  inner_join(valid_rankings_long, by = c("model_id", "modNum", "AIC")) %>%
  group_by(model_id, generic_node, LisaName) %>%
  summarise(
    Total_Selections_In_Ensemble = n(),
    Akaike_Importance_Weight     = round(sum(akaike_weight, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  mutate(
    Node2 = case_when(
      grepl("X02|X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", LisaName) ~ "Competitor",
      TRUE ~ as.character(generic_node)
    ),
    order = case_when(
      Node2 == "Abundance"  ~ 0,
      Node2 == "Growth"     ~ 0.5,
      Node2 == "Competitor" ~ 0.75,
      Node2 == "PreyNCC"    ~ 1,
      Node2 == "PredNCC"    ~ 2,
      Node2 == "PreyAK"     ~ 3,
      Node2 == "PredAK"     ~ 4,
      TRUE                  ~ 99
    )
  ) %>%
  arrange(model_id, order, desc(Akaike_Importance_Weight))

write_csv(variable_importance_long, file.path(data_out_dir, "obj2_data_variable_importance_long.csv"))

# Wide dataset (summary across all 4 configurations)
variable_importance_wide <- variable_importance_long %>%
  mutate(
    Region = if_else(grepl("^X0[1-9]", LisaName), "NCC", "AK"),
    Trophic = case_when(
      Node2 == "Abundance"  ~ "Abundance",
      Node2 == "Growth"     ~ "Growth",
      Node2 == "Competitor" ~ "Competitor",
      Node2 == "PreyNCC"    ~ "Prey",
      Node2 == "PredNCC"    ~ "Predator",
      Node2 == "PreyAK"     ~ "Prey",
      Node2 == "PredAK"     ~ "Predator",
      TRUE                  ~ "Other"
    )
  ) %>%
  select(model_id, LisaName, Region, Trophic, order, Akaike_Importance_Weight) %>%
  pivot_wider(
    names_from = model_id,
    values_from = Akaike_Importance_Weight
  ) %>%
  arrange(order, LisaName)

print(variable_importance_wide,n=Inf)
write_csv(variable_importance_wide, file.path(data_out_dir, "obj2_data_variable_importance_wide.csv"))

# ------------------------------------------------------------------------------
# STEP 4: FILTER BY max(importance) >= 0.05 FOR FIGURE
# ------------------------------------------------------------------------------
df_long <- variable_importance_wide %>%
  filter(order > 0) %>%
  pivot_longer(
    cols      = any_of(c("DAG1A_long", "DAG1B_long", "DAG1A_short", "DAG1B_short")),
    names_to  = c("variable", "type"),
    names_sep = "_",
    values_to = "importance"
  ) %>%
  group_by(LisaName) %>%
  filter(max(importance, na.rm = TRUE) >= 0.05) %>%
  ungroup() %>%
  mutate(
    Trophic    = fct_reorder(Trophic, order),
    Region     = fct_reorder(Region, order),
    plot_names = fct_rev(fct_reorder(LisaName, order))
  )

missing_labels <- df_long %>%
  filter(is.na(importance)) %>%
  left_join(completeness_df, by = c("plot_names" = "LisaName")) %>%
  filter(type == "long" & complete == FALSE) %>%
  group_by(plot_names, type, Trophic, Region) %>%
  summarize(label = "Missing data", .groups = "drop")

# ------------------------------------------------------------------------------
# STEP 5: GENERATE & EXPORT BARCHART
# ------------------------------------------------------------------------------
importance_barchart <- ggplot(df_long, aes(x = plot_names, y = importance, fill = variable)) +
  geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(
    values = c("DAG1A" = "#2E4053", "DAG1B" = "#D35400"),
    labels = c("DAG 1A (Sequential)", "DAG 1B (Direct Growth)")
  ) +
  geom_text(
    data        = missing_labels, 
    aes(x = plot_names, y = 0, label = label), 
    inherit.aes = FALSE,
    hjust       = -0.1, 
    color       = "grey40", 
    fontface    = "italic",
    size        = 3
  ) +
  coord_flip() + 
  facet_nested(
    Region + Trophic ~ type, 
    scales    = "free_y", 
    space     = "free_y",
    nest_line = element_line(color = "black"),
    resect    = unit(2, "mm")
  ) + 
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1.10),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal(base_size = 11) +
  labs(
    x    = "Ecosystem Indicator", 
    y    = "Relative Variable Importance (Akaike Weight across Valid Models)",
    fill = "DAG Topology:" 
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(), 
    strip.text.y       = element_text(angle = 0, face = "bold"),
    strip.text.x       = element_text(face = "bold"),
    panel.spacing.x    = unit(1.5, "lines"), 
    panel.spacing.y    = unit(0.8, "lines"),
    strip.background   = element_blank(),
    legend.position    = "bottom"
  )

ggsave(
  file.path(fig_out_dir, "obj2_fig1_variable_importance_barchart.png"),
  plot   = importance_barchart,
  width  = 10, 
  height = 8, 
  dpi    = 300
)

ggsave(
  file.path(fig_out_dir, "obj2_fig1_variable_importance_barchart.pdf"),
  plot   = importance_barchart,
  width  = 10, 
  height = 8
)

cat("Objective 2 products saved to Outputs/figures and Outputs/data!\n")
# ------------------------------------------------------------------------------
# STEP 6: EXPORT FORMATTED SUPPLEMENTAL TABLE (Table S3 / S-Variable Importance)
# ------------------------------------------------------------------------------
cat("Building Formatted Supplemental Table for Variable Importance...\n")

table_supp_vi <- variable_importance_wide %>%
  filter(order > 0) %>%
  select(Trophic, LisaName, DAG1A_long, DAG1A_short, DAG1B_long, DAG1B_short) %>%
  gt(groupname_col = "Trophic") %>%
  tab_header(
    title = md("**Table S3. Relative Variable Importance (Akaike Weights) Across All Candidate Subguild Indicators**"),
    subtitle = "Cumulative Akaike model weights (w_i) summed across valid structural equation models by DAG configuration."
  ) %>%
  tab_spanner(label = "Sequential Mediation (DAG 1A)", columns = c(DAG1A_long, DAG1A_short)) %>%
  tab_spanner(label = "Direct Growth (DAG 1B)",        columns = c(DAG1B_long, DAG1B_short)) %>%
  cols_label(
    LisaName    = md("**Indicator Subguild Index**"),
    DAG1A_long  = md("**Long (24-Yr)**"),
    DAG1A_short = md("**Short (19-Yr)**"),
    DAG1B_long  = md("**Long (24-Yr)**"),
    DAG1B_short = md("**Short (19-Yr)**")
  ) %>%
  fmt_number(
    columns = c(DAG1A_long, DAG1A_short, DAG1B_long, DAG1B_short),
    decimals = 3
  ) %>%
  sub_missing(columns = everything(), missing_text = "—") %>%
  cols_align(align = "left", columns = LisaName) %>%
  cols_align(align = "center", columns = starts_with("DAG")) %>%
  cols_width(
    LisaName ~ px(280),
    starts_with("DAG") ~ px(95)
  ) %>%
  tab_source_note(
    source_note = md("Weights represent the sum of Akaike weights across all statistically valid models within each configuration. Em-dash ('—') denotes series not evaluated in the 24-year complete-case dataset due to incomplete sampling duration.")
  ) %>%
  tab_options(
    table.font.size            = px(10.5),
    heading.title.font.size    = px(12.5),
    heading.subtitle.font.size = px(10.0),
    column_labels.font.size    = px(10.5),
    column_labels.font.weight  = "bold",
    row_group.font.size        = px(11.0),
    row_group.font.weight      = "bold",
    row_group.padding          = px(3),
    data_row.padding           = px(2),
    source_notes.font.size     = px(8.5)
  )

gtsave(table_supp_vi, file.path(tbl_out_dir, "obj2_tableS3_variable_importance_wide.pdf"))
gtsave(table_supp_vi, file.path(tbl_out_dir, "obj2_tableS3_variable_importance_wide.png"), vwidth = 720, vheight = 1300)
