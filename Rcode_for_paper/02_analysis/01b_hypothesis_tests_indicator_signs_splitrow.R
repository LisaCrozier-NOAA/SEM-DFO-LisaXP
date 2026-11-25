
#Summary
# Total number of models
# model_id    n_all
# 1 DAG1A_long  15552
# 2 DAG1A_short 72576
# 3 DAG1B_long  15552
# 4 DAG1B_short 72576

# Total number of valid models
#   model_id    n_all
# 1 DAG1A_long    139
# 2 DAG1A_short  2301
# 3 DAG1B_long    137
# 4 DAG1B_short  1983

#Number of models w/ dAIC<=5
# valid_keys_daic5 %>%   count(model_id, name = "n_after_exclusions")
#   model_id    n_after_exclusions
# 1 DAG1A_long                  26
# 2 DAG1A_short                 15
# 3 DAG1B_long                   2
# 4 DAG1B_short                 20


# ==============================================================================
# Script: 02_analysis/01b_hypothesis_tests_indicator_signs.R
# Purpose: Generate Main Table 2 and Supplemental Table S2 under dAIC <= 5.0
# Outputs:
#   - Outputs/tables/obj1_table2_hypo_path_tests.[pdf|png]
#   - Outputs/tables/obj1_tableS2_hypo_path_tests_by_indicator.[pdf|png]
#   - Outputs/data/obj1_datasummary_table2_hypo_path_tests_by_hypo.csv
#   - Outputs/data/obj1_data_table2_hypo_path_tests_by_indicator.csv
# ==============================================================================

library(tidyverse)
library(gt)
library(stringr)

# ------------------------------------------------------------------------------
# STEP 0: DIRECTORIES & CROSSWALK SETUP
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
doug_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26"

master_out   <- file.path(proj_dir, "Outputs")
tbl_out_dir  <- file.path(master_out, "tables")
data_out_dir <- file.path(master_out, "data")

dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

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

target_paths <- tribble(
  ~lhs_var,      ~rhs_var,     ~path_label,             ~H_code,
  "Growth",      "PreyNCC",    "PreyNCC -> Growth",     "H1a",
  "Abundance",   "Growth",     "Growth -> CPUE",        "H1b",
  "logSAR",      "Growth",     "Growth -> SAR",         "H1c",
  "logSAR",      "PreyAK",     "PreyAK -> SAR",         "H1e",
  "Abundance",   "PredNCC",    "PredNCC -> CPUE",       "H2a",
  "logSAR",      "PredAK",     "PredAK -> SAR",         "H2c",
  "logSAR",      "Abundance",  "CPUE -> SAR",           "H4a"
)

# ------------------------------------------------------------------------------
# STEP 1: LOAD RANKINGS & FILTER VALID MODELS (dAIC <= 5.0)
# ------------------------------------------------------------------------------
cat("Loading rankings and applying valid model filter...\n")
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

# Total number of models
rankings_all %>%
  anti_join(cr_seal_keys, by = c("model_id", "modNum")) %>%
  group_by(model_id) %>%
  mutate(dAIC = AIC - min(AIC, na.rm = TRUE)) %>%
  select(model_id, modNum, dAIC) %>%
  ungroup()%>%
  count(model_id, name = "n_all")
      # model_id    n_all
      # 1 DAG1A_long  15552
      # 2 DAG1A_short 72576
      # 3 DAG1B_long  15552
      # 4 DAG1B_short 72576

# Total number of valid models
rankings_all %>%
  anti_join(cr_seal_keys, by = c("model_id", "modNum")) %>%
  filter(
    !is.na(AIC),
    p_Chi2 >= 0.05,
    CFI    >= 0.80,
    SRMR   <= 0.10,
    AGFI   >= 0.80
  ) %>%
  group_by(model_id) %>%
  mutate(dAIC = AIC - min(AIC, na.rm = TRUE)) %>%
  select(model_id, modNum, dAIC) %>%
  ungroup()%>%
  count(model_id, name = "n_all")
#   model_id    n_all
# 1 DAG1A_long    139
# 2 DAG1A_short  2301
# 3 DAG1B_long    137
# 4 DAG1B_short  1983


# Valid model filter + dAIC <= 5.0 relative to best valid model
valid_keys_daic5 <- rankings_all %>%
  anti_join(cr_seal_keys, by = c("model_id", "modNum")) %>%
  filter(
    !is.na(AIC),
    p_Chi2 >= 0.05,
    CFI    >= 0.80,
    SRMR   <= 0.10,
    AGFI   >= 0.80
  ) %>%
  group_by(model_id) %>%
  mutate(dAIC = AIC - min(AIC, na.rm = TRUE)) %>%
  filter(dAIC <= 5.0) %>%
  select(model_id, modNum, dAIC) %>%
  ungroup()

#Number of models w/ dAIC<=5
valid_keys_daic5 %>%
  count(model_id, name = "n_after_exclusions")#the number of models in pool (w/o extra terms)
#   model_id    n_after_exclusions
# 1 DAG1A_long                  26
# 2 DAG1A_short                 15
# 3 DAG1B_long                   2
# 4 DAG1B_short                 20

# ------------------------------------------------------------------------------
# STEP 2: LOAD PARAMETERS & MERGE EVALUATED POOL
# ------------------------------------------------------------------------------
cat("Loading parameter estimates...\n")
estimates_all <- map_dfr(model_names, function(name) {
  f_path <- file.path(doug_dir, name, "parameterEstimates.csv")
  if (!file.exists(f_path)) return(tibble())
  read_csv(f_path, show_col_types = FALSE) %>% mutate(model_id = name)
})

estimates_base <- estimates_all %>%
  filter(op == "~") %>%
  inner_join(target_paths, by = c("lhs" = "lhs_var", "rhs" = "rhs_var")) %>%
  inner_join(rankings_long, by = c("model_id", "modNum", "rhs" = "generic_node")) %>%
  inner_join(valid_keys_daic5, by = c("model_id", "modNum")) %>%
  mutate(
    is_significant = if_else(pvalue < 0.05, 1, 0),
    is_positive    = if_else(est > 0, 1, 0),
    is_negative    = if_else(est < 0, 1, 0),
    is_competitor  = grepl("X02|X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", LisaName),
    Expected_Sign  = case_when(
      is_competitor                    ~ "Negative (-)",
      rhs %in% c("PredNCC", "PredAK")  ~ "Negative (-)",
      TRUE                             ~ "Positive (+)"
    ),
    Model_Supports = case_when(
      Expected_Sign == "Positive (+)" & is_positive == 1 ~ TRUE,
      Expected_Sign == "Negative (-)" & is_negative == 1 ~ TRUE,
      TRUE                                               ~ FALSE
    )
  )

# ------------------------------------------------------------------------------
# STEP 3: BUILD DATA & DISPLAY FOR TABLE 2 (SPLIT ROWS FOR PREY PATHS)
# ------------------------------------------------------------------------------
cat("Building Main Table 2 with split prey/competitor sub-rows...\n")

# 1. Tag individual estimates for split paths
estimates_split <- estimates_base %>%
  mutate(
    is_split_path = H_code %in% c("H1a", "H1e"),
    sub_guild = case_when(
      is_split_path & is_competitor  ~ "Competitor (-)",
      is_split_path & !is_competitor ~ "Prey (+)",
      TRUE                           ~ NA_character_
    ),
    Expected_Sign = case_when(
      sub_guild == "Competitor (-)" ~ "Negative (-)",
      sub_guild == "Prey (+)"       ~ "Positive (+)",
      H_code %in% c("H2a", "H2c")   ~ "Negative (-)",
      TRUE                          ~ "Positive (+)"
    ),
    Model_Supports = case_when(
      Expected_Sign == "Positive (+)" & is_positive == 1 ~ TRUE,
      Expected_Sign == "Negative (-)" & is_negative == 1 ~ TRUE,
      TRUE                                               ~ FALSE
    )
  )

# 2. Summarize (groups by sub_guild so H1a and H1e generate 2 rows each)
path_summary <- estimates_split %>%
  group_by(model_id, H_code, path_label, is_split_path, sub_guild, Expected_Sign) %>%
  summarise(
    N_Models      = n(),
    Pct_Sig       = round((sum(is_significant) / N_Models) * 100, 1),
    Pct_Positive  = round((sum(is_positive) / N_Models) * 100, 1),
    Pct_Supported = round((sum(Model_Supports) / N_Models) * 100, 1),
    Mean_Beta     = round(mean(est, na.rm = TRUE), 3),
    .groups       = "drop"
  ) %>%
  mutate(
    # Set Baseline_Expectation to show "Prey (+)" / "Competitor (-)" for split rows
    Baseline_Expectation = if_else(is_split_path, sub_guild, Expected_Sign),
    Verdict = case_when(
      Pct_Sig >= 50 & Pct_Supported >= 50 ~ "Supported",
      Pct_Sig >= 50 & Pct_Supported < 50  ~ "Flipped Sign",
      TRUE                                ~ "NS"
    ),
    Cell_Display = paste0(
      "<b>", Verdict, "</b><br>",
      Pct_Sig, "% Sig | ", Pct_Positive, "% +<br>",
      "<span style='font-size:9.5px; color:#555;'>(N=", N_Models, ")</span>"
    )
  )

write_csv(path_summary, file.path(data_out_dir, "obj1_datasummary_table2_hypo_path_tests_by_hypo_all_valid.csv"))
#write_csv(path_summary, file.path(data_out_dir, "obj1_datasummary_table2_hypo_path_tests_by_hypo_all_valid.csv"))

# 3. Pivot wide across models while preserving original column structure
t2_wide <- path_summary %>%
  select(H_code, path_label, Baseline_Expectation, is_split_path, model_id, Cell_Display) %>%
  pivot_wider(
    names_from = model_id, 
    values_from = Cell_Display, 
    values_fill = "Path Excluded"
  ) %>%
  mutate(
    # Order: H1a (Prey then Competitor), H1b, H1c, H1e (Prey then Competitor), H2a, H2c, H4a
    h_rank = match(H_code, c("H1a", "H1b", "H1c", "H1e", "H2a", "H2c", "H4a")),
    sub_rank = if_else(grepl("Competitor", Baseline_Expectation), 2, 1)
  ) %>%
  arrange(h_rank, sub_rank) %>%
  select(-h_rank, -sub_rank, -is_split_path)

# 4. Identify the second sub-row of split paths (where text in col 1 & 2 should visually merge)
dup_rows <- which(t2_wide$H_code == lag(t2_wide$H_code, default = ""))

# Blank out duplicate text in the first two columns so it acts like a vertically merged block
t2_display <- t2_wide
t2_display$H_code[dup_rows]     <- ""
t2_display$path_label[dup_rows] <- ""

cols_to_style <- intersect(model_names, names(t2_display))

all_labels <- list(
  H_code               = md("**Code**"),
  path_label           = md("**Structural Path**"),
  Baseline_Expectation = md("**Expected**"),
  DAG1A_long           = md("**Long (24-Yr)**"),
  DAG1A_short          = md("**Short (19-Yr)**"),
  DAG1B_long           = md("**Long (24-Yr)**"),
  DAG1B_short          = md("**Short (19-Yr)**")
)

# ------------------------------------------------------------------------------
# STEP 3B: RENDER PUBLICATION GT TABLE (RE-PROPORTIONED WIDTHS)
# ------------------------------------------------------------------------------
table2_gt <- t2_display %>%
  gt() %>%
  tab_header(
    title = md("**Table 2. Evaluation of Structural Path Hypotheses Across Plausible Models (&Delta;AIC &le; 5.0)**"),
    subtitle = "Directional verdict, % significant (p < 0.05), % positive coefficients, and candidate model count (N)."
  ) %>%
  tab_spanner(label = "Sequential Mediation (DAG 1A)", columns = any_of(c("DAG1A_long", "DAG1A_short"))) %>%
  tab_spanner(label = "Direct Growth (DAG 1B)",        columns = any_of(c("DAG1B_long", "DAG1B_short"))) %>%
  cols_label(.list = all_labels[names(all_labels) %in% names(t2_display)]) %>%
  fmt_markdown(columns = c(all_of(cols_to_style), "path_label")) %>%
  cols_align(align = "center", columns = c(H_code, Baseline_Expectation, all_of(cols_to_style))) %>%
  cols_align(align = "left", columns = path_label) %>%
  
  # --- TIGHTENED COLUMN WIDTHS ---
  cols_width(
    H_code               ~ px(42),
    path_label           ~ px(125),
    Baseline_Expectation ~ px(88),
    any_of(cols_to_style) ~ px(115)
  ) %>%
  
  # Remove the horizontal border between sub-rows to visually create the merged cell
  tab_style(
    style = cell_borders(sides = "top", color = "#E0E0E0", weight = px(0.5)),
    locations = cells_body(rows = dup_rows)
  ) %>%
  tab_style(
    style = cell_borders(sides = "top", color = "transparent"),
    locations = cells_body(columns = c(H_code, path_label), rows = dup_rows)
  ) %>%
  tab_source_note(
    source_note = md("Evaluated across valid models with &Delta;AIC &le; 5.0 relative to the top valid model per configuration. Supported: &ge;50% significant in hypothesized direction; Flipped Sign: &ge;50% significant in counter-intuitive direction; NS: <50% significant.")
  ) %>%
  tab_options(
    table.font.size            = px(11.5),
    heading.title.font.size    = px(13.5),
    heading.subtitle.font.size = px(11.0),
    column_labels.font.size    = px(11.5),
    column_labels.font.weight  = "bold",
    data_row.padding           = px(3.5),
    source_notes.font.size     = px(9.0)
  )

# Shading for Supported and Flipped Sign cells
for (col in cols_to_style) {
  table2_gt <- table2_gt %>%
    tab_style(
      style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#1B5E20")),
      locations = cells_body(columns = all_of(col), rows = str_detect(.data[[col]], "Supported"))
    ) %>%
    tab_style(
      style = list(cell_fill(color = "#FFEBEE"), cell_text(color = "#B71C1C")),
      locations = cells_body(columns = all_of(col), rows = str_detect(.data[[col]], "Flipped Sign"))
    )
}

# Export with wider viewport (900px) so the last column has plenty of clearance
gtsave(table2_gt, file.path(tbl_out_dir, "obj1_table2_hypo_path_tests.pdf"))
gtsave(table2_gt, file.path(tbl_out_dir, "obj1_table2_hypo_path_tests.png"), vwidth = 900, vheight = 620)
gtsave(table2_gt, file.path(tbl_out_dir, "obj1_table2_hypo_path_tests_all_valid.pdf"))
gtsave(table2_gt, file.path(tbl_out_dir, "obj1_table2_hypo_path_tests_all_valid.png"), vwidth = 900, vheight = 620)

cat("Table 2 re-exported with tightened path columns and full 4th column visibility.\n")

# ------------------------------------------------------------------------------
# STEP 4: BUILD DATA & DISPLAY FOR TABLE S2 (INDICATOR DIAGNOSTICS)
# ------------------------------------------------------------------------------
cat("Building Supplemental Table S2...\n")

tableS2_data <- estimates_base %>%
  group_by(H_code, path_label, rhs, LisaName) %>%
  summarise(
    N_Models     = n(),
    N_Supp       = sum(Model_Supports == TRUE),
    Pct_Sig_Supp = if_else(N_Supp > 0, round((sum(Model_Supports == TRUE & is_significant == 1) / N_Supp) * 100, 1), NA_real_),
    Beta_Supp    = if_else(N_Supp > 0, round(mean(est[Model_Supports == TRUE], na.rm = TRUE), 3), NA_real_),
    N_Flip       = sum(Model_Supports == FALSE),
    Pct_Sig_Flip = if_else(N_Flip > 0, round((sum(Model_Supports == FALSE & is_significant == 1) / N_Flip) * 100, 1), NA_real_),
    Beta_Flip    = if_else(N_Flip > 0, round(mean(est[Model_Supports == FALSE], na.rm = TRUE), 3), NA_real_),
    .groups      = "drop"
  ) %>%
  rename(SEM_Node = rhs, Indicator = LisaName) %>%
  arrange(H_code, SEM_Node, Indicator)

write_csv(tableS2_data, file.path(data_out_dir, "obj1_data_table2_hypo_path_tests_by_indicator.csv"))

tableS2_gt <- tableS2_data %>%
  select(H_code, Indicator, N_Models, Pct_Sig_Supp, Beta_Supp, Pct_Sig_Flip, Beta_Flip) %>%
  gt(groupname_col = "H_code") %>%
  tab_header(
    title = md("**Table S2. Individual Indicator Diagnostics Across Structural Paths (&Delta;AIC &le; 5.0)**"),
    subtitle = "Frequency and mean effect size (&beta;) split by directional agreement across valid models."
  ) %>%
  tab_spanner(label = md("**Supported Direction**"), columns = c(Pct_Sig_Supp, Beta_Supp)) %>%
  tab_spanner(label = md("**Flipped Sign Direction**"), columns = c(Pct_Sig_Flip, Beta_Flip)) %>%
  cols_label(
    Indicator    = md("**Indicator Subguild Index**"),
    N_Models     = md("**N**"),
    Pct_Sig_Supp = md("**% Sig**"),
    Beta_Supp    = md("**Mean &beta;**"),
    Pct_Sig_Flip = md("**% Sig**"),
    Beta_Flip    = md("**Mean &beta;**")
  ) %>%
  cols_align(align = "left", columns = Indicator) %>%
  cols_align(align = "center", columns = c(N_Models, Pct_Sig_Supp, Beta_Supp, Pct_Sig_Flip, Beta_Flip)) %>%
  cols_width(
    Indicator    ~ px(290),
    N_Models     ~ px(45),
    Pct_Sig_Supp ~ px(65),
    Beta_Supp    ~ px(75),
    Pct_Sig_Flip ~ px(65),
    Beta_Flip    ~ px(75)
  ) %>%
  sub_missing(columns = everything(), missing_text = "-") %>%
  tab_style(
    style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#2E7D32", weight = "bold")),
    locations = cells_body(columns = c(Pct_Sig_Supp, Beta_Supp), rows = !is.na(Pct_Sig_Supp) & Pct_Sig_Supp >= 50)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#FFEBEE"), cell_text(color = "#C62828", weight = "bold")),
    locations = cells_body(columns = c(Pct_Sig_Flip, Beta_Flip), rows = !is.na(Pct_Sig_Flip) & Pct_Sig_Flip >= 50)
  ) %>%
  tab_options(
    table.font.size            = px(11.5),
    heading.title.font.size    = px(13.5),
    heading.subtitle.font.size = px(11.0),
    column_labels.font.size    = px(11.5),
    column_labels.font.weight  = "bold",
    row_group.font.size        = px(12.0),
    row_group.font.weight      = "bold",
    row_group.padding          = px(2),
    data_row.padding           = px(2),
    source_notes.font.size     = px(9.5)
  )

gtsave(tableS2_gt, file.path(tbl_out_dir, "obj1_tableS2_hypo_path_tests_by_indicator.pdf"))
gtsave(tableS2_gt, file.path(tbl_out_dir, "obj1_tableS2_hypo_path_tests_by_indicator.png"), vwidth = 680, vheight = 1150)
# gtsave(tableS2_gt, file.path(tbl_out_dir, "obj1_tableS2_hypo_path_tests_by_indicator_all_valid.pdf"))
# gtsave(tableS2_gt, file.path(tbl_out_dir, "obj1_tableS2_hypo_path_tests_by_indicator_all_valid.png"), vwidth = 680, vheight = 1150)

cat("\nTightened tables generated in Outputs/tables!\n")


