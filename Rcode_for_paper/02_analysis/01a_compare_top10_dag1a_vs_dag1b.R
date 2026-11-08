

# ==============================================================================
# Script: 02_analysis/01a_compare_top10_dag1a_vs_dag1b.R
# Purpose: Compare Top 10 Valid Models for DAG 1A vs DAG 1B under liberal fit thresholds
# Output: Rcode_for_paper/Outputs/tables/obj1_table1_aic_dag1a_v_dag1b.[pdf|png]
# ==============================================================================

library(tidyverse)
library(gt)

# ------------------------------------------------------------------------------
# 1. Directory Setup (Long & Short Sweeps)
# ------------------------------------------------------------------------------
base_sweep_dir <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26"

dir_1a_long  <- file.path(base_sweep_dir, "DAG1A_long")
dir_1b_long  <- file.path(base_sweep_dir, "DAG1B_long")
dir_1a_short <- file.path(base_sweep_dir, "DAG1A_short")
dir_1b_short <- file.path(base_sweep_dir, "DAG1B_short")

out_dir     <- "Rcode_for_paper/Outputs"
tbl_out_dir <- file.path(out_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)

# Read CSV files
res_1a_long_raw  <- read.csv(file.path(dir_1a_long,  "SEMresultsByClus.csv"), stringsAsFactors = FALSE)
res_1b_long_raw  <- read.csv(file.path(dir_1b_long,  "SEMresultsByClus.csv"), stringsAsFactors = FALSE)
res_1a_short_raw <- read.csv(file.path(dir_1a_short, "SEMresultsByClus.csv"), stringsAsFactors = FALSE)
res_1b_short_raw <- read.csv(file.path(dir_1b_short, "SEMresultsByClus.csv"), stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 2. Filter Valid Models (Liberal Thresholds)
#    - p_Chi2 >= 0.05 (Bollen 1989)
#    - CFI    >= 0.80 (Bentler 1990; Marsh et al. 2004)
#    - SRMR   <= 0.10 (Hu & Bentler 1998; Kline 2016)
#    - AGFI   >= 0.80 (Cole 1987; Schermelleh-Engel et al. 2003)
# ------------------------------------------------------------------------------
filter_valid_sem <- function(df) {
  df %>%
    filter(
      !is.na(AIC),
      p_Chi2 >= 0.05,
      CFI    >= 0.80,
      SRMR   <= 0.10,
      AGFI   >= 0.80
    ) %>%
    arrange(AIC)
}

valid_1a_long  <- filter_valid_sem(res_1a_long_raw)
valid_1b_long  <- filter_valid_sem(res_1b_long_raw)
valid_1a_short <- filter_valid_sem(res_1a_short_raw)
valid_1b_short <- filter_valid_sem(res_1b_short_raw)

cat(sprintf("LONG Valid Models:  DAG 1A = %d | DAG 1B = %d\n", nrow(valid_1a_long), nrow(valid_1b_long)))
cat(sprintf("SHORT Valid Models: DAG 1A = %d | DAG 1B = %d\n", nrow(valid_1a_short), nrow(valid_1b_short)))

# ------------------------------------------------------------------------------
# 3. Build Formatted Panel Data Frame (Dropping Model Index Column)
# ------------------------------------------------------------------------------
prep_panel_df <- function(valid_a, valid_b, panel_name, top_n = 10) {
  top_a <- valid_a %>% slice_head(n = top_n)
  top_b <- valid_b %>% slice_head(n = top_n)
  
  tibble(
    Panel                 = panel_name,
    Rank                  = 1:top_n,
    `1A AIC`              = round(top_a$AIC, 2),
    `1A CFI`              = round(top_a$CFI, 3),
    `1A SRMR`             = round(top_a$SRMR, 3),
    `1A AGFI`             = round(top_a$AGFI, 3),
    `1A p_Chi2`           = round(top_a$p_Chi2, 3),
    `1B AIC`              = round(top_b$AIC, 2),
    `1B CFI`              = round(top_b$CFI, 3),
    `1B SRMR`             = round(top_b$SRMR, 3),
    `1B AGFI`             = round(top_b$AGFI, 3),
    `1B p_Chi2`           = round(top_b$p_Chi2, 3),
    `Delta AIC (1A - 1B)` = round(top_a$AIC - top_b$AIC, 2)
  )
}

df_long_panel  <- prep_panel_df(valid_1a_long,  valid_1b_long,  "Panel A: Long / 24-Year Sweep (1998–2021)")
df_short_panel <- prep_panel_df(valid_1a_short, valid_1b_short, "Panel B: Short / 19-Year Sweep (2003–2021)")
stacked_df     <- bind_rows(df_long_panel, df_short_panel)

# ------------------------------------------------------------------------------
# 4. Construct Unified gt Table with Tight Column Packing
# ------------------------------------------------------------------------------
unified_gt <- stacked_df %>%
  gt(groupname_col = "Panel") %>%
  tab_header(
    title = md("**Goodness-of-Fit and Information-Theoretic Comparison of Top 10 Models: DAG 1A vs. DAG 1B**")
  ) %>%
  tab_spanner(label = "DAG 1A (Sequential Mediation)", columns = starts_with("1A")) %>%
  tab_spanner(label = "DAG 1B (Direct Growth)",        columns = starts_with("1B")) %>%
  cols_label(
    Rank                  = md("**Rank**"),
    `1A AIC`              = md("**AIC**"),
    `1A CFI`              = md("**CFI**"),
    `1A SRMR`             = md("**SRMR**"),
    `1A AGFI`             = md("**AGFI**"),
    `1A p_Chi2`           = md("**p(&chi;<sup>2</sup>)**"),
    `1B AIC`              = md("**AIC**"),
    `1B CFI`              = md("**CFI**"),
    `1B SRMR`             = md("**SRMR**"),
    `1B AGFI`             = md("**AGFI**"),
    `1B p_Chi2`           = md("**p(&chi;<sup>2</sup>)**"),
    `Delta AIC (1A - 1B)` = md("**&Delta;AIC<br>(1A &minus; 1B)**")
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  fmt_number(columns = c(`1A AIC`, `1B AIC`, `Delta AIC (1A - 1B)`), decimals = 2) %>%
  fmt_number(columns = c(`1A CFI`, `1A SRMR`, `1A AGFI`, `1A p_Chi2`, `1B CFI`, `1B SRMR`, `1B AGFI`, `1B p_Chi2`), decimals = 3) %>%
  cols_width(
    Rank ~ px(45),
    starts_with("1A") ~ px(58),
    starts_with("1B") ~ px(58),
    `Delta AIC (1A - 1B)` ~ px(75)
  ) %>%
  tab_source_note(
    source_note = md(
      "**Filtering Criteria:** Top 10 candidate models filtered using liberal fit criteria: p(&chi;<sup>2</sup>) &ge; 0.05, CFI &ge; 0.80, SRMR &le; 0.10, AGFI &ge; 0.80. Positive &Delta;AIC values favor the direct growth topology of DAG 1B."
    )
  ) %>%
  tab_options(
    table.font.size            = px(12),
    heading.title.font.size    = px(13.5),
    row_group.font.size        = px(12.5),
    row_group.font.weight      = "bold",
    row_group.padding          = px(3),
    column_labels.font.size    = px(12),
    column_labels.font.weight  = "bold",
    data_row.padding           = px(2.5),
    source_notes.font.size     = px(9.5),
    source_notes.padding       = px(4)
  )

# ------------------------------------------------------------------------------
# 5. Export Unified Table (Tight Canvas)
# ------------------------------------------------------------------------------
base_filename <- "obj1_table1_aic_dag1a_v_dag1b"

pdf_target <- file.path(tbl_out_dir, paste0(base_filename, ".pdf"))
png_target <- file.path(tbl_out_dir, paste0(base_filename, ".png"))

gtsave(unified_gt, filename = pdf_target)
gtsave(unified_gt, filename = png_target, vwidth = 780, vheight = 800)