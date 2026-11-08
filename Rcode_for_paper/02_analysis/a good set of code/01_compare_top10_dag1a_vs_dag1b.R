

# ==============================================================================
# Script: 02_analysis/01_compare_top10_dag1a_vs_dag1b.R
# Purpose: Compare Top 10 Valid Models for DAG 1A vs DAG 1B under liberal fit thresholds
# Output: Rcode_for_paper/Routput_for_paper/tables/
# ==============================================================================


#RESULTS-------
#DAG1b wins easily in the long model comparison (daic 3-7 for top 10 models), 
#so we will go with dag1B
# and they are essentially identical in the short comparison (daic -0.5 in top model, but >1 for models 2-6)

# LONG Valid Models:  DAG 1A = 171 | DAG 1B = 173
# SHORT Valid Models: DAG 1A = 2643 | DAG 1B = 2219
#total N models:  n_1b_long 18144 | n_1a_short 80640


# > print(as.data.frame(res_long$df))
# Rank 1A Mod 1A AIC 1A CFI 1A SRMR 1A AGFI 1A p_Chi2 1B Mod 1B AIC 1B CFI 1B SRMR 1B AGFI 1B p_Chi2 Delta AIC (1A - 1B)
# 1     1   8027 116.74  0.954   0.065   0.933     0.143   7861 109.03  0.940   0.076   0.919     0.082                7.71
# 2     2   8011 117.35  0.923   0.074   0.906     0.056   7863 110.33  0.946   0.077   0.925     0.102                7.01
# 3     3   7855 117.70  0.930   0.074   0.910     0.071  13909 110.97  0.969   0.049   0.907     0.208                6.74
# 4     4   8181 118.36  0.958   0.058   0.932     0.165  13911 112.26  0.973   0.048   0.911     0.231                6.09
# 5     5  14075 118.67  0.923   0.047   0.902     0.059   7853 114.44  0.940   0.079   0.936     0.093                4.23
# 6     6   7675 118.76  0.959   0.071   0.926     0.171   8349 114.51  0.935   0.071   0.918     0.079                4.25
# 7     7   9019 118.77  0.933   0.072   0.918     0.080   8011 115.14  0.945   0.065   0.923     0.109                3.62
# 8     8   7691 119.06  0.964   0.062   0.938     0.193   9039 115.24  0.971   0.066   0.937     0.227                3.82
# 9     9   8515 119.16  0.945   0.072   0.918     0.115   8013 115.25  0.979   0.058   0.944     0.278                3.90
# 10   10   9037 119.28  0.931   0.095   0.913     0.077   9037 115.83  0.967   0.068   0.935     0.207                3.45
# > cat("\n========================= TOP 10 SHORT SWEEP ========================\n")
# 
# ========================= TOP 10 SHORT SWEEP ========================
#   > print(as.data.frame(res_short$df))
# Rank 1A Mod 1A AIC 1A CFI 1A SRMR 1A AGFI 1A p_Chi2 1B Mod 1B AIC 1B CFI 1B SRMR 1B AGFI 1B p_Chi2 Delta AIC (1A - 1B)
# 1     1  30693  56.01  0.928   0.025   0.956     0.060  30713  56.53  1.000   0.049   0.982     0.442               -0.53
# 2     2  33509  57.53  1.000   0.025   0.976     0.473  30697  57.09  0.988   0.050   0.975     0.338                0.44
# 3     3  31253  58.31  0.923   0.027   0.954     0.054  31229  58.25  0.952   0.061   0.964     0.136                0.06
# 4     4  32885  58.91  0.967   0.032   0.965     0.208  30673  58.49  0.956   0.052   0.972     0.152                0.42
# 5     5  35189  59.13  0.987   0.039   0.973     0.331  31253  58.50  0.921   0.062   0.967     0.051                0.64
# 6     6  34069  59.57  0.962   0.045   0.962     0.184  33509  58.99  0.988   0.060   0.971     0.344                0.59
# 7     7  30709  59.62  0.978   0.026   0.966     0.271    473  60.20  0.999   0.065   0.981     0.432               -0.58
# 8     8  34629  59.81  0.946   0.042   0.967     0.117  32885  60.51  0.950   0.056   0.961     0.134               -0.70
# 9     9  32389  59.84  0.972   0.031   0.969     0.239  35193  60.61  0.986   0.065   0.972     0.330               -0.77
# 10   10  32949  59.94  0.988   0.026   0.972     0.342  30709  60.62  0.967   0.059   0.962     0.210               -0.68


#==================================
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

out_dir     <- "Rcode_for_paper/Routput_for_paper"
tbl_out_dir <- file.path(out_dir, "tables")
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
# Total N models sweeps
n_1a_long  <- nrow(res_1a_long_raw)
n_1b_long  <- nrow(res_1b_long_raw)
n_1a_short <- nrow(res_1a_short_raw)
n_1b_short <- nrow(res_1b_short_raw)


# Filter sweeps
valid_1a_long  <- filter_valid_sem(res_1a_long_raw)
valid_1b_long  <- filter_valid_sem(res_1b_long_raw)
valid_1a_short <- filter_valid_sem(res_1a_short_raw)
valid_1b_short <- filter_valid_sem(res_1b_short_raw)

cat(sprintf("LONG Valid Models:  DAG 1A = %d | DAG 1B = %d\n", nrow(valid_1a_long), nrow(valid_1b_long)))
cat(sprintf("SHORT Valid Models: DAG 1A = %d | DAG 1B = %d\n", nrow(valid_1a_short), nrow(valid_1b_short)))

# ------------------------------------------------------------------------------
# 3. Helper Function to Build and Format gt Tables
# ------------------------------------------------------------------------------
generate_comparison_table <- function(valid_a, valid_b, sweep_label, top_n = 10) {
  
  top_a <- valid_a %>% slice_head(n = top_n)
  top_b <- valid_b %>% slice_head(n = top_n)
  
  comp_df <- tibble(
    Rank                  = 1:top_n,
    `1A Mod`              = top_a$modNum,
    `1A AIC`              = round(top_a$AIC, 2),
    `1A CFI`              = round(top_a$CFI, 3),
    `1A SRMR`             = round(top_a$SRMR, 3),
    `1A AGFI`             = round(top_a$AGFI, 3),
    `1A p_Chi2`           = round(top_a$p_Chi2, 3),
    `1B Mod`              = top_b$modNum,
    `1B AIC`              = round(top_b$AIC, 2),
    `1B CFI`              = round(top_b$CFI, 3),
    `1B SRMR`             = round(top_b$SRMR, 3),
    `1B AGFI`             = round(top_b$AGFI, 3),
    `1B p_Chi2`           = round(top_b$p_Chi2, 3),
    `Delta AIC (1A - 1B)` = round(top_a$AIC - top_b$AIC, 2)
  )
  
  gt_tbl <- comp_df %>%
    gt() %>%
    tab_header(
      title = md(sprintf("**Comparison of Top 10 Valid Models: DAG 1A vs. DAG 1B (%s Sweep)**", sweep_label)),
      subtitle = md("*Filtered by liberal fit criteria: p(&chi;<sup>2</sup>) &ge; 0.05, CFI &ge; 0.80, SRMR &le; 0.10, AGFI &ge; 0.80*")
    ) %>%
    tab_spanner(label = "DAG 1A (Main Line)", columns = starts_with("1A")) %>%
    tab_spanner(label = "DAG 1B (Converging)", columns = starts_with("1B")) %>%
    cols_label(
      Rank                  = md("**Rank**"),
      `1A Mod`              = md("**Model**"),
      `1A AIC`              = md("**AIC**"),
      `1A CFI`              = md("**CFI**"),
      `1A SRMR`             = md("**SRMR**"),
      `1A AGFI`             = md("**AGFI**"),
      `1A p_Chi2`           = md("**p(&chi;<sup>2</sup>)**"),
      `1B Mod`              = md("**Model**"),
      `1B AIC`              = md("**AIC**"),
      `1B CFI`              = md("**CFI**"),
      `1B SRMR`             = md("**SRMR**"),
      `1B AGFI`             = md("**AGFI**"),
      `1B p_Chi2`           = md("**p(&chi;<sup>2</sup>)**"),
      `Delta AIC (1A - 1B)` = md("**&Delta;AIC<br>(1A &minus; 1B)**")
    ) %>%
    cols_align(align = "center", columns = everything()) %>%
    tab_source_note(
      source_note = md("Threshold references: p(&chi;<sup>2</sup>) (Bollen 1989); CFI (Marsh et al. 2004); SRMR (Hu & Bentler 1998; Kline 2016); AGFI (Cole 1987; Schermelleh-Engel et al. 2003).")
    ) %>%
    tab_options(
      table.font.size = px(9.5),
      heading.title.font.size = px(12),
      heading.subtitle.font.size = px(9.5),
      column_labels.font.weight = "bold",
      data_row.padding = px(3)
    )
  
  return(list(df = comp_df, gt = gt_tbl))
}

# ------------------------------------------------------------------------------
# 4. Generate Tables for Long & Short
# ------------------------------------------------------------------------------
res_long  <- generate_comparison_table(valid_1a_long,  valid_1b_long,  "Long / 24-Year")
res_short <- generate_comparison_table(valid_1a_short, valid_1b_short, "Short Sweep")

# Console display
cat("\n========================= TOP 10 LONG SWEEP =========================\n")
print(as.data.frame(res_long$df))

cat("\n========================= TOP 10 SHORT SWEEP ========================\n")
print(as.data.frame(res_short$df))

# ------------------------------------------------------------------------------
# 5. Export Tables (Word .docx, PNG, and HTML)
# ------------------------------------------------------------------------------
# Export Long Sweep Tables
gtsave(res_long$gt, file.path(tbl_out_dir, "Table_DAG1A_vs_DAG1B_Top10_Long.docx"))
gtsave(res_long$gt, file.path(tbl_out_dir, "Table_DAG1A_vs_DAG1B_Top10_Long.png"), vwidth = 1100, vheight = 650)
gtsave(res_long$gt, file.path(tbl_out_dir, "Table_DAG1A_vs_DAG1B_Top10_Long.html"))

# Export Short Sweep Tables
gtsave(res_short$gt, file.path(tbl_out_dir, "Table_DAG1A_vs_DAG1B_Top10_Short.docx"))
gtsave(res_short$gt, file.path(tbl_out_dir, "Table_DAG1A_vs_DAG1B_Top10_Short.png"), vwidth = 1100, vheight = 650)
gtsave(res_short$gt, file.path(tbl_out_dir, "Table_DAG1A_vs_DAG1B_Top10_Short.html"))

cat("\nAll Long & Short comparison tables successfully saved to: ", tbl_out_dir, "\n")
