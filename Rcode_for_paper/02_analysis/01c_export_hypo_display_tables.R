# ==============================================================================
# Script: 02_analysis/01c_export_hypo_display_tables.R
# Purpose: Export LaTeX-ready display CSVs for Main Table & Supplemental Table S2
#          from previously saved summary CSVs (no model re-running needed)
# ==============================================================================

library(tidyverse)

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
data_out_dir <- file.path(proj_dir, "Outputs", "data")

# ------------------------------------------------------------------------------
# PART A: MAIN-TEXT TABLE (Hypothesis Path Tests Across DAGs)
# ------------------------------------------------------------------------------
path_summary <- read_csv(
  file.path(data_out_dir, "obj1_datasummary_table2_hypo_path_tests_by_hypo_all_valid.csv"),
  show_col_types = FALSE
)

model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short")

t2_wide <- path_summary %>%
  mutate(
    # Format multi-line LaTeX cell with green/red highlighting
    Cell_LaTeX = case_when(
      Verdict == "Supported" ~ paste0(
        "\\textcolor[HTML]{1B5E20}{\\textbf{Supported}} \\newline ",
        Pct_Sig, "\\% Sig | ", Pct_Positive, "\\% + \\newline (N=", N_Models, ")"
      ),
      Verdict == "Flipped Sign" ~ paste0(
        "\\textcolor[HTML]{B71C1C}{\\textbf{Flipped Sign}} \\newline ",
        Pct_Sig, "\\% Sig | ", Pct_Positive, "\\% + \\newline (N=", N_Models, ")"
      ),
      TRUE ~ paste0(
        "\\textbf{NS} \\newline ",
        Pct_Sig, "\\% Sig | ", Pct_Positive, "\\% + \\newline (N=", N_Models, ")"
      )
    ),
    path_label = gsub("->", "$\\\\rightarrow$", path_label)
  ) %>%
  select(H_code, path_label, Baseline_Expectation, model_id, Cell_LaTeX) %>%
  pivot_wider(
    names_from  = model_id,
    values_from = Cell_LaTeX,
    values_fill = "Path Excluded"
  ) %>%
  mutate(
    h_rank   = match(H_code, c("H1a", "H1b", "H1c", "H1e", "H2a", "H2c", "H4a")),
    sub_rank = if_else(grepl("Competitor", Baseline_Expectation), 2, 1)
  ) %>%
  arrange(h_rank, sub_rank) %>%
  select(H_code, path_label, Baseline_Expectation, any_of(model_names))

# Blank out duplicate Code & Structural Path on the 2nd sub-row of split paths (H1a, H1e)
dup_rows <- which(t2_wide$H_code == lag(t2_wide$H_code, default = ""))
t2_wide$H_code[dup_rows]     <- ""
t2_wide$path_label[dup_rows] <- ""

write_csv(t2_wide, file.path(data_out_dir, "obj1_table2_display_main.csv"))

# ------------------------------------------------------------------------------
# PART B: SUPPLEMENTARY TABLE S2 (Individual Indicator Diagnostics)
# ------------------------------------------------------------------------------
tableS2_data <- read_csv(
  file.path(data_out_dir, "obj1_data_table2_hypo_path_tests_by_indicator.csv"),
  show_col_types = FALSE
)

s2_formatted <- tableS2_data %>%
  arrange(H_code, SEM_Node, Indicator) %>%
  mutate(
    Indicator_Clean = gsub("_", " ", Indicator),
    N_Str           = as.character(N_Models),
    Pct_Supp_Str    = if_else(is.na(Pct_Sig_Supp), "-", sprintf("%.1f", Pct_Sig_Supp)),
    Beta_Supp_Str   = if_else(is.na(Beta_Supp),    "-", sprintf("%.3f", Beta_Supp)),
    Pct_Flip_Str    = if_else(is.na(Pct_Sig_Flip), "-", sprintf("%.1f", Pct_Sig_Flip)),
    Beta_Flip_Str   = if_else(is.na(Beta_Flip),    "-", sprintf("%.3f", Beta_Flip)),
    
    # Apply green (>=50% supported sig) and red (>=50% flipped sig) text styling
    Pct_Supp_Disp = if_else(!is.na(Pct_Sig_Supp) & Pct_Sig_Supp >= 50,
                            paste0("\\textcolor[HTML]{1B5E20}{\\textbf{", Pct_Supp_Str, "}}"), Pct_Supp_Str),
    Beta_Supp_Disp = if_else(!is.na(Pct_Sig_Supp) & Pct_Sig_Supp >= 50,
                             paste0("\\textcolor[HTML]{1B5E20}{\\textbf{", Beta_Supp_Str, "}}"), Beta_Supp_Str),
    Pct_Flip_Disp = if_else(!is.na(Pct_Sig_Flip) & Pct_Sig_Flip >= 50,
                            paste0("\\textcolor[HTML]{B71C1C}{\\textbf{", Pct_Flip_Str, "}}"), Pct_Flip_Str),
    Beta_Flip_Disp = if_else(!is.na(Pct_Sig_Flip) & Pct_Sig_Flip >= 50,
                             paste0("\\textcolor[HTML]{B71C1C}{\\textbf{", Beta_Flip_Str, "}}"), Beta_Flip_Str)
  )

# Insert full-width bold H_code subheader rows
h_codes <- unique(s2_formatted$H_code)
s2_grouped <- map_dfr(h_codes, function(hc) {
  p_lbl <- gsub("->", "$\\\\rightarrow$", s2_formatted$path_label[s2_formatted$H_code == hc][1])
  header_row <- tibble(
    Indicator_Clean = paste0("\\multicolumn{6}{l}{\\textbf{", hc, " (", p_lbl, ")}}"),
    N_Str = "", Pct_Supp_Disp = "", Beta_Supp_Disp = "", Pct_Flip_Disp = "", Beta_Flip_Disp = ""
  )
  body_rows <- s2_formatted %>%
    filter(H_code == hc) %>%
    select(Indicator_Clean, N_Str, Pct_Supp_Disp, Beta_Supp_Disp, Pct_Flip_Disp, Beta_Flip_Disp)
  bind_rows(header_row, body_rows)
})

write_csv(s2_grouped, file.path(data_out_dir, "obj1_tableS2_display_supp.csv"))
cat("Saved obj1_table2_display_main.csv and obj1_tableS2_display_supp.csv!\n")
