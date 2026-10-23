# ------------------------------------------------------------------------------
# STEP 3 (FIXED): RENDER GT PUBLICATION SUMMARY TABLE
# ------------------------------------------------------------------------------
gt_comp_df <- comparison_table %>%
  select(
    Indicator,
    Baseline_Status,
    Baseline_Beta,
    Selected_Gates,
    MinMax_Beta_Peff,
    MinMax_AIC,
    Delta_AIC 
  )

gt_comp <- gt(gt_comp_df) %>%
  tab_header(
    title = md("**Comparative Diagnostic: Additive Table 4 SEMs vs. Min-Max Gated Models**"),
    subtitle = "Evaluating whether spatial overlap / alternate prey gating resolves baseline flipped predator signs."
  ) %>%
  tab_spanner(
    label = md("**Previous Baseline (Table 4 Additive)**"),
    columns = c(Baseline_Status, Baseline_Beta)
  ) %>%
  tab_spanner(
    label = md("**New Min-Max Attenuation Model**"),
    columns = c(Selected_Gates, MinMax_Beta_Peff, MinMax_AIC, Delta_AIC)
  ) %>%
  cols_label(
    Indicator        = md("**Predator Indicator**"),
    Baseline_Status  = md("**Status**"),
    Baseline_Beta    = md("**Base β**"),
    Selected_Gates   = md("**Selected Gate(s)**"),
    MinMax_Beta_Peff = md("**Gated β**"),
    MinMax_AIC       = md("**AIC**"),
    Delta_AIC        = md("**ΔAIC**")
  ) %>%
  # --- ROUND ALL NUMERIC COLUMNS TO 2 DECIMAL PLACES ---
  fmt_number(
    columns = where(is.numeric),
    decimals = 2
  ) %>%
  cols_align(align = "left", columns = c(Indicator, Selected_Gates)) %>%
  cols_align(align = "center", columns = c(Baseline_Status, Baseline_Beta, MinMax_Beta_Peff, MinMax_AIC, Delta_AIC)) %>%
  sub_missing(columns = everything(), missing_text = "-") %>%
  tab_style(
    style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#2E7D32", weight = "bold")),
    locations = cells_body(
      columns = c(MinMax_Beta_Peff, Selected_Gates),
      rows = Baseline_Status == "Flipped (+)" & MinMax_Beta_Peff < 0
    )
  ) %>%
  tab_options(
    table.font.size = px(10),
    heading.title.font.size = px(12),
    column_labels.font.weight = "bold",
    data_row.padding = px(3)
  )

print(gt_comp)

gtsave(gt_comp, file.path(tbl_out_dir, "Table4_vs_MinMax_Predator_Comparison.html"))
message("Comparison complete! Results written to: ", tbl_out_dir)

library(webshot2) # Required by gt to render PNGs

# Save as high-res PNG image (adjust vwidth/vheight as needed)
gtsave(
  data = gt_comp,
  filename = file.path(tbl_out_dir, "Table4_vs_MinMax_Predator_Comparison.png"),
  vwidth = 1000,
  vheight = 800,
  zoom = 2 # Increases image crispness for slides
)
