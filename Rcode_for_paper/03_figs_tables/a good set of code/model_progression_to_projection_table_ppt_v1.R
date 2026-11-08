# ==============================================================================
# Script: 03_figs_tables/model_progression_to_projection_table_ppt_v1.R
# Purpose: Generate Slide-Ready Stepped Model Comparison Table
# ==============================================================================

library(tidyverse)
library(gt)

progression_data <- tribble(
  ~Step, ~Model_Description, ~Topology, ~SAR_R2, ~SAR_AIC, ~Mechanistic_Takeaway,
  
  "Baseline (No AK Preds)", 
  "Early Marine Only", 
  "CPUE ~ Hake<br>SAR ~ CPUE", 
  "~0.38", "-22.6", 
  "Early ocean abundance provides modest signal; late-marine mortality remains entirely unmodeled.",
  
  "Raw AK Baseline", 
  "Raw Predator Only", 
  "CPUE ~ Hake<br>SAR ~ CPUE + Pred_raw", 
  "0.47 – 0.48", "-25.6 to -26.4", 
  "Predator biomass alone lacks interannual explanatory power (often yields counter-intuitive positive signs).",
  
  "Single Arena (Ocean)", 
  "Thermal & Forage Gating", 
  "CPUE ~ Hake<br>SAR ~ CPUE + P_ocean_eff", 
  "0.53 – 0.56", "-28.8 to -29.8", 
  "Joint SST gating and alternate prey buffering consistently resolve sign flips and deliver significant negative impacts (p < 0.005).",
  
  "Dual Arena (Final A)", 
  "Sleeper Shark + Capelin + embSSL", 
  "CPUE ~ Hake<br>SAR ~ CPUE + P_ocean_eff + embSSL", 
  "<b>0.610</b>", "<b>-30.1</b>", 
  "<b>Top Model:</b> Joint subarctic marine arena (p = 0.0006) paired with freshwater adult return bottleneck (p = 0.12).",
  
  "Dual Arena (Final B)", 
  "Sablefish (3+) + Forage + embSSL", 
  "CPUE ~ Hake<br>SAR ~ CPUE + P_ocean_eff + embSSL", 
  "<b>0.577</b>", "<b>-28.6</b>", 
  "<b>Viable Parallel:</b> Demersal groundfish arena buffered by pollock & herring (p = 0.0021) plus adult return bottleneck."
)

final_progression_gt <- progression_data %>%
  gt() %>%
  tab_header(
    title = md("**Stepped Evolution of Chinook Marine Survival Models (24-Year Dataset)**"),
    subtitle = "Progressing from phenomenological baselines to a 3-stage mechanistic framework"
  ) %>%
  cols_label(
    Step                 = md("**Step**"),
    Model_Description    = md("**Hypothesis Tested**"),
    Topology             = md("**Model Structure**"),
    SAR_R2               = md("**SAR R²**"),
    SAR_AIC              = md("**AIC**"),
    Mechanistic_Takeaway = md("**Ecological Takeaway**")
  ) %>%
  fmt_markdown(columns = c(Topology, SAR_R2, SAR_AIC, Mechanistic_Takeaway)) %>%
  cols_align(align = "left", columns = c(Step, Model_Description, Topology, Mechanistic_Takeaway)) %>%
  cols_align(align = "center", columns = c(SAR_R2, SAR_AIC)) %>%
  cols_width(
    Step                 ~ px(130),
    Model_Description    ~ px(160),
    Topology             ~ px(200),
    SAR_R2               ~ px(85),
    SAR_AIC              ~ px(85),
    Mechanistic_Takeaway ~ px(300)
  ) %>%
  tab_options(
    table.width                = px(960),
    table.font.size            = px(13),
    heading.title.font.size    = px(18),
    heading.subtitle.font.size = px(12),
    column_labels.font.size    = px(14),
    column_labels.font.weight  = "bold",
    data_row.padding           = px(6)
  )

print(final_progression_gt)
gtsave(
  data     = final_progression_gt,
  filename = file.path(tbl_out_dir, "Goal5_Table_Verified_Progression_Models_simple.png"),
  vwidth   = 1000,
  vheight  = 600
)
