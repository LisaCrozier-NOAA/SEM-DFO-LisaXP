
#Table 2 Template (Methods): A Priori Ecological Hypotheses
#Defines directional expectations for path coefficients across all DAG model configurations.

library(gt)
library(tibble)

hypotheses_df <- tribble(
  ~Hypothesis, ~Mechanism, ~Target_Path, ~Expected_Sign, ~Primary_Reference,
  "H1 (Bottom-up)", "Early marine growth promotes juvenile survival", "NCC Growth -> NCC Survival", "+", "Crozier et al. (2021)",
  "H2 (Bottom-up)", "Broad marine productivity enhances overall return", "AK Productivity -> SAR", "+", "Wells et al. (2020)",
  "H3 (Top-down)", "Early ocean predator pressure reduces survival", "NCC Predators -> NCC Survival", "-", "Wells et al. (2025)",
  "H4 (Top-down)", "Late marine predation creates survival bottleneck", "AK Predators -> SAR", "-", "Crozier et al. (2025)",
  "H5 (Buffer Prey)", "Alternate prey buffers predator consumption of salmon", "Pred_AK * Prey_AK -> SAR", "+", "Wells et al. (2023)"
)

table2_gt <- hypotheses_df %>%
  gt() %>%
  tab_header(
    title = md("**Table 2. Proposed Ecological Hypotheses and Path Expectations**")
  ) %>%
  cols_label(
    Hypothesis = md("**Hypothesis**"),
    Mechanism = md("**Ecological Mechanism**"),
    Target_Path = md("**Target SEM Path**"),
    Expected_Sign = md("**Expected Sign**"),
    Primary_Reference = md("**Supporting Literature**")
  ) %>%
  cols_align(align = "center", columns = c(Expected_Sign)) %>%
  tab_options(table.font.size = px(12))

gtsave(table2_gt, "output/tables/Table2_Hypotheses.html")