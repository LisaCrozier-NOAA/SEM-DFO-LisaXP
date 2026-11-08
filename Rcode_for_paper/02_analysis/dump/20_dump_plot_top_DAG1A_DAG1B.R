library(tidyverse)
library(lubridate)
library(lavaan)
library(tidySEM)
library(ggplot2)

# ==============================================================================
# 1. READ ROW 2 FOR DAG1A AND DAG1B (LONG)
# ==============================================================================
sweep_paths <- list(
  "DAG1A_long" = "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1A_long/SEMresultsByClus.csv",
  "DAG1B_long" = "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1B_long/SEMresultsByClus.csv"
)

step1_second_models <- imap_dfr(sweep_paths, function(path, name) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df %>%
    slice(2) %>% # Select second best model
    mutate(model_case = name) %>%
    select(model_case, modNum, AIC, p_Chi2, CFI, SRMR, RMSEA,
           PreyNCCindNames, PredNCCindNames, GrowthIndNames, 
           AbundanceIndNames, PreyAKindNames, PredAKindNames)
})

cat("--- 2nd Top Models for DAG 1A & 1B (Long) ---\n")
print(t(step1_second_models))


# ==============================================================================
# 2. LOAD DATA & METADATA
# ==============================================================================
clus_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DFA/clusDataDFA.rds"
cw_path   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/crosswalk_names.csv"

clusData <- readRDS(clus_path)
cw       <- read.csv(cw_path, stringsAsFactors = FALSE)

clusData_sub <- clusData %>%
  filter(minYear == 1998) %>%
  mutate(year = year(ymd(date))) %>%
  filter(year >= 1998, year <= 2021)

# Helpers
pull_std_node <- function(cdat, var_node, ind_name) {
  df_node <- cdat %>%
    filter(SEMlatent == var_node, shortName == ind_name) %>%
    select(date, finalVal) %>%
    arrange(date)
  df_node[[var_node]] <- as.vector(scale(df_node$finalVal))
  df_node %>% select(date, all_of(var_node))
}

pull_std_sar <- function(cdat) {
  cdat %>%
    filter(SEMlatent == "SAR") %>%
    select(date, finalVal) %>%
    arrange(date) %>%
    mutate(logSAR = as.vector(scale(finalVal))) %>%
    select(date, logSAR)
}

build_explicit_label <- function(node_title, ind_key, cw_df) {
  row_match <- cw_df %>%
    filter(DFAname == ind_key | dfa_cols == ind_key)
  lisa_val <- if (nrow(row_match) > 0 && !is.na(row_match$LisaName[1]) && row_match$LisaName[1] != "") {
    row_match$LisaName[1]
  } else {
    ind_key
  }
  paste0(node_title, "\n[", lisa_val, "]")
}

# ==============================================================================
# 3. FIT & PLOT: 2ND TOP MODEL DAG 1A LONG (24yr)
# ==============================================================================
row_1a_2nd <- step1_second_models %>% filter(model_case == "DAG1A_long")

dag1a_2nd_inds <- c(
  "PreyNCC"   = row_1a_2nd$PreyNCCindNames[1],
  "PredNCC"   = row_1a_2nd$PredNCCindNames[1],
  "Growth"    = row_1a_2nd$GrowthIndNames[1],
  "Abundance" = row_1a_2nd$AbundanceIndNames[1],
  "PreyAK"    = row_1a_2nd$PreyAKindNames[1],
  "PredAK"    = row_1a_2nd$PredAKindNames[1]
)

df_dag1a_2nd <- pull_std_node(clusData_sub, "PreyNCC",   dag1a_2nd_inds["PreyNCC"]) %>%
  full_join(pull_std_node(clusData_sub, "PredNCC",   dag1a_2nd_inds["PredNCC"]),   by = "date") %>%
  full_join(pull_std_node(clusData_sub, "Growth",    dag1a_2nd_inds["Growth"]),    by = "date") %>%
  full_join(pull_std_node(clusData_sub, "Abundance", dag1a_2nd_inds["Abundance"]), by = "date") %>%
  full_join(pull_std_node(clusData_sub, "PreyAK",    dag1a_2nd_inds["PreyAK"]),    by = "date") %>%
  full_join(pull_std_node(clusData_sub, "PredAK",    dag1a_2nd_inds["PredAK"]),    by = "date") %>%
  full_join(pull_std_sar(clusData_sub), by = "date") %>%
  filter(!is.na(logSAR)) %>%
  arrange(date)

syntax_1a <- '
  Growth    ~ PreyNCC
  Abundance ~ Growth + PredNCC
  logSAR    ~ Abundance + PreyAK + PredAK
'
fit_1a_2nd <- sem(syntax_1a, data = df_dag1a_2nd)

cat("\n--- DAG 1A (2nd Top Model) Parameter Estimates ---\n")
parameterEstimates(fit_1a_2nd, standardized = TRUE, ci = TRUE) %>%
  filter(op == "~") %>%
  select(lhs, op, rhs, est, se, z, pvalue, std.all, ci.lower, ci.upper) %>%
  print()

labels_1a_2nd <- c(
  "PreyNCC"   = build_explicit_label("Prey (NCC)", dag1a_2nd_inds["PreyNCC"], cw),
  "PredNCC"   = build_explicit_label("Predators (NCC)", dag1a_2nd_inds["PredNCC"], cw),
  "Growth"    = build_explicit_label("Salmon Growth", dag1a_2nd_inds["Growth"], cw),
  "Abundance" = build_explicit_label("Salmon Abund.", dag1a_2nd_inds["Abundance"], cw),
  "PreyAK"    = build_explicit_label("Prey (AK)", dag1a_2nd_inds["PreyAK"], cw),
  "PredAK"    = build_explicit_label("Predators (AK)", dag1a_2nd_inds["PredAK"], cw),
  "logSAR"    = "Salmon SAR\n[X16_SAR]"
)

layout_1a <- get_layout(
  "PreyNCC", "", "",          "", "",       "", "PreyAK",
  "",        "", "Growth",    "", "",       "", "",
  "",        "", "Abundance", "", "logSAR", "", "",
  "PredNCC", "", "",          "", "",       "", "PredAK",
  rows = 4
)

p_1a_2nd <- Lisa_sem_graph.v3(
  fit               = fit_1a_2nd,
  title             = paste0("DAG 1A Long (24yr) - 2nd Top Model (Mod #", row_1a_2nd$modNum, ", AIC = ", round(row_1a_2nd$AIC, 2), ")"),
  layout            = layout_1a,
  node_labels       = labels_1a_2nd,
  r2_nodes          = c("Growth", "Abundance", "logSAR"),
  p_threshold       = 0.05,
  nonsig_linetype   = "dashed",
  font_family       = "serif",
  rect_width        = 2.4, rect_height = 0.9,
  ellipses_width    = 2.4, ellipses_height = 0.9,
  text_size         = 3.2, variance_diameter = 0.6, curvature = 10,
  width             = 14, height = 8
)

# ==============================================================================
# 4. FIT & PLOT: 2ND TOP MODEL DAG 1B LONG (24yr)
# ==============================================================================
row_1b_2nd <- step1_second_models %>% filter(model_case == "DAG1B_long")

dag1b_2nd_inds <- c(
  "PreyNCC"   = row_1b_2nd$PreyNCCindNames[1],
  "PredNCC"   = row_1b_2nd$PredNCCindNames[1],
  "Growth"    = row_1b_2nd$GrowthIndNames[1],
  "Abundance" = row_1b_2nd$AbundanceIndNames[1],
  "PreyAK"    = row_1b_2nd$PreyAKindNames[1],
  "PredAK"    = row_1b_2nd$PredAKindNames[1]
)

df_dag1b_2nd <- pull_std_node(clusData_sub, "PreyNCC",   dag1b_2nd_inds["PreyNCC"]) %>%
  full_join(pull_std_node(clusData_sub, "PredNCC",   dag1b_2nd_inds["PredNCC"]),   by = "date") %>%
  full_join(pull_std_node(clusData_sub, "Growth",    dag1b_2nd_inds["Growth"]),    by = "date") %>%
  full_join(pull_std_node(clusData_sub, "Abundance", dag1b_2nd_inds["Abundance"]), by = "date") %>%
  full_join(pull_std_node(clusData_sub, "PreyAK",    dag1b_2nd_inds["PreyAK"]),    by = "date") %>%
  full_join(pull_std_node(clusData_sub, "PredAK",    dag1b_2nd_inds["PredAK"]),    by = "date") %>%
  full_join(pull_std_sar(clusData_sub), by = "date") %>%
  filter(!is.na(logSAR)) %>%
  arrange(date)

# Symmetrical DAG 1B converging structure: Growth and Abundance both connect directly to logSAR[cite: 1]
syntax_1b <- '
  Growth    ~ PreyNCC
  Abundance ~ PredNCC
  logSAR    ~ Growth + Abundance + PreyAK + PredAK
'
fit_1b_2nd <- sem(syntax_1b, data = df_dag1b_2nd)

cat("\n--- DAG 1B (2nd Top Model) Parameter Estimates ---\n")
parameterEstimates(fit_1b_2nd, standardized = TRUE, ci = TRUE) %>%
  filter(op == "~") %>%
  select(lhs, op, rhs, est, se, z, pvalue, std.all, ci.lower, ci.upper) %>%
  print()

labels_1b_2nd <- c(
  "PreyNCC"   = build_explicit_label("Prey (NCC)", dag1b_2nd_inds["PreyNCC"], cw),
  "PredNCC"   = build_explicit_label("Predators (NCC)", dag1b_2nd_inds["PredNCC"], cw),
  "Growth"    = build_explicit_label("Salmon Growth", dag1b_2nd_inds["Growth"], cw),
  "Abundance" = build_explicit_label("Salmon Abund.", dag1b_2nd_inds["Abundance"], cw),
  "PreyAK"    = build_explicit_label("Prey (AK)", dag1b_2nd_inds["PreyAK"], cw),
  "PredAK"    = build_explicit_label("Predators (AK)", dag1b_2nd_inds["PredAK"], cw),
  "logSAR"    = "Salmon SAR\n[X16_SAR]"
)

layout_1b <- get_layout(
  "PreyNCC", "", "",          "", "",       "", "PreyAK",
  "",        "", "Growth",    "", "",       "", "",
  "",        "", "",          "", "logSAR", "", "",
  "",        "", "Abundance", "", "",       "", "",
  "PredNCC", "", "",          "", "",       "", "PredAK",
  rows = 5
)

p_1b_2nd <- Lisa_sem_graph.v3(
  fit               = fit_1b_2nd,
  title             = paste0("DAG 1B Long (24yr) - 2nd Top Model (Mod #", row_1b_2nd$modNum, ", AIC = ", round(row_1b_2nd$AIC, 2), ")"),
  layout            = layout_1b,
  node_labels       = labels_1b_2nd,
  r2_nodes          = c("Growth", "Abundance", "logSAR"),
  p_threshold       = 0.05,
  nonsig_linetype   = "dashed",
  font_family       = "serif",
  rect_width        = 2.4, rect_height = 0.9,
  ellipses_width    = 2.4, ellipses_height = 0.9,
  text_size         = 3.2, variance_diameter = 0.6, curvature = 10,
  width             = 14, height = 9
)
