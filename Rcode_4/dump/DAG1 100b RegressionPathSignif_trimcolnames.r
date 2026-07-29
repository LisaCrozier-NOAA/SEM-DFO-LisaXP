library(tidyverse)

# --- 1. Define Paths & Runs ---
path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/"
model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short", "DAG1C_long", "DAG1C_short")

# --- 2. Define the Paths You Want to Track ---
target_paths <- tribble(
  ~lhs_var,      ~rhs_var,   ~path_label,
  "Abundance",   "Growth",   "Growth -> Abundance",
  "logSAR",      "Growth",   "Growth -> SAR",
  "logSAR",      "Abundance", "Abundance -> SAR",
  "logSAR",      "PredNCC", "PredNCC -> SAR",
  "Abundance",   "PreyNCC",  "PreyNCC -> Abundance"
)

# --- 3. Read & Combine Model Rankings (To find the Top 100, excluding HCI) ---
cat("Loading model rankings...\n")
rankings_all <- map_dfr(model_names, function(name) {
  file_path <- paste0(path, name, "/SEMresultsByClus.csv")
  
  if (!file.exists(file_path)) {
    warning(paste("File missing:", file_path))
    return(tibble())
  }
  
  read_csv(file_path, show_col_types = FALSE) %>%
    mutate(model_id = name)
})

# Identify the top 100 cleanest models (excluding Habitat Compression Index)
top_100_keys <- rankings_all %>%
  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>%
  group_by(model_id) %>%
  slice_min(order_by = AIC, n = 100, with_ties = FALSE) %>%
  select(model_id, modNum) %>%
  ungroup()

# --- 4. Read & Combine Parameter Estimates ---
cat("Loading parameter estimates...\n")
estimates_all <- map_dfr(model_names, function(name) {
  file_path <- paste0(path, name, "/parameterEstimates.csv")
  
  if (!file.exists(file_path)) {
    warning(paste("File missing:", file_path))
    return(tibble())
  }
  
  read_csv(file_path, show_col_types = FALSE) %>%
    mutate(model_id = name)
})

# --- 5. Extract and Analyze Path Significance & Direction ---
targeted_estimates <- estimates_all %>%
  inner_join(target_paths, by = c("lhs" = "lhs_var", "rhs" = "rhs_var")) %>%
  mutate(
    is_significant = if_else(pvalue < 0.05, 1, 0),
    is_positive    = if_else(est > 0, 1, 0),
    is_top_100     = if_else(paste(model_id, modNum) %in% paste(top_100_keys$model_id, top_100_keys$modNum), TRUE, FALSE)
  )

# --- 6. Generate the Clean Summary Table with Short Column Names ---
summary_table_noHCI <- targeted_estimates %>%
  group_by(model_id, path_label) %>%
  summarise(
    # Overall Model Stats
    N_All      = n(),
    Sig_All    = round((sum(is_significant, na.rm = TRUE) / N_All) * 100, 1),
    Pos_All    = round((sum(is_positive, na.rm = TRUE) / N_All) * 100, 1),
    Neg_All    = 100 - Pos_All,
    
    # Top 100 Model Stats (Omitting count column)
    Sig_100    = round((sum(is_significant[is_top_100], na.rm = TRUE) / sum(is_top_100)) * 100, 1),
    Pos_100    = round((sum(is_positive[is_top_100], na.rm = TRUE) / sum(is_top_100)) * 100, 1),
    Neg_100    = 100 - Pos_100,
    
    .groups = "drop"
  ) %>%
  arrange(path_label, model_id)

# --- 7. Add Ecological Hypotheses and Confirmations ---
summary_table_noHCI<-summary_table_noHCI %>% rename(model_id=Model, path_label=Path)
summary_table_noHCI <- summary_table_noHCI %>%
  mutate(
    Hypothesis = case_when(
      path_label == "Abundance -> SAR"     ~ "Early marine survival predicts SAR",
      path_label == "Growth -> Abundance"  ~ "Bigger is better",
      path_label == "Growth -> SAR"        ~ "Bigger is better",
      path_label == "PreyNCC -> Abundance" ~ "Alternative prey improves salmon survival",
      path_label == "PredNCC -> SAR" ~ "NCC predators reduce survival beyond JSOES cpue",
      TRUE                                 ~ " "
    ),
    Support = case_when(
      Pos_100 > 50 & Sig_100 > 50  ~ "Yes",
      Pos_100 > 50 & Sig_100 > 10 & Sig_100 <= 50 ~ "Weak Support (Low Sig)",
      Pos_100 > 50 & Sig_100 <= 10 ~ "No (Not Significant)",
      Pos_100 <= 50 & Sig_100 <= 10 ~ "No (NS & Flipped Sign)",
      Pos_100 <= 50 & Sig_100 > 10 ~ "Yes (but reduced main effect of abundance)",
      TRUE                         ~ " "
    )
  ) %>%
  # Reorganize for tight column packing
  select(
    Model = model_id,
    Path = path_label,
    N_All,
    Sig_All,
    Pos_All,
    Neg_All,
    Sig_100,
    Pos_100,
    Neg_100,
    Hypothesis,
    Support
  )

# --- 8. Save and View Output ---
print(as_tibble(summary_table_noHCI), n = Inf)
write_csv(summary_table_noHCI, "LisaXP/outputs_4/MainRegression_Significance_Summary_noHCI_shortnames.csv")



