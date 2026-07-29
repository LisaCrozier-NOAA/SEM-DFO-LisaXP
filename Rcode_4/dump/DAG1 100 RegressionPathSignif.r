library(tidyverse)
model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short", "DAG1C_long", "DAG1C_short")

# --- 2. Define the Paths You Want to Track ---
# You can add or modify any pathways here (e.g., Growth -> cpue/Abundance, Growth -> SAR)
target_paths <- tribble(
  ~lhs_var,      ~rhs_var,   ~path_label,
  "Abundance",   "Growth",   "Growth -> Abundance",
  "logSAR",     "Growth",    "Growth -> SAR",
  "logSAR",   "Abundance",   "Abundance -> SAR",
  "Abundance",   "PreyNCC",  "PreyNCC -> Abundance"
)

# --- 3. Read & Combine Model Rankings (To find the Top 100) ---
cat("Loading model rankings...\n")
rankings_all <- map_dfr(model_names, function(name) {
  file_path <- paste0(path, name, "/SEMresultsByClus.csv")

# --- 1. Define Paths & Runs ---
path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/"
  
  if (!file.exists(file_path)) {
    warning(paste("File missing:", file_path))
    return(tibble())
  }
  
  read_csv(file_path, show_col_types = FALSE) %>%
    mutate(model_id = name)
})

# Identify the top 100 models for each run based on AIC
top_100_keys <- rankings_all %>%
  filter(PreyNCCindNames != "01.ZooPreyNCC_JSOES_1_smoothed") %>%
  group_by(model_id) %>%
  # Sort by AIC (ascending) and grab the top 100
  slice_min(order_by = AIC, n = 100, with_ties = FALSE) %>%
#  slice_min(order_by = AIC, n = 25, with_ties = FALSE) %>%
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

# --- 5. Extract and Analyze Path Significance ---
# Filter estimates strictly to our target pathways
targeted_estimates <- estimates_all %>%
  inner_join(target_paths, by = c("lhs" = "lhs_var", "rhs" = "rhs_var")) %>%
  mutate(
    is_significant = if_else(pvalue < 0.05, 1, 0),
    # Add indicator flag if this model iteration is in the top 100 for its group
    is_top_100 = if_else(paste(model_id, modNum) %in% paste(top_100_keys$model_id, top_100_keys$modNum), TRUE, FALSE)
  )

# --- 6. Generate the Clean Summary Table ---
summary_table_noHCI <- targeted_estimates %>%
  #summary_table <- targeted_estimates %>%
  group_by(model_id, path_label) %>%
  summarise(
    # Overall Model Stats
    Total_Models_Evaluated = n(),
    Signif_Overall = sum(is_significant, na.rm = TRUE),
    Pct_Signif_Overall = round((Signif_Overall / Total_Models_Evaluated) * 100, 1),
    
    # Top 100 Model Stats
    Top_100_Evaluated = sum(is_top_100),
    Signif_In_Top_100 = sum(is_significant[is_top_100], na.rm = TRUE),
    Pct_Signif_In_Top_100 = round((Signif_In_Top_100 / Top_100_Evaluated) * 100, 1),
    
    .groups = "drop"
  ) %>%
  arrange(path_label, model_id)


summary_table_noHCI <- summary_table_noHCI %>%
  mutate(Hypothesis_tested=case_when(
    path_label=="Abundance -> SAR" ~ "Early marine survival predicts SAR",
    path_label=="Growth -> Abundance" ~ "Bigger is better",
    path_label=="Growth -> SAR" ~ "Bigger is better",
    path_label=="PreyNCC -> Abundance" ~ "Alternative prey for NCC predators improves salmon survival",
    TRUE ~ " "
    
  ))



# --- 7. Save and View Output ---
print(summary_table_noHCI, n = Inf)
write_csv(summary_table_noHCI, "LisaXP/outputs_4/MainRegression_Significance_Summary_noHCI.csv")
#summary_table_top25<-summary_table
#write_csv(summary_table_top25, paste0("LisaXP/outputs_4/MainRegression_Significance_Summary_top25.csv"))
