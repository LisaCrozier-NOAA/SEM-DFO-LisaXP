
#Conclusion:
#AIC penalty for missing years: 1-3 (1 for 3-4 param to 3 for 1 param models)
#see 

library(tidyverse)
library(lavaan)

# --- 1. Preparation Function to Parse Formulas ---
# This helper automatically converts text strings like "X16_SAR ~ var1 + var2" into lavaan code
run_sim_model <- function(formula_str, model_data) {
  fit <- tryCatch({
    sem(formula_str, data = model_data, missing = "ML")
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(NA)
  return(fitMeasures(fit)[["aic"]])
}

# Ensure your dataset is sorted chronologically by your time variable (e.g., year/date)
guild_sorted <- guild.dfas1 %>% arrange(year)

# ==============================================================================
# EXPERIMENT 1: 24 Complete Years down to 19 Years------
# ==============================================================================

# 1. Identify the top 10 complete models for each parameter count (1 to 4)
top_24yr_models <- dredge_results %>%
  filter(Complete_Yrs == 24) %>%
  group_by(Num_Params) %>%
  slice_min(order_by = Raw_AIC, n = 10) %>%
  ungroup() %>%
  mutate(Model_ID = row_number())

sim_results_24to19 <- list()

# 2. Run the chronological truncation loop
for (drop_count in 0:5) {
  current_complete_yrs <- 24 - drop_count
  
  # Systematically slice rows from the beginning of the time series forward
  truncated_data <- guild_sorted %>% slice((drop_count + 1):n())
  
  for (i in 1:nrow(top_24yr_models)) {
    f_text <- top_24yr_models$Formula[i]
    k_val  <- top_24yr_models$Num_Params[i]
    m_id   <- top_24yr_models$Model_ID[i]
    
    # Calculate baseline absolute AIC score
    base_aic <- top_24yr_models$Raw_AIC[i]
    
    # Recalculate model fit on truncated timeline
    new_aic <- run_sim_model(f_text, truncated_data)
    
    sim_results_24to19[[length(sim_results_24to19) + 1]] <- tibble(
      Experiment = "24 to 19",
      Model_ID = m_id,
      Num_Params = k_val,
      Target_Yrs = current_complete_yrs,
      Yrs_Removed = drop_count,
      Baseline_AIC = base_aic,
      Recalculated_AIC = new_aic,
      # Track how much absolute AIC dropped due to the missing year
      AIC_Deflation = base_aic - new_aic 
    )
  }
}

df_sim_24to19 <- bind_rows(sim_results_24to19)

# ==============================================================================
# EXPERIMENT 2: 19 Inherently Missing Years down to 14 Years-------
# ==============================================================================


top_19yr_models <- dredge_results %>%
  filter(Complete_Yrs == 19) %>%
  group_by(Num_Params) %>%
  slice_min(order_by = Raw_AIC, n = 10) %>%
  ungroup() %>%
  mutate(Model_ID = row_number())

sim_results_19to14 <- list()

for (drop_count in 0:5) {
  current_complete_yrs <- 19 - drop_count
  
  for (i in 1:nrow(top_19yr_models)) {
    f_text <- top_19yr_models$Formula[i]
    k_val  <- top_19yr_models$Num_Params[i]
    m_id   <- top_19yr_models$Model_ID[i]
    base_aic <- top_19yr_models$Raw_AIC[i]
    
    # 1. Identify which variables are in this specific formula
    vars_in_model <- lavaanify(f_text) %>% 
      filter(op == "~") %>% 
      select(lhs, rhs) %>% 
      pivot_longer(cols = everything()) %>% 
      distinct() %>% pull(value)
    
    # 2. Isolate rows where this specific model actually has data
    model_specific_data <- guild_sorted %>% 
      filter(complete.cases(guild_sorted[, vars_in_model]))
    
    # 3. Drop years from the beginning of THIS specific model's active timeline
    truncated_data <- model_specific_data %>% slice((drop_count + 1):n())
    
    # 4. Re-fit on the truly nested subset
    new_aic  <- run_sim_model(f_text, truncated_data)
    
    sim_results_19to14[[length(sim_results_19to14) + 1]] <- tibble(
      Experiment = "19 to 14",
      Model_ID = m_id,
      Num_Params = k_val,
      Target_Yrs = current_complete_yrs,
      Yrs_Removed = drop_count,
      Baseline_AIC = base_aic,
      Recalculated_AIC = new_aic,
      AIC_Deflation = base_aic - new_aic
    )
  }
}


df_sim_19to14 <- bind_rows(sim_results_19to14)
tail(df_sim_19to14)

# Combine both simulation experiments into a master output dataset--------
master_simulation_profiles <- bind_rows(df_sim_24to19, df_sim_19to14)

head(master_simulation_profiles)
nrow(master_simulation_profiles) #450 = 10models * 5yrs_rm * 4param
summary(master_simulation_profiles)

master_simulation_profiles2<- master_simulation_profiles %>% filter(!Yrs_Removed==0) %>%
  mutate(daic=AIC_Deflation/Yrs_Removed)
tail(master_simulation_profiles2)


#lm model fit----------
fit<-lm(daic~Experiment + Num_Params, data=master_simulation_profiles2);summary(fit)
boxplot(daic~Experiment + Num_Params, data=master_simulation_profiles2)
fit2<-lm(daic~Experiment * Num_Params, data=master_simulation_profiles2);summary(fit2)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                    4.32281    0.21853  19.781  < 2e-16 ***
#   Experiment24 to 19            -0.71239    0.28165  -2.529 0.011839 *  
#   Num_Params                    -0.98301    0.07527 -13.059  < 2e-16 ***
#   Experiment24 to 19:Num_Params  0.36797    0.09937   3.703 0.000245 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.026 on 371 degrees of freedom
# Multiple R-squared:  0.4265,	Adjusted R-squared:  0.4219 
# F-statistic: 91.97 on 3 and 371 DF,  p-value: < 2.2e-16


# Ensure Num_Params is treated as a discrete group (factor) for proper spacing
master_simulation_profiles2 <- master_simulation_profiles2 %>%
  mutate(Num_Params = factor(Num_Params))



#Grouped Boxplot Code with Grand Mean Labels--------
#This script dynamically calculates the overall mean for each parameter number (pooling both experiments together), builds a clean label string, and positions it directly on the chart.
library(ggplot2)
library(dplyr)

# 1. Calculate the Grand Means (pooling Experiments) to use as labels
grand_means <- master_simulation_profiles2 %>%
  group_by(Num_Params) %>%
  summarize(
    Mean_Val = mean(daic, na.rm = TRUE),
    # Construct a clean text label
    Label_Text = paste0("Grand Mean:\n", round(Mean_Val, 2))
  ) %>%
  mutate(
    # Ensure Num_Params matches your factor structure
    Num_Params = factor(Num_Params, levels = c("1", "2", "3", "4"))
  )

# 2. Prep main plotting data structure
plot_data_force <- master_simulation_profiles2 %>%
  mutate(
    Num_Params = factor(Num_Params, levels = c("1", "2", "3", "4")),
    Experiment = factor(Experiment, levels = c("24 to 19", "19 to 14"))
  )

# 3. Create the Visualization

boxplot_aic_penalty_per_param_per_yr<-
  ggplot(plot_data_force, aes(x = Num_Params, y = daic)) +
  # Map 'fill' inside the boxplot specifically to keep the experiments side-by-side
  geom_boxplot(aes(fill = Experiment), outlier.color = "red", outlier.size = 2, alpha = 0.75, width = 0.6) +
  
  # Add a visual baseline reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  
  # OVERLAY TEXT: Explicitly pass the separate grand_means summary data
  geom_text(
    data = grand_means,
    #    aes(x = Num_Params, y = Mean_Val, label = Label_Text),
    aes(x = Num_Params, y = -1.5, label = Label_Text),
    # Place text slightly above the calculated mean value to clear the box centers
    vjust = -1.2, 
    color = "black",
    fontface = "bold",
    size = 3.8
  ) +
  
  # Ensure the x-axis scale preserves all parameter slots
  scale_x_discrete(drop = FALSE) +
  
  scale_fill_brewer(palette = "Set1", name = "Timeline Experiment") +
  theme_bw(base_size = 13) +
  labs(
    title = "Empirical Distribution of Annual AIC Deflation",
    subtitle = "Text callouts represent the combined Grand Mean for each parameter pool",
    x = "Number of Parameters in Model (k)",
    y = "Annual AIC Deflation (ΔAIC per year removed)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(boxplot_aic_penalty_per_param_per_yr,file="LisaXP/outputs_4/boxplot_aic_penalty_per_param_per_yr.png")


library(dplyr)
library(tidyr)

#annual_shift_summary table---------
annual_shift_summary <- master_simulation_profiles %>%
  # 1. Group by Experiment and Model to calculate step-by-step changes chronologically
  group_by(Experiment, Model_ID) %>%
  arrange(desc(Target_Yrs), .by_groups = TRUE) %>% # Order from highest years to lowest
  
  # 2. Calculate the exact difference in AIC from the previous step
  # (e.g., AIC at 23 years minus AIC at 24 years)
  mutate(Annual_AIC_Drop = Recalculated_AIC - lag(Recalculated_AIC)) %>%
  ungroup() %>%
  
  # 3. Filter out the very first year (where lag is NA) to isolate the active drops
  filter(!is.na(Annual_AIC_Drop)) %>%
  
  # 4. Group by your key structural criteria to see if parameter count matters
  group_by(Experiment, Num_Params) %>%
  summarize(
    Mean_Annual_Drop = round(mean(Annual_AIC_Drop, na.rm = TRUE), 3),
    SD_Annual_Drop   = round(sd(Annual_AIC_Drop, na.rm = TRUE), 3),
    Total_Steps_Run  = n(),
    .groups = "drop"
  )

# Print the clean summary matrix
print(annual_shift_summary)


#   Experiment Num_Params Mean_Annual_Drop SD_Annual_Drop Total_Steps_Run
# 1 19 to 14            1           -3.39            1.91              25
# 2 19 to 14            2           -2.56            1.76              50
# 3 19 to 14            3           -1.55            2.34              50
# 4 19 to 14            4           -1.36            2.20              50
# 5 24 to 19            1           -3.03            2.28              50
# 6 24 to 19            2           -3.17            3.96              50
# 7 24 to 19            3           -2.32            2.39              50
# 8 24 to 19            4           -0.999           1.44              50

#SAVE rdata file-----
save(dredge_results,
     top_24yr_models,
     top_19yr_models,
     annual_shift_summary,
     master_simulation_profiles,
     master_simulation_profiles2,
     boxplot_aic_penalty_per_param_per_yr,
     file="LisaXP/outputs_4/aicpenalty.rdata")
