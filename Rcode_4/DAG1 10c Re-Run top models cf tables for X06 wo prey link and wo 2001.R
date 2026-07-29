#Conclusions:
# DFA_IGF_mu still chosen by B, even without preyNCC predicting it
# StomFull_May beats Lmu and IGF for A, even recognizing its advantage because it was missing 1998 -- it has a natural 2-4 adv, but the daic was 7
# Lmu_IntSprJunH won C by a landslide, even without 2001





#Now we have a harder question, which is whether the same top variables would have been chosen under slightly different conditions. 
#First I want to focus on the selection of the X06 node. 
#I think there might have been one very influential year, 2001, and I want to rerun the models without that year. 
#In the current comparison, X06_DFA_IGF_mu was chosen in 1A and 1B, all forms, while X06_Lmu_IntSprJunH and X06_StomFull_May were chosen in the 1C models. 
#There is an alternative explanation for this switch, which is that in A and B models, I forced prey into the X06 node, and there are other variables that are better at predicting X06_DFA_IGF_mu than X06_Lmu_IntSprJunH or X06_StomFull_May. 
#I want to test those 2 hypotheses. 
#Can you set up a code block that 
#1) removes the X06_DFA_IGF_mu ~ X05_DFA_abundSardine link in all of the ones where it existed, 
#2) anywhere else that any of the X06 variables appear, test all 3 of the options in that position, and 
#3) repeat that full test without 2001. I want to compare dAIC for all 3 sar response variables, like we did before.

#This is an excellent way to stress-test your ecological model. You are systematically checking if your structural pathways are driven by a single high-leverage year (2001) or if the forced link in your A and B frameworks artificially constrained which X06 proxy variable floated to the top.To test both of these hypotheses simultaneously, we need to:Drop the forced link: Remove X06_DFA_IGF_mu ~ X05_DFA_abundSardine entirely.Loosen the node slots: Identify where any X06 variable was acting as a predictor or target, and programmatically cycle all 3 options (X06_DFA_IGF_mu, X06_Lmu_IntSprJunH, X06_StomFull_May) through that position.Split the timeline: Run this full combinatorial test on your complete dataset, and then repeat it entirely on a subsetted dataset where the year 2001 is removed.Here is the complete script to generate, fit, and pivot the $\Delta\text{AIC}$ results for this analysis.


library(tidyverse)
library(lavaan)
library(glue)

# --- 1. Define Base Syntaxes (Forced X05 link removed, X06 replaced with a token) ---
# We use {X06_VAR} and {SAR_VAR} as placeholders to dynamically reconstruct models.
base_dags <- list(
  DAG1C_long = '
    X07_DFA_cpue_IntSprJunHW ~ X05_herring_GAM  + X09_DFA_HakeAge5Plus 
    {SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + {X06_VAR} + X05_herring_GAM + X09_DFA_HakeAge5Plus + X12_copepodCom_EGoA + X15_salmonSharkGoA_predAK
  ',
  DAG1C_long_reduced = '
    X07_DFA_cpue_IntSprJunHW ~  X09_DFA_HakeAge5Plus 
    {SAR_VAR} ~ {X06_VAR} + X05_herring_GAM + X09_DFA_HakeAge5Plus + X12_copepodCom_EGoA + X15_salmonSharkGoA_predAK
  ',
  DAG1C_short = '
    X07_DFA_cpue_IntSprJunHW ~ X01_habCompInd + X08_commonMurre_JSOES 
    {SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + {X06_VAR} + X01_habCompInd + X08_commonMurre_JSOES + X14_pinkSalmon + X10_Harbor_seal_CR_2yrLead
  ',
  DAG1C_short_reduced = '
    X07_DFA_cpue_IntSprJunHW ~  X08_commonMurre_JSOES 
    {SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + {X06_VAR} + X01_habCompInd + X08_commonMurre_JSOES + X14_pinkSalmon + X10_Harbor_seal_CR_2yrLead
  ',
  DAG1A_long = '
    X07_DFA_cpue_IntSprJunHW ~ {X06_VAR} + X09_DFA_HakeAge5Plus 
    {SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X12_DFA_biomassEuphShelfSum + X15_DFA_sleeperSharkBSAI_predAK  
  ',
  DAG1A_long_reduced = '
    X07_DFA_cpue_IntSprJunHW ~  X09_DFA_HakeAge5Plus 
    {SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X12_DFA_biomassEuphShelfSum + X15_DFA_sleeperSharkBSAI_predAK  
  ',
  DAG1A_short = '
    X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES + {X06_VAR}
    {SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X12_copepodBiomass_WGoA + X10_Harbor_seal_CR_2yrLead
  ',
  DAG1A_short_reduced = '
    X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES
    {SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + X12_copepodBiomass_WGoA + X10_Harbor_seal_CR_2yrLead
  ',
  DAG1B_long = '
    X07_DFA_cpue_IntSprJunHW ~ X09_DFA_HakeAge5Plus 
    {SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + {X06_VAR} + X12_DFA_biomassEuphShelfSum + X15_PacificCodBiomass_predAK
  ',
  DAG1B_long_reduced = '
    X07_DFA_cpue_IntSprJunHW ~ X09_DFA_HakeAge5Plus 
    {SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + {X06_VAR} + X12_DFA_biomassEuphShelfSum + X15_PacificCodBiomass_predAK
  ',
  DAG1B_short = '
    X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES 
    {SAR_VAR} ~ X07_DFA_cpue_IntSprJunHW + {X06_VAR} + X13_pollockBiomassGoAage3plus_predAK + X15_DFA_sleeperSharkBSAI_predAK
  ',
  DAG1B_short_reduced = '
    X07_DFA_cpue_IntSprJunHW ~ X08_commonMurre_JSOES 
    {SAR_VAR} ~ {X06_VAR} + X13_pollockBiomassGoAage3plus_predAK + X15_DFA_sleeperSharkBSAI_predAK
  '
)

# --- 2. Iteration Parameters ---
x06_options <- c("X06_DFA_IGF_mu", "X06_Lmu_IntSprJunH", "X06_StomFull_May")
target_responses <- c("X16_SAR", "sarSR", "sarUC")

# Create two data objects: full data, and data excluding 2001
# (Assumes guild.dfas1 has a column named 'date' or 'year' representing the timeline)
data_full <- guild.dfas1
data_no_2001 <- guild.dfas1 %>% filter(year != 2001)

datasets <- list(With_2001 = data_full, Without_2001 = data_no_2001)
raw_results <- list()

# --- 3. The Execution Loop ---
for (data_name in names(datasets)) {
  current_data <- datasets[[data_name]]
  
  for (dag in names(base_dags)) {
    for (x06_var in x06_options) {
      
      # If a model doesn't even contain an X06 slot (like DAG1A_long_reduced/short_reduced),
      # we evaluate it exactly once per response loop to keep our tables balanced.
      has_x06_slot <- grepl("\\{X06_VAR\\}", base_dags[[dag]])
      x06_loop_targets <- if (has_x06_slot) x06_var else x06_options[1]
      
      for (resp in target_responses) {
        
        # Build the functional syntax text dynamically
        model_syntax <- glue(base_dags[[dag]], X06_VAR = x06_var, SAR_VAR = resp)
        
        fit <- tryCatch({
          sem(model_syntax, data = current_data, std.lv = TRUE, missing = "ML")
        }, error = function(e) NULL)
        
        if (!is.null(fit)) {
          measures <- fitMeasures(fit)
          
          raw_results[[length(raw_results) + 1]] <- tibble(
            Timeline = data_name,
            DAG = dag,
            X06_Proxy = if (has_x06_slot) x06_var else "No X06 Slot",
            Response = resp,
            AIC = measures[["aic"]]
          )
        }
      }
      if (!has_x06_slot) break # Skip redundant variable variations if the model doesn't use X06
    }
  }
}

# Combine into a master results dataframe
combined_results <- bind_rows(raw_results)

# --- 4. Process into a Clean Pivot Layout comparing dAIC ---
final_pivot <- combined_results %>%
  # Group locally within each unique combination of Timeline and Response
  group_by(Timeline, Response) %>%
  # Calculate localized dAIC inside that specific group block
  mutate(dAIC = round(AIC - min(AIC, na.rm = TRUE), 1)) %>%
  ungroup() %>%
  # Format into rows for each DAG/Proxy configuration and columns for each Response
  select(Timeline, DAG, X06_Proxy, Response, dAIC) %>%
  pivot_wider(
    names_from = Response,
    values_from = dAIC,
    names_prefix = "dAIC_"
  ) %>%
  # Organize logically to evaluate the leverage of 2001 side-by-side
  arrange(Timeline, dAIC_X16_SAR)

# View the full nested pivot matrix
print(final_pivot, n = Inf)

write.csv(final_pivot, "LisaXP/outputs_4/X06comp_daic.csv",row.names = FALSE)


plot(1:nrow(final_pivot),final_pivot$dAIC_X16_SAR,type='l')
lines(1:nrow(final_pivot),final_pivot$dAIC_sarSR,col=2)
lines(1:nrow(final_pivot),final_pivot$dAIC_sarUC,col=3)


#heatmap of X06 aic comparison============
library(ggplot2)

# 1. Clean up and prep data for plotting (uses the 'combined_results' df from your loop)
plot_data_aic <- combined_results %>%
  group_by(Timeline, Response) %>%
  mutate(dAIC = AIC - min(AIC, na.rm = TRUE)) %>%
  ungroup() %>%
  # Simplify labels slightly for a cleaner plot window
  mutate(
    X06_Proxy = gsub("X06_", "", X06_Proxy),
    Response = factor(Response, levels = c("X16_SAR", "sarSR", "sarUC"))
  )

# 2. Build the Heatmap Matrix
ggplot(plot_data_aic, aes(x = X06_Proxy, y = DAG, fill = dAIC)) +
  geom_tile(color = "white", lwd = 0.4) +
  geom_text(aes(label = round(dAIC, 1)), color = "black", size = 3) +
  # Using a standard clean perceptual color scheme (darker = lower dAIC / better fit)
  scale_fill_viridis_c(direction = -1, option = "mako", name = "delta AIC") +
  # Split columns by data subsetting and rows by target time-series types
  facet_grid(Response ~ Timeline) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Sensitivity Analysis Matrix: X06 Selection Profiles",
    subtitle = "Comparing model fitting performance metrics (dAIC) across structural subsets",
    x = "X06 Proxy Alternatives Evaluated",
    y = "Structural DAG Architecture Base"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill = "gray95", color = "white"),
    strip.text = element_text(face = "bold")
  )


#Replace plot with re-ordered models and excluding the _reduced models===========
library(ggplot2)
library(dplyr)

# 1. Filter out reduced models and strictly order the remaining DAGs
plot_data_clean <- combined_results %>%
  # Eliminate all "_reduced" variants
  filter(!grepl("_reduced", DAG)) %>%
  # Calculate localized dAIC within each Timeline and Response block
  group_by(Timeline, Response) %>%
  mutate(dAIC = AIC - min(AIC, na.rm = TRUE)) %>%
  ungroup() %>%
  # Clean up labels for scannability
  mutate(
    X06_Proxy = gsub("X06_", "", X06_Proxy),
    Response = factor(Response, levels = c("X16_SAR", "sarSR", "sarUC")),
    # Strictly define the order of the remaining base models
    DAG = factor(DAG, levels = c(
      "DAG1A_short", "DAG1A_long",
      "DAG1B_short", "DAG1B_long",
      "DAG1C_short", "DAG1C_long"
    ))
  )

# 2. Generate the streamlined Heatmap Matrix
daic_plot_rmpreyNCC<-
ggplot(plot_data_clean, aes(x = X06_Proxy, y = DAG, fill = dAIC)) +
  geom_tile(color = "white", lwd = 0.5) +
  geom_text(aes(label = round(dAIC, 1)), color = "black", size = 3.5) +
  # Darker colors = lower dAIC (better fit)
  scale_fill_viridis_c(direction = -1, option = "mako", name = "delta AIC") +
  # Timeline on columns, Response types on rows
  facet_grid(Response ~ Timeline) +
  theme_minimal(base_size = 12) +
  labs(
    title = "X06 Model Selection Sensitivity Matrix",
    subtitle = "Comparing full structural architectures (Ordered A -> B -> C)",
    x = "X06 Proxy Alternatives Evaluated",
    y = "Structural Model Framework"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill = "gray93", color = "white"),
    strip.text = element_text(face = "bold", size = 11)
  )


ggsave("LisaXP/outputs_4/daic_plot_rmpreyNCC_rm2001.png",daic_plot_rmpreyNCC)
write.csv(plot_data_clean, "LisaXP/outputs_4/X06comp_daic_reduced.csv",row.names = FALSE)
