library(tidyverse)

# Here is a complete, ready-to-run diagnostic script that does exactly what you described:
#   
#   Calculates the pairwise correlations to define the Growth Proxies and Abundance Proxies.
# 
# Evaluates the top 100 models (p4my0) to see what percentage of them contain at least one proxy for both.
# 
# Builds a tidy, symmetric Co-occurrence (Pairing) Matrix to see which variables are working together.
# 
# Extracts and sorts the heaviest-hitting indicator pairs (those occurring together in multiple models).

# Results p4my0: both proxies are in 92% of all models
# HYPOTHESIS TEST: 92.0% of the Top 100 models contain both
# a Growth proxy and an Abundance proxy standing in for them.

      # proxy.ma<-as.data.frame(round(cor_matrix_full[,c("X06_Lmu_IntSprJunH","X07_DFA_cpue_IntSprJunHW") ],2))
      # proxy.ma %>% arrange(X06_Lmu_IntSprJunH) %>% filter(abs(X06_Lmu_IntSprJunH)>0.5 |abs(X07_DFA_cpue_IntSprJunHW)>0.5)
      #                                             X06_Lmu_IntSprJunH X07_DFA_cpue_IntSprJunHW
      # X13_ammod_WAI                                     -0.93                    -0.04
      # X10_Harbour_s_2yrLead_WS                          -0.79                    -0.16
      # X15_ArrowtoothFlounderBiomass_predAK              -0.71                     0.29
      # X05_herring_GAM                                    0.56                    -0.07
      # X14_pinkSalmonAsia                                 0.62                    -0.23
      # X06_Lmu_IntSprJunH                                 1.00                     0.01

#summary(hypothesis_testing[,c("Has_Growth_Proxy", "Has_Abundance_Proxy", "Meets_Both")])
      # Has_Growth_Proxy Has_Abundance_Proxy Meets_Both     
      # Mode :logical    Mode :logical       Mode :logical  
      # FALSE:1          FALSE:7             FALSE:8        
      # TRUE :99         TRUE :93            TRUE :92 

# 
# The Proxy & Co-occurrence Diagnostic Script

# --- 1. Isolate Your Target Group of 100 Models ---
p4my0 <- complex_models_only %>%
  filter(Num_Params == 4) %>%
  filter(Missing_Yrs == 0) %>%
  arrange(adjaic2) %>%
  head(100)

# Extract all unique predictors appearing in these 100 models
unique_predictors <- p4my0 %>%
  separate_rows(Predictors, sep = " \\+ ") %>%
  distinct(Predictors) %>%
  pull(Predictors)

# Ensure the core reference mediators are included in the pool for comparison
reference_mediators <- c("X06_Lmu_IntSprJunH", "X07_DFA_cpue_IntSprJunHW")
all_calc_vars <- unique(c(unique_predictors, reference_mediators))
all_calc_vars <- all_calc_vars[all_calc_vars %in% colnames(guild.dfas1)]

# --- 2. Calculate the Pairwise Correlation Matrix ---
cor_matrix_full <- cor(guild.dfas1[, all_calc_vars], use = "pairwise.complete.obs")

# Define stand-in correlation threshold (e.g., |r| >= 0.5)
# You can change this threshold to make the test more or less conservative
cor_threshold <- 0.50

growth_proxies <- names(which(abs(cor_matrix_full["X06_Lmu_IntSprJunH", ]) >= cor_threshold))
abundance_proxies <- names(which(abs(cor_matrix_full["X07_DFA_cpue_IntSprJunHW", ]) >= cor_threshold))

proxy.ma<-as.data.frame(round(cor_matrix_full[,c("X06_Lmu_IntSprJunH","X07_DFA_cpue_IntSprJunHW") ],2))
proxy.ma %>% arrange(X06_Lmu_IntSprJunH) %>% filter(abs(X06_Lmu_IntSprJunH)>0.5 |abs(X07_DFA_cpue_IntSprJunHW)>0.5)
#                                             X06_Lmu_IntSprJunH X07_DFA_cpue_IntSprJunHW
# X13_ammod_WAI                                     -0.93                    -0.04
# X10_Harbour_s_2yrLead_WS                          -0.79                    -0.16
# X15_ArrowtoothFlounderBiomass_predAK              -0.71                     0.29
# X05_herring_GAM                                    0.56                    -0.07
# X14_pinkSalmonAsia                                 0.62                    -0.23
# X06_Lmu_IntSprJunH                                 1.00                     0.01

# Print the proxy list to console for validation
cat("\n--- Growth Proxies (|r| >= ", cor_threshold, " with X06_Lmu) ---\n", sep = "")
print(growth_proxies)
# [1] "X13_ammod_WAI"                        "X06_Lmu_IntSprJunH"                   "X15_ArrowtoothFlounderBiomass_predAK" "X10_Harbour_s_2yrLead_WS"            
# [5] "X14_pinkSalmonAsia"                   "X05_herring_GAM"                     


cat("\n--- Abundance Proxies (|r| >= ", cor_threshold, " with X07_cpue) ---\n", sep = "")
print(abundance_proxies)
# [1] "X09_DFA_HakeAge5Plus"     "X07_DFA_cpue_IntSprJunHW"


# --- 3. Test the "Stand-in" Hypothesis ---
hypothesis_testing <- p4my0 %>%
  mutate(
    # Split the formula's predictors into a list of individual variables
    Pred_List = str_split(Predictors, " \\+ "),
    
    # Check if this model contains a Growth proxy (or Growth itself)
    Has_Growth_Proxy = map_lgl(Pred_List, ~ any(.x %in% growth_proxies)),
    
    # Check if this model contains an Abundance proxy (or Abundance itself)
    Has_Abundance_Proxy = map_lgl(Pred_List, ~ any(.x %in% abundance_proxies)),
    
    # Does it contain both?
    Meets_Both = Has_Growth_Proxy & Has_Abundance_Proxy
  )

# Calculate support percentage
percent_supported <- mean(hypothesis_testing$Meets_Both) * 100
summary(hypothesis_testing[,c("Has_Growth_Proxy", "Has_Abundance_Proxy", "Meets_Both")])
# Has_Growth_Proxy Has_Abundance_Proxy Meets_Both     
# Mode :logical    Mode :logical       Mode :logical  
# FALSE:1          FALSE:7             FALSE:8        
# TRUE :99         TRUE :93            TRUE :92 


cat("\n=========================================================\n")
cat(sprintf("HYPOTHESIS TEST: %.1f%% of the Top 100 models contain both\n", percent_supported))
cat("a Growth proxy and an Abundance proxy standing in for them.\n")
cat("=========================================================\n\n")

# HYPOTHESIS TEST: 92.0% of the Top 100 models contain both
# a Growth proxy and an Abundance proxy standing in for them.


# --- 4. Build the Pairwise Co-occurrence Matrix ---
# Create an empty square matrix for all active predictors in the top 100
active_predictors <- unique(unique_predictors)
n_active <- length(active_predictors)

co_occurrence <- matrix(0, 
                        nrow = n_active, 
                        ncol = n_active, 
                        dimnames = list(active_predictors, active_predictors))

# Populate the matrix by counting pairings across the 100 models
for (i in 1:nrow(p4my0)) {
  preds <- str_split(p4my0$Predictors[i], " \\+ ")[[1]]
  if (length(preds) > 1) {
    # Find all pairwise combinations of predictors in this model
    pairs <- combn(preds, 2, simplify = FALSE)
    for (pair in pairs) {
      p1 <- pair[1]
      p2 <- pair[2]
      co_occurrence[p1, p2] <- co_occurrence[p1, p2] + 1
      co_occurrence[p2, p1] <- co_occurrence[p2, p1] + 1
    }
  }
}


# --- 5. Extract the Heaviest-Hitting (Highly Paired) Couples ---
tidy_pairings <- as.data.frame(as.table(co_occurrence)) %>%
  rename(Var1 = Var1, Var2 = Var2, Models_Shared = Freq) %>%
  # Ensure we don't count self-pairing and remove mirror duplicates (Var1 < Var2)
  filter(as.character(Var1) < as.character(Var2)) %>%
  filter(Models_Shared > 0) %>%
  arrange(desc(Models_Shared))



#Order table
global_frequencies <- tidy_pairings %>%
  pivot_longer(cols = c(Var1, Var2), names_to = "Role", values_to = "Variable") %>%
  group_by(Variable) %>%
  summarise(Global_Frequency = sum(Models_Shared), .groups = "drop") %>%
  arrange(desc(Global_Frequency))

# Convert to a named vector for quick lookup (e.g., Variable -> Score)
freq_lookup <- deframe(global_frequencies)

# --- 2. Rearrange Var1 and Var2 so the Most Frequent is Always on the Left ---
ordered_pairings <- tidy_pairings %>%
  mutate(
    # Get global frequency scores for both sides
    Freq_V1 = freq_lookup[as.character(Var1)],
    Freq_V2 = freq_lookup[as.character(Var2)],
    
    # Keep track of original names to swap them if Var2 is more globally frequent
    New_Var1 = if_else(Freq_V1 >= Freq_V2, as.character(Var1), as.character(Var2)),
    New_Var2 = if_else(Freq_V1 >= Freq_V2, as.character(Var2), as.character(Var1)),
    
    # Store their numeric priority (for sorting)
    Sort_Priority_V1 = if_else(Freq_V1 >= Freq_V2, Freq_V1, Freq_V2),
    Sort_Priority_V2 = if_else(Freq_V1 >= Freq_V2, Freq_V2, Freq_V1)
  ) %>%
  # --- 3. Sort by Left-Side Frequency, then Right-Side Frequency, then Shared Models ---
  arrange(
    desc(Sort_Priority_V1), # Most frequent globally on the left gets top billing
    desc(Sort_Priority_V2), # Second most frequent globally gets second billing
    desc(Models_Shared)     # Tie-breaker
  ) %>%
  # Clean up and select final columns
  select(
    Var1 = New_Var1,
    Var2 = New_Var2,
    Models_Shared
  )

# --- 4. Print and Save ---
print(as_tibble(ordered_pairings), n=Inf)

ordered_pairings %>%
  filter(Models_Shared>5)
print(head(ordered_pairings, 20))
write_csv(ordered_pairings, "LisaXP/outputs_4/Ordered_Tight_Indicator_Pairings_Top100.csv")
