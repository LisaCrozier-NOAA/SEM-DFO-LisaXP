# Step 1: The Automated Layout Mapper Function
# This function takes the active model fit object and dynamically assigns the right names into your layout shapes based on your SEMnode classifications.
# Step 1 (Revised): Dynamic Node Extractor
# This function looks at the variables that survived in your specific model object (m_fit), matches them against your master lookup table, and outputs a clean, simple character vector where the names of the vector are your abstract layout roles (like "PreyNCC", "PredNCC", or "Growth"), and the values are the exact long variable names R needs to draw the boxes.


library(tidySEM)
library(ggplot2)
library(dplyr)

extract_nodes_for_plot_layout_fxn <- function(m_fit, lookup_df) {
#extract_active_layout_nodes <- function(m_fit, lookup_df) {
  # 1. Extract all observed variables that survived in this specific model
  active_vars <- lavaan::lavNames(m_fit, type = "ov")
  
  # 2. Filter your lookup table to ONLY include the surviving variables
  surviving_lookup <- lookup_df %>%
    filter(var %in% active_vars | Lisaname %in% active_vars)
  
  # 3. Build a dynamic named vector map
  # We match whatever long string lavaan is using to its abstract role (SEMnode)
  node_map <- surviving_lookup %>%
    mutate(actual_name = if_else(var %in% active_vars, var, Lisaname)) %>%
    select(SEMnode, actual_name) %>%
    distinct()
  
  # Convert to a named vector for effortless lookup downstream
  # Example output: c("X05" = "X05_DFA_abundSardine", "X16" = "X16_SAR")
  mapped_vector <- setNames(node_map$actual_name, node_map$SEMnode)
  
  return(mapped_vector)
}

# Test the function with your first model
test_map <- extract_active_layout_nodes(DAG1A_short_topmodel, var_lookup_NCC_AK)

# View the named vector results
print(test_map)
