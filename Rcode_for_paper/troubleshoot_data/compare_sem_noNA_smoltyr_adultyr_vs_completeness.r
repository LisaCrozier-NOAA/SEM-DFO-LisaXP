library(tidyverse)

new_sem_data_file <- read.csv(file.path(proj_dir, "metadata/sem_data_noNA_1998_2021.csv"), stringsAsFactors = FALSE)
sem_orig     <- read.csv("data_Lisa/sem_master_data.csv", stringsAsFactors = FALSE)
pred_lookup <- read.csv(file.path(proj_dir, "metadata/all_pred_dfa_altprey.csv"), stringsAsFactors = FALSE)

# new_sem_data_file <- new_sem_data_file %>%
#   mutate(across(
#     .cols = -Year, 
#     .fns  = ~ .x - 5
#   ))

# Universal Min-Max Helper Function [0, 1]
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# Apply [0, 1] Min-Max Scaling across ALL numeric indicators EXCEPT targets & Year
new_sem_data_file <- new_sem_data_file %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), min_max_scale))

sem_orig <- sem_orig %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), min_max_scale))

# 1. Extract all unique target column names referenced in pred_lookup (excluding NA)
target_cols <- pred_lookup %>%
  select(pred_data_col, altprey1_data_col, altprey2_data_col, altprey3_data_col) %>%
  pivot_longer(cols = everything(), values_to = "target_var") %>%
  filter(!is.na(target_var) & target_var != "<NA>") %>%
  pull(target_var) %>%
  unique()

# Filter target_cols to those actually present in new_sem_data_file
target_cols <- intersect(target_cols, names(new_sem_data_file))

# 2. Build mapping dictionary to pair new_sem names with sem_orig names
mapping_df <- tibble(new_col = target_cols) %>%
  mutate(
    # Strip suffix (_smoltyr, _adultyr, or _pred)
    base_name = str_remove(new_col, "_smoltyr$|_adultyr$"),
    # Fix trailing underscores if any (e.g., matching X02_DFA_SeaNettle_)
    orig_candidate = case_when(
      base_name %in% names(sem_orig) ~ base_name,
      paste0(base_name, "_") %in% names(sem_orig) ~ paste0(base_name, "_"),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(orig_candidate)) # Keep only valid matches found in sem_orig

cat("Matched", nrow(mapping_df), "variable pairs out of", length(target_cols), "requested targets.\n")

# 3. Standardize Year column names and pivot to long format for ggplot
df_new_long <- new_sem_data_file %>%
  select(Year, any_of(mapping_df$new_col)) %>%
  pivot_longer(-Year, names_to = "new_col", values_to = "val_new")

df_orig_long <- sem_orig %>%
  # Handle capitalization of 'Year' vs 'year'
  rename_with(~"Year", matches("^year$", ignore.case = TRUE)) %>%
  select(Year, any_of(mapping_df$orig_candidate)) %>%
  pivot_longer(-Year, names_to = "orig_candidate", values_to = "val_orig")

# 4. Join dataset values together
comparison_df <- mapping_df %>%
  inner_join(df_new_long, by = "new_col") %>%
  inner_join(df_orig_long, by = c("Year", "orig_candidate")) %>%
  filter(!is.na(val_new) | !is.na(val_orig)) # Remove years where both are missing

# 5. Plot comparisons (Faceting in batches to keep plots readable)
# Create a plot label showing both column names being compared
comparison_df <- comparison_df %>%
  mutate(var_pair = paste0("New: ", new_col, "\nOrig: ", orig_candidate))

# Plot in batches of 9 pairs per page
unique_pairs <- unique(comparison_df$var_pair)
batch_size <- 9
num_batches <- ceiling(length(unique_pairs) / batch_size)

for (i in 1:num_batches) {
  pairs_sub <- unique_pairs[((i - 1) * batch_size + 1):min(i * batch_size, length(unique_pairs))]
  
  p <- comparison_df %>%
    filter(var_pair %in% pairs_sub) %>%
    ggplot(aes(x = Year)) +
    # Draw original line (dashed dark blue)
    geom_line(aes(y = val_orig, color = "sem_orig"), linewidth = 1, linetype = "dashed") +
    geom_point(aes(y = val_orig, color = "sem_orig"), size = 1.5, alpha = 0.7) +
    # Draw new line (solid red)
    geom_line(aes(y = val_new, color = "new_sem_data_file"), linewidth = 0.7, alpha = 0.8) +
    geom_point(aes(y = val_new, color = "new_sem_data_file"), size = 1.2, alpha = 0.8) +
    facet_wrap(~var_pair, scales = "free_y", ncol = 3) +
    scale_color_manual(
      values = c("sem_orig" = "darkblue", "new_sem_data_file" = "firebrick"),
      name = "Data Source"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 7),
      legend.position = "top"
    ) +
    labs(
      title = paste0("Time Series Validation (Batch ", i, " of ", num_batches, ")"),
      x = "Year",
      y = "Variable Value"
    )
  
  print(p)
}

# Find the missing target columns
missing_targets <- setdiff(target_cols, mapping_df$new_col)

cat("The following", length(missing_targets), "target columns had no match in sem_orig:\n")
print(missing_targets)

new_sem_data_file[,missing_targets]

# [1] "X05_eulachon_during_chinook_smoltyr"    "X05_eulachon_during_chinook_adultyr"    "X11_ssl_seak_pup_pred"                 
# [4] "X13_sitkaHerring_EGoA_smoltyr"          "X13_mid_il_capelin_smoltyr"             "X12_egoa_krill_smoltyr"                
# [7] "X13_pollock_age1plus_adultyr"           "X13_sitkaHerring_EGoA_adultyr"          "X13_pollock_age1plus_smoltyr"          
# [10] "X15_DFA_sleeperSharks_adultyr"          "X15_DFA_sleeperSharks_smoltyr"          "X10_AllSeaLionsEMB_2025_raw_sc_adultyr"