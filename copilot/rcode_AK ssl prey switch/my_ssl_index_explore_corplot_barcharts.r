library(tidyverse)
library(corrr)
library(ggcorrplot)
library(patchwork)

#SEE COOL INTERACTIVE PLOT AT THE END OF THIS SCRIPT



#raw data plot-------
library(tidyverse)

# 1. Identify all index columns (excluding metadata & raw sst)
all_index_cols <- setdiff(names(data_base), c("year", "i_ssl", "i_shark", "sst_wgoa_coastwatch_raw"))

# 2. Pivot into long format and scale unscaled shark indices dynamically
db_long <- data_base %>%
  select(year, all_of(all_index_cols)) %>%
  pivot_longer(
    cols = -year,
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    category = case_when(
      str_detect(variable, "(?i)ssl") ~ "1. Sea Lion (SSL) Indices",
      str_detect(variable, "(?i)(clim|sst|pdo|enso)") ~ "2. Temperature & Climate Indices",
      str_detect(variable, "(?i)(herr|herring)") ~ "3. Herring Indices",
      str_detect(variable, "(?i)(cap|capelin)") ~ "4. Capelin Indices",
      str_detect(variable, "(?i)shark") ~ "5. Shark Indices",
      TRUE ~ "Other"
    )
  ) %>%
  filter(category != "Other") %>%
  # Group by variable to standardise (z-score) any unscaled series
  group_by(variable) %>%
  mutate(
    # If a variable is unscaled (e.g., mean != 0 or standard deviation != 1), scale it
    value = if_else(
      category == "5. Shark Indices" & (abs(mean(value, na.rm = TRUE)) > 0.1 | sd(value, na.rm = TRUE) > 1.5),
      as.vector(scale(value)),
      value
    )
  ) %>%
  ungroup()

# 3. Faceted Plot with Legend in Empty 6th Slot
p_facet <- ggplot(db_long, aes(x = year, y = value, color = variable, group = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 1.2, alpha = 0.85) +
  facet_wrap(~ category, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(db_long$year, na.rm = TRUE), max(db_long$year, na.rm = TRUE), by = 2)) +
  theme_minimal(base_size = 10) +
  labs(
    title = "Input Index Dynamics by Ecological Category",
    subtitle = "All series standard-scaled (mean = 0, SD = 1) for direct trajectory comparison",
    x = "Year",
    y = "Scaled Value (Standard Deviations)",
    color = "Candidate Index"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.18), # Fits inside the empty 6th panel
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 6.5),
    legend.title = element_text(size = 8, face = "bold"),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.margin = margin(4, 4, 4, 4),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(ncol = 2, override.aes = list(linewidth = 1.2)))

print(p_facet)

ggsave("ssl_component_data_plots.png", p_facet)

# ==============================================================================
# 1. DEFINE FUNCTIONS-----
# ==============================================================================

ssl_index_fxn <- function(ssl_scaled, sst_scaled, herring_scaled, capelin_scaled,
                          scalar_ssl = 1, scalar_sst = -0.3,
                          scalar_herring = 2, scalar_capelin = 1) {
  I_SSL_raw = -1 * ( scalar_ssl * ssl_scaled 
                     + scalar_sst * (ssl_scaled * sst_scaled) 
                     + scalar_herring * (ssl_scaled * herring_scaled) 
                     + scalar_capelin * (ssl_scaled * capelin_scaled) )
  return(I_SSL_raw)
}

shark_index_fxn.v2 <- function(shark_scaled, 
                               temp_raw_Mt,
                               temp_ref_Mt = 0, 
                               temp_Ot,
                               Q10 = 2, 
                               overlap_form = "logistic", 
                               overlap_slope = 1) {
  shark_vec <- shark_scaled
  shark_z   <- shark_vec
  shark_z_roll2 = rowMeans(cbind(shark_z, dplyr::lag(shark_z, 1)), na.rm = TRUE)
  
  T_shark_raw <- temp_raw_Mt
  Tref_shark  <- temp_ref_Mt
  M_t         <- Q10^((T_shark_raw - Tref_shark) / 10)
  
  T_shark <- temp_Ot
  
  # Cleaned up nrow(data_base) to length(temp_Ot) for execution safety
  O_t <- if (overlap_form == "constant") {
    rep(1, length(temp_Ot))
  } else if (overlap_form == "linear") {
    pmax(0.1, 1 + overlap_slope * T_shark)
  } else if (overlap_form == "logistic") {
    plogis(overlap_slope * T_shark) * 2
  } else {
    return(NULL)
  }
  
  I_shark_z = shark_z * M_t * O_t
  return(I_shark_z)
}

# ==============================================================================
# 2. COLUMN CLASSIFICATION & CATEGORY CORRELATION MATRIX (GRAPH 1)--------
# ==============================================================================

# Classify variable names dynamically by category
all_cols <- setdiff(names(data_base), c("year", "i_ssl", "i_shark"))

ssl_cols    <- all_cols[str_detect(all_cols, "(?i)ssl")]
# Exclude "raw" SST from temperature options for SSL calculation
temp_cols   <- all_cols[str_detect(all_cols, "(?i)(clim|sst|pdo|enso)") & !str_detect(all_cols, "raw")]
herring_cols <- all_cols[str_detect(all_cols, "(?i)(herr|herring)")]
capelin_cols <- all_cols[str_detect(all_cols, "(?i)(cap|capelin)")]
shark_cols   <- all_cols[str_detect(all_cols, "(?i)shark")]

# Grouped correlation subset
cat_cols <- c(ssl_cols, temp_cols, herring_cols, capelin_cols, shark_cols)
corr_subset <- data_base %>% select(all_of(cat_cols))

# Compute Pearson correlations
cor_mat <- cor(corr_subset, use = "pairwise.complete.obs")

# Plot 1: Within & Between Category Correlation Heatmap
plot_cat_cor <- ggcorrplot(
  cor_mat, 
  hc.order = FALSE, 
  type = "full",
  lab = FALSE, 
  outline.col = "white",
  colors = c("#6D9EC1", "white", "#E46726"),
  title = "Correlation Structure Within & Between Index Categories"
) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8)
  )

print(plot_cat_cor)
ggsave("CorPlot_ssl_shark_sst_prey.png",plot_cat_cor)
write.csv(cor_mat,"CorPlot_ssl_shark_sst_prey.csv")
# ==============================================================================
# 3. PERMUTATION ENGINE (Generating all valid SSL index combinations)-----
# ==============================================================================

# Build grid of all variable combinations -- 1280 rows each specifying 1 permutation of the 4 variables in columns 
perm_grid <- expand.grid(
  ssl_var     = ssl_cols,
  sst_var     = temp_cols,
  herring_var = herring_cols,
  capelin_var = capelin_cols,
  stringsAsFactors = FALSE
)

# Function to compute predicted SSL index across time for a single permutation
calc_perm_series <- function(ssl_v, sst_v, herr_v, cap_v, df) {
  ssl_index_fxn(
    ssl_scaled     = df[[ssl_v]],
    sst_scaled     = df[[sst_v]],
    herring_scaled = df[[herr_v]],
    capelin_scaled = df[[cap_v]]
  )
}

# Compute predicted index matrix (Rows = Years, 1280 Columns = Permutation IDs)
perm_matrix <- map_dfc(1:nrow(perm_grid), function(i) {
  row <- perm_grid[i, ]
  pred <- calc_perm_series(row$ssl_var, row$sst_var, row$herring_var, row$capelin_var, data_base)
  tibble(!!paste0("perm_", i) := pred)
})

matplot(1998:2021,perm_matrix,type='l')
write.csv(perm_matrix,file="copilot/outputs_5/ssl_permutation_matrix_yrrow.csv")
write.csv(perm_grid,file="copilot/outputs_5/ssl_permutation_grid_4var_allcomb.csv")

# ==============================================================================
# 4. CORRELATIONS WITH BASELINE SSL AND SHARK INDEX (GRAPH 2)-------
# ==============================================================================

# Calculate correlation of each permutation vs. Baseline (i_ssl) & Shark (i_shark)
perm_results <- perm_grid %>%
  mutate(
    perm_id = paste0("perm_", row_number()),
    cor_with_i_ssl = map_dbl(perm_id, ~ cor(perm_matrix[[.x]], data_base$i_ssl, use = "pairwise.complete.obs")),
    cor_with_i_shark = map_dbl(perm_id, ~ cor(perm_matrix[[.x]], data_base$i_shark, use = "pairwise.complete.obs"))
  )

# Plot 2A: Distribution of Permutation Correlations
p_dist <- ggplot(perm_results) +
  geom_density(aes(x = cor_with_i_ssl, fill = "vs Baseline I_SSL"), alpha = 0.5) +
  geom_density(aes(x = cor_with_i_shark, fill = "vs Baseline I_Shark"), alpha = 0.5) +
  scale_fill_manual(values = c("vs Baseline I_SSL" = "#2b5c8f", "vs Baseline I_Shark" = "#d95f02")) +
  theme_minimal(base_size = 11) +
  labs(
    title = "A. Distribution of Correlations Across All SSL Index Permutations",
    x = "Pearson Correlation Coefficient (r)",
    y = "Density",
    fill = "Comparison Target"
  ) +
  theme(legend.position = "bottom")

# Plot 2B: SSL Index Permutation vs. Shark Index Coupling
p_scatter <- ggplot(perm_results, aes(x = cor_with_i_ssl, y = cor_with_i_shark, color = ssl_var)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black") +
  theme_minimal(base_size = 11) +
  labs(
    title = "B. Coupling: Correlation with I_SSL vs. Correlation with I_Shark",
    subtitle = "Color-coded by choice of SSL input variable",
    x = "Correlation with Baseline I_SSL",
    y = "Correlation with Baseline I_Shark",
    color = "SSL Variable Used"
  ) +
  theme(legend.position = "bottom")

# Combine Graph 2 panel using patchwork
plot_permutations <- p_dist / p_scatter + plot_layout(heights = c(1, 1.2))
print(plot_permutations)


#deconstruction---------
matplot(1998:2021,perm_matrix,type='l')

#How to See Which Variables Drive Positive vs. Negative Years
#To fix the issue where matplot gives you "spaghetti lines" without variable labels, here is an R script that:
#Deconstructs each year's prediction into variable-level feature contributions.
#Plots a Heatmap / Parallel Coordinate Plot showing exactly which variable choices (Capelin $A$ vs. Capelin $B$, SST $A$ vs. SST $B$) drive positive vs. negative $I_{\text{SSL}}$ outputs in key divergence years.R


library(tidyverse)

# 1. Tidy up the permutation matrix to track variable origins
perm_long <- perm_results %>%
  select(perm_id, ssl_var, sst_var, herring_var, capelin_var, cor_with_i_shark) %>%
  inner_join(
    perm_matrix %>% 
      mutate(year = data_base$year) %>% 
      pivot_longer(-year, names_to = "perm_id", values_to = "i_ssl_pred"),
    by = "perm_id"
  )

# Classify permutations by their relationship with I_Shark
perm_long <- perm_long %>%
  mutate(shark_coupling = if_else(cor_with_i_shark >= 0, "Positive Correlation with Shark", "Negative Correlation with Shark"))

# 2. Plot 1: Ribbon / Spaghetti Plot by Shark Coupling Group--------
# This clearly separates the trajectory of positively vs negatively correlated SSL indices over time
ggplot(perm_long, aes(x = year, y = i_ssl_pred, group = perm_id, color = capelin_var)) +
  geom_line(alpha = 0.4, linewidth = 0.8) +
  facet_wrap(~ shark_coupling, ncol = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Trajectory of I_SSL Predictions Over Time (1998–2021)",
    subtitle = "Separated by directional relationship with I_Shark & colored by Capelin Index choice",
    x = "Year",
    y = "Predicted Raw I_SSL Index",
    color = "Capelin Variable Used"
  ) +
  theme(legend.position = "bottom")


# 3. Plot 2: Variable Driver Heatmap in Divergent Years (e.g., years with highest variance)------
# Identify the 5 years where the permutations differ the most
high_var_years <- perm_long %>%
  group_by(year) %>%
  summarise(var_pred = var(i_ssl_pred, na.rm = TRUE)) %>%
  slice_max(var_pred, n = 6) %>%
  pull(year)

# Plot a heatmap of permutations vs year for high-variance years
ggplot(perm_long %>% filter(year %in% high_var_years), 
       aes(x = factor(year), y = reorder(perm_id, cor_with_i_shark), fill = i_ssl_pred)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2b5c8f", mid = "white", high = "#d95f02", midpoint = 0) +
  facet_grid(capelin_var ~ sst_var, scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 10) +
  labs(
    title = "Drivers of I_SSL Directionality in High-Divergence Years",
    subtitle = "Faceted by Capelin (Rows) and SST (Columns) variable choices",
    x = "High Divergence Year",
    y = "Permutation (Ordered by Correlation with I_Shark)",
    fill = "I_SSL Value"
  ) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank())
#better ribbon plot-------
library(tidyverse)

# 1. Attach Raw SST to the long permutation predictions
perm_long_sst <- perm_results %>%
  select(perm_id, ssl_var, sst_var, herring_var, capelin_var) %>%
  inner_join(
    perm_matrix %>% 
      mutate(
        year = data_base$year,
        sst_raw = data_base$sst_wgoa_coastwatch_raw
      ) %>% 
      pivot_longer(
        cols = starts_with("perm_"), 
        names_to = "perm_id", 
        values_to = "i_ssl_pred"
      ),
    by = "perm_id"
  )

# 2. Categorize capelin variable and EXPLICITLY ORDER factor levels
perm_long_sst <- perm_long_sst %>%
  mutate(
    capelin_group = if_else(
      capelin_var == "wgoa_bio_pav_capelin_cpue",
      "Pavlof Capelin Index (wgoa_bio_pav_capelin_cpue)",
      "All Other Capelin Indices"
    ),
    # Enforce exact panel order: Pavlof first (Top), Others second (Bottom)
    capelin_group = factor(
      capelin_group, 
      levels = c(
        "Pavlof Capelin Index (wgoa_bio_pav_capelin_cpue)",
        "All Other Capelin Indices"
      )
    )
  )

# 3. Plot: Predicted I_SSL vs Raw SST
ggplot(
  perm_long_sst, 
  aes(x = sst_raw, y = i_ssl_pred, group = perm_id, color = ssl_var)
) +
  geom_point(alpha = 0.3, size = 1.2) +
  geom_line(alpha = 0.4, linewidth = 0.7) +
  facet_wrap(~ capelin_group, ncol = 1, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 12) +
  labs(
    title = "I_SSL Predictions vs. Raw SST Controlled by Capelin Index Choice",
    subtitle = "Top Panel: Pavlof CPUE | Bottom Panel: Other Capelin Indices | Colors: SSL Count Index",
    x = "Raw Water Temperature (sst_wgoa_coastwatch_raw)",
    y = "Predicted Raw I_SSL Index",
    color = "SSL Variable Used"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 11)
  )
# Note: If you prefer coloring lines by the temperature index instead of the ssl index, 
# simply replace `color = ssl_var` with `color = sst_var` in the aes() mapping above.

# 4. Plot: Predicted I_SSL vs Raw SST w/ bands= SST index
ggplot(
  perm_long_sst, 
  aes(x = sst_raw, y = i_ssl_pred, group = perm_id, color = sst_var)
) +
  geom_point(alpha = 0.3, size = 1.2) +
  geom_line(alpha = 0.4, linewidth = 0.7) +
  facet_wrap(~ capelin_group, ncol = 1, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 12) +
  labs(
    title = "I_SSL Predictions vs. Raw SST Controlled by Capelin Index Choice",
    subtitle = "Top vs. Bottom panels isolate Pavlof CPUE; Color highlights SST Index choice",
    x = "Raw Water Temperature (sst_wgoa_coastwatch_raw)",
    y = "Predicted Raw I_SSL Index",
    color = "SST Variable Used"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 11)
  )


#Facet grid of SSL index vs capelin----------
library(tidyverse)

# Facet grid: SSL variables as rows, Capelin variables as columns
ggplot(perm_long_sst, aes(x = sst_raw, y = i_ssl_pred, group = perm_id, color = sst_var)) +
  geom_line(alpha = 0.5, linewidth = 0.6) +
  facet_grid(ssl_var ~ capelin_var) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal(base_size = 10) +
  labs(
    title = "I_SSL Predictions vs. Raw SST across SSL and Capelin Variable Combinations",
    subtitle = "Rows = SSL Index Input | Columns = Capelin Index Input | Colors = SST Variable Input",
    x = "Raw Water Temperature (sst_wgoa_coastwatch_raw)",
    y = "Predicted Raw I_SSL Index",
    color = "SST Variable"
  ) +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0, size = 7), # Horizontal labels for readability
    strip.text.x = element_text(size = 7),
    strip.background = element_rect(fill = "grey92", color = NA)
  )


#table of SSL index predictions for 5 most variables years-----
library(tidyverse)

# 1. Identify the 5 most variable years based on predicted I_SSL variance
top5_years <- perm_long_sst %>%
  group_by(year) %>%
  summarise(var_pred = var(i_ssl_pred, na.rm = TRUE)) %>%
  slice_max(var_pred, n = 5) %>%
  pull(year)

cat("The 5 most variable years are:", paste(top5_years, collapse = ", "), "\n\n")

# 2. Compute Mean and SD for each SSL index across all other variable permutations
ssl_summary_table <- perm_long_sst %>%
  filter(year %in% top5_years) %>%
  group_by(year, ssl_var) %>%
  summarise(
    mean_prediction = mean(i_ssl_pred, na.rm = TRUE),
    sd_prediction   = sd(i_ssl_pred, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Round for clean display
  mutate(
    mean_prediction = round(mean_prediction, 3),
    sd_prediction   = round(sd_prediction, 3)
  )

# 3. Format as a clean wide table for easy comparison across years
ssl_summary_wide <- ssl_summary_table %>%
  pivot_wider(
    names_from = year,
    values_from = c(mean_prediction, sd_prediction),
    names_glue = "Yr_{year}_{.value}"
  )

# Print formatted table to console
print(ssl_summary_wide %>% select(ssl_var,contains("mean")), n = Inf)
print(ssl_summary_wide %>% select(ssl_var,contains("sd")), n = Inf)

#barchart of ssl predictions-------
library(tidyverse)

# 1. Identify the 5 most variable years
top5_years <- perm_long_sst %>%
  group_by(year) %>%
  summarise(var_pred = var(i_ssl_pred, na.rm = TRUE)) %>%
  slice_max(var_pred, n = 5) %>%
  pull(year)

# 2. Calculate the mean prediction for each SSL index within each of those 5 years
ssl_means <- perm_long_sst %>%
  filter(year %in% top5_years) %>%
  group_by(year, ssl_var) %>%
  summarise(
    mean_prediction = mean(i_ssl_pred, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Plot: Grouped Bar Chart with Line Overlays
ggplot(ssl_means, aes(x = factor(year), y = mean_prediction, fill = ssl_var)) +
  # Grouped bars to show mean values side-by-side per year
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
  # Connecting lines for each SSL index across the 5 years
  geom_line(
    aes(group = ssl_var, color = ssl_var), 
    position = position_dodge(width = 0.8), 
    linewidth = 1
  ) +
  geom_point(
    aes(group = ssl_var, color = ssl_var), 
    position = position_dodge(width = 0.8), 
    size = 2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_brewer(palette = "Set2", name = "SSL Variable") +
  scale_color_brewer(palette = "Set2", name = "SSL Variable") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Mean Predicted I_SSL across the 5 Most Variable Years",
    subtitle = "Bars show mean prediction per SSL index; Lines connect trajectory across years",
    x = "Year (Ordered by High Variance)",
    y = "Mean Predicted I_SSL Index"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

#rank order years by global mean prediction-------
library(tidyverse)

# 1. Identify and rank the top 5 years by GLOBAL mean prediction (across ALL permutations)
top5_ranked_years <- perm_long_sst %>%
  group_by(year) %>%
  summarise(
    global_mean = mean(i_ssl_pred, na.rm = TRUE),
    global_sd   = sd(i_ssl_pred, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  slice_max(global_mean, n = 5) %>%
  arrange(desc(global_mean)) %>%
  pull(year)

cat("Top 5 years ranked by highest global mean prediction:", paste(top5_ranked_years, collapse = ", "), "\n\n")

# 2. Calculate the mean prediction for each SSL index within those ranked years
ssl_means_ranked <- perm_long_sst %>%
  filter(year %in% top5_ranked_years) %>%
  group_by(year, ssl_var) %>%
  summarise(
    mean_prediction = mean(i_ssl_pred, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Set year as an ordered factor matching the global mean rank
  mutate(year_factor = factor(year, levels = top5_ranked_years))

# 3. Plot: Ranked Bar Chart with Line Overlays
ggplot(ssl_means_ranked, aes(x = year_factor, y = mean_prediction, fill = ssl_var)) +
  # Grouped bars to show individual SSL index means side-by-side per year
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
  # Connecting lines for each SSL index across the ranked years
  geom_line(
    aes(group = ssl_var, color = ssl_var), 
    position = position_dodge(width = 0.8), 
    linewidth = 1
  ) +
  geom_point(
    aes(group = ssl_var, color = ssl_var), 
    position = position_dodge(width = 0.8), 
    size = 2.5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_brewer(palette = "Set2", name = "SSL Variable") +
  scale_color_brewer(palette = "Set2", name = "SSL Variable") +
  theme_minimal(base_size = 12) +
  labs(
    title = "I_SSL Predictions by SSL Index Across Top 5 Ranked Years",
    subtitle = "Years ranked from left to right by highest overall global mean I_SSL prediction",
    x = "Year (Ranked High to Low by Global Mean)",
    y = "Mean Predicted I_SSL Index"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  )
#try to, barchart of most variable years-----
library(tidyverse)

# 1. Identify the 5 years with the HIGHEST VARIANCE
high_var_years <- perm_long_sst %>%
  group_by(year) %>%
  summarise(var_pred = var(i_ssl_pred, na.rm = TRUE)) %>%
  slice_max(var_pred, n = 5) %>%
  pull(year)

# 2. Rank those 5 high-variance years by their GLOBAL MEAN prediction
ranked_high_var_years <- perm_long_sst %>%
  filter(year %in% high_var_years) %>%
  group_by(year) %>%
  summarise(global_mean = mean(i_ssl_pred, na.rm = TRUE)) %>%
  arrange(desc(global_mean)) %>%
  pull(year)

cat("Top 5 highest variance years, ordered by global mean:", paste(ranked_high_var_years, collapse = ", "), "\n\n")

# 3. Calculate mean predictions for each SSL index in these ordered years
ssl_means_ordered <- perm_long_sst %>%
  filter(year %in% ranked_high_var_years) %>%
  group_by(year, ssl_var) %>%
  summarise(
    mean_prediction = mean(i_ssl_pred, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Convert year to an ordered factor matching the rank order
  mutate(year_factor = factor(year, levels = ranked_high_var_years))

# 4. Plot: Bar Chart with Line Overlays
ggplot(ssl_means_ordered, aes(x = year_factor, y = mean_prediction, fill = ssl_var)) +
  # Grouped bars to show individual SSL index means side-by-side
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
  # Connecting lines for each SSL index across the ordered years
  geom_line(
    aes(group = ssl_var, color = ssl_var), 
    position = position_dodge(width = 0.8), 
    linewidth = 1
  ) +
  geom_point(
    aes(group = ssl_var, color = ssl_var), 
    position = position_dodge(width = 0.8), 
    size = 2.5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_brewer(palette = "Set2", name = "SSL Variable") +
  scale_color_brewer(palette = "Set2", name = "SSL Variable") +
  theme_minimal(base_size = 12) +
  labs(
    title = "I_SSL Predictions Across the 5 Most Variable Years",
    subtitle = "Subset: 5 highest-variance years | X-Axis Order: Ranked by highest to lowest global mean",
    x = "Year (Highest Variance, Ordered by Global Mean)",
    y = "Mean Predicted I_SSL Index"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

#pick the anomalous index----------
# Find which SSL index flipped positive in 2020 and 2021
perm_long_sst %>%
  filter(year %in% c(2020, 2021)) %>%
  group_by(year, ssl_var) %>%
  summarise(
    mean_prediction = mean(i_ssl_pred, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = year, values_from = mean_prediction, names_prefix = "Yr_") %>%
  mutate(is_outlier = if_else(Yr_2020 > 0 | Yr_2021 > 0, "<- OUTLIER INDEX", ""))


# ssl_var           Yr_2020 Yr_2021 is_outlier        
#1 ssl_cent_np_pred   -7.64   -6.88  ""                
# 2 ssl_cent_pup_pred  -8.51   -8.74  ""                
# 3 ssl_count_eric    NaN     NaN      NA               
# 4 ssl_east_np_pred   -0.383   0.511 "<- OUTLIER INDEX"
# 5 ssl_east_pup_pred  -5.36   -6.10  ""                
# 6 ssl_model_eric     -9.99  -10.7   ""                
# 7 ssl_seak_np_pred   -4.05   -3.67  ""                
# 8 ssl_seak_pup_pred  -3.93   -3.66  ""                
# 9 ssl_west_np_pred   -6.16   -6.29  ""                
# 10 ssl_west_pup_pred  -7.15   -7.13  ""                

#all years------------

library(tidyverse)

# 1. Rank ALL years by their GLOBAL mean prediction across all variable permutations
ranked_all_years <- perm_long_sst %>%
  group_by(year) %>%
  summarise(global_mean = mean(i_ssl_pred, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(global_mean)) %>%
  pull(year)

# 2. Calculate mean prediction for each SSL index within every year
ssl_means_all <- perm_long_sst %>%
  group_by(year, ssl_var) %>%
  summarise(
    mean_prediction = mean(i_ssl_pred, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Convert year to an ordered factor using the global mean rank order
  mutate(year_factor = factor(year, levels = ranked_all_years))

# 3. Plot: Full Time Series Ranked by Global Mean Output
i_ssl_all_years_plot<-
ggplot(ssl_means_all, aes(x = year_factor, y = mean_prediction, fill = ssl_var)) +
  # Grouped bars to show individual SSL index means side-by-side per year
  geom_col(position = position_dodge(width = 0.85), width = 0.75, alpha = 0.85) +
  # Connecting lines for each SSL index across all ranked years
  geom_line(
    aes(group = ssl_var, color = ssl_var), 
    position = position_dodge(width = 0.85), 
    linewidth = 0.8
  ) +
  geom_point(
    aes(group = ssl_var, color = ssl_var), 
    position = position_dodge(width = 0.85), 
    size = 1.8
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_brewer(palette = "Set2", name = "SSL Variable") +
  scale_color_brewer(palette = "Set2", name = "SSL Variable") +
  theme_minimal(base_size = 11) +
  labs(
    title = "I_SSL Predictions Across All Years Ranked by Global Mean",
    subtitle = "X-Axis: Ranked from highest overall average I_SSL prediction to lowest | Bars: SSL Index Variant Mean",
    x = "Year (Ranked High to Low by Global Mean I_SSL)",
    y = "Mean Predicted I_SSL Index"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave("ssl_index_by_mean_all_years_plot.png",i_ssl_all_years_plot)

#normal year plot of i_ssl---------
library(tidyverse)

# 1. Calculate the mean prediction for each SSL index in chronological year order
ssl_dots_data <- perm_long_sst %>%
  group_by(year, ssl_var) %>%
  summarise(
    mean_prediction = mean(i_ssl_pred, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Calculate yearly min/max to draw subtle vertical connectors for each year
yearly_ranges <- ssl_dots_data %>%
  group_by(year) %>%
  summarise(
    min_pred = min(mean_prediction),
    max_pred = max(mean_prediction),
    .groups = "drop"
  )

# 3. Plot: Chronological Dot Timeline
i_ssl_all_years_chron_order_plot <-
ggplot(ssl_dots_data, aes(x = year, y = mean_prediction)) +
  # Light grey vertical segments connecting the dots within each year to show the spread
  geom_segment(
    data = yearly_ranges,
    aes(x = year, xend = year, y = min_pred, yend = max_pred),
    color = "grey75", linewidth = 0.8, linetype = "solid"
  ) +
  # Zero baseline
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  # Colored dots for each SSL index variant
  geom_point(aes(color = ssl_var), size = 3, alpha = 0.9) +
  scale_color_brewer(palette = "Set2", name = "SSL Variable") +
  scale_x_continuous(breaks = seq(min(ssl_dots_data$year), max(ssl_dots_data$year), by = 1)) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Chronological Timeline of Predicted I_SSL Index Means",
    subtitle = "Points show mean prediction per SSL Index Variant; Grey lines highlight annual spread",
    x = "Year",
    y = "Mean Predicted I_SSL Index",
    color = "SSL Variable Used"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(i_ssl_all_years_chron_order_plot)

ggsave("ssl_index_all_years_chron_order_plot.png",i_ssl_all_years_chron_order_plot)


#herring explains the amplification of ssl index divergence in 2020 and 2021-------
library(tidyverse)

# Check relationship between Herring levels and Index Disagreement Range
herring_vs_range <- perm_long_sst %>%
  group_by(year) %>%
  summarise(
    mean_herring_sd = mean(data_base$egoa_bio_stka_herr_matbiom[data_base$year == first(year)], na.rm = TRUE),
    index_range     = max(i_ssl_pred, na.rm = TRUE) - min(i_ssl_pred, na.rm = TRUE),
    .groups         = "drop"
  )

ggplot(herring_vs_range, aes(x = mean_herring_sd, y = index_range)) +
  geom_point(size = 3, color = "#2b5c8f") +
  geom_text(aes(label = year), vjust = -0.7, size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40", linetype = "dashed") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Herring Abundance vs. Annual Prediction Disagreement Range",
    x = "Herring Scaled Value (SD)",
    y = "Prediction Range (Max I_SSL - Min I_SSL)"
  )



#show temp index effects on all years-----------
library(tidyverse)

# 1. Filter for 1998–2018 and lock down to `ssl_seak_pup_pred`
dot_band_data <- perm_long_sst %>%
  filter(
    year >= 1998 & year <= 2018,
    ssl_var == "ssl_seak_pup_pred"
  )

# 2. Calculate yearly stats for background context
ggplot(dot_band_data, aes(x = factor(year), y = i_ssl_pred)) +
  # Light boxplot underneath to show the interquartile range for each year
  geom_boxplot(
    outlier.shape = NA, 
    alpha = 0.25, 
    fill = "grey80", 
    color = "grey50", 
    width = 0.5
  ) +
  # Zero baseline
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  # Staggered/jittered points colored by Temperature Index
  # position_jitterdodge or position_jitter keeps points from overlapping vertically
  geom_point(
    aes(color = sst_var), 
    position = position_jitter(width = 0.22, height = 0, seed = 123),
    size = 2.2, 
    alpha = 0.85
  ) +
  scale_color_brewer(palette = "Dark2", name = "Temperature Index") +
  theme_minimal(base_size = 11) +
  labs(
    title = "I_SSL Predictions Across Temperature Indices (1998–2018)",
    subtitle = "Fixed SSL Variable: ssl_seak_pup_pred | Points staggered to show index bands",
    x = "Year",
    y = "Predicted Raw I_SSL Index",
    color = "Temperature Index Used"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) # Wraps legend entries neatly


#try again-------
library(tidyverse)

# 1. Filter dataset: 
#    - Exclude Pavlof capelin index
#    - Years 1998-2018
#    - Fixed SSL variable: ssl_seak_pup_pred
clean_dot_data <- perm_long_sst %>%
  filter(
    capelin_var != "wgoa_bio_pav_capelin_cpue",
    year >= 1998 & year <= 2018,
    ssl_var == "ssl_seak_pup_pred"
  )

# 2. Plot: Ordered, dodge-aligned bands by Temperature Index
ggplot(clean_dot_data, aes(x = factor(year), y = i_ssl_pred)) +
  # Zero baseline reference
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  # Light boxplot in background for overall spread per year
  geom_boxplot(
    outlier.shape = NA, 
    alpha = 0.2, 
    fill = "grey80", 
    color = "grey60", 
    width = 0.6
  ) +
  # Line up points neatly into organized bands per temperature index
  geom_point(
    aes(color = sst_var, group = sst_var), 
    position = position_dodge(width = 0.7),
    size = 2.2, 
    alpha = 0.9
  ) +
  scale_color_brewer(palette = "Dark2", name = "Temperature Index") +
  theme_minimal(base_size = 11) +
  labs(
    title = "I_SSL Predictions Across Temperature Indices (1998–2018)",
    subtitle = "Excludes wgoa_bio_pav_capelin_cpue | Fixed SSL: ssl_seak_pup_pred | Temperature Indices Dodged",
    x = "Year",
    y = "Predicted Raw I_SSL Index",
    color = "Temperature Index Used"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


#still clumping eg 2014-----
ggplot(clean_dot_data, aes(x = factor(year), y = i_ssl_pred)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_point(
    aes(color = sst_var, shape = herring_var, group = interaction(sst_var, herring_var)), 
    position = position_dodge(width = 0.8),
    size = 2.5, 
    alpha = 0.9
  ) +
  scale_color_brewer(palette = "Dark2", name = "Temperature Index") +
  scale_shape_discrete(name = "Herring Index") +
  theme_minimal(base_size = 11) +
  labs(
    title = "I_SSL Predictions (1998–2018): Temperature (Color) & Herring (Shape)",
    subtitle = "Excludes Pavlof Capelin | Fixed SSL: ssl_seak_pup_pred",
    x = "Year",
    y = "Predicted Raw I_SSL Index"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.minor = element_blank()
  )


#INTERACTIVE PLOT! html hover over any point for exact permutation id-----------
library(tidyverse)
library(plotly)

# 1. Filter dataset: Exclude Pavlof capelin, lock SSL to ssl_seak_pup_pred, years 1998-2018
clean_dot_data <- perm_long_sst %>%
  filter(
    capelin_var != "wgoa_bio_pav_capelin_cpue",
    capelin_var != "capelin_avg",
    year >= 1998 & year <= 2018,
    ssl_var == "ssl_seak_pup_pred"
  )

# 2. Build interactive plot
p <- ggplot(
  clean_dot_data, 
  aes(
    x = factor(year), 
    y = i_ssl_pred, 
    color = sst_var,
    # Custom hover text to inspect any point
    text = paste0(
      "Year: ", year,
      "<br>Permutation ID: ", perm_id,
      "<br>SST Var: ", sst_var,
      "<br>Herring Var: ", herring_var,
      "<br>Capelin Var: ", capelin_var,
      "<br>Predicted Index: ", round(i_ssl_pred, 3)
    )
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(
    aes(group = interaction(sst_var, herring_var, capelin_var)),
    position = position_dodge(width = 0.6),
    size = 2.5,
    alpha = 0.85
  ) +
  scale_color_brewer(palette = "Dark2", name = "Temperature Index") +
  theme_minimal(base_size = 11) +
  labs(
    title = "Interactive I_SSL Predictions (1998–2018)",
    subtitle = "Hover over points to inspect exact variable combinations",
    x = "Year",
    y = "Predicted Raw I_SSL Index"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "bottom"
  )

# 3. Convert to interactive HTML object
interactive_plot <- ggplotly(p, tooltip = "text")
interactive_plot

#what was pos vs neg in 2014?---------
#wgoa_bio_mid_il_capelin was pos, wgoa_cap_pcod_trend  was negative
y2014 <- perm_long_sst %>%
  filter(
    capelin_var != "wgoa_bio_pav_capelin_cpue",
    capelin_var != "capelin_avg",
    year == 2014,
    ssl_var == "ssl_seak_pup_pred"
  )

y2014
names(y2014)

neg<-y2014 %>% filter(i_ssl_pred<0)
pos<-y2014 %>% filter(i_ssl_pred>0)

rbind(table(neg$herring_var),table(pos$herring_var))
rbind(table(neg$capelin_var),table(pos$capelin_var))
rbind(table(neg$sst_var),table(pos$sst_var))
rbind(table(neg$ssl_var),table(pos$ssl_var))


#ok, I need to move on!!!!!---------