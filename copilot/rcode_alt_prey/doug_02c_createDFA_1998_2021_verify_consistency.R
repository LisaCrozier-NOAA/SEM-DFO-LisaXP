# ==============================================================================
# Script: 03_verify_dfa_recreation.R (Matches X-prefixed DFA names & YMD dates)
# ==============================================================================

library(tidyverse)
library(lubridate)

rootdir <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP"

new_dfa_path   <- file.path(rootdir, "copilot/outputs_10/clusDataDFA.rds")
doug_comp_path <- file.path(rootdir, "2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/DAG1A_long/completeness.csv")

new_dfa   <- readRDS(new_dfa_path)
doug_comp <- read.csv(doug_comp_path)

# ------------------------------------------------------------------------------
# CLEAN & RESHAPE DOUG'S COMPLETENESS MATRIX
# ------------------------------------------------------------------------------
doug_comp_clean <- doug_comp %>%
  mutate(Year = year(ymd(date))) %>%
  select(-date) %>%
  pivot_longer(cols = -Year, names_to = "DFA_node", values_to = "orig_val") %>%
  filter(!is.na(orig_val)) %>%
  # Remove leading 'X' added by R's read.csv if present (e.g., X01.Zoo -> 01.Zoo)
  mutate(DFA_node = str_remove(DFA_node, "^X(?=\\d)"))

# ------------------------------------------------------------------------------
# CLEAN NEW DFA TRENDS
# ------------------------------------------------------------------------------
new_clean <- new_dfa %>%
  mutate(Year = year(date)) %>%
  select(Year, DFA_node = shortName, new_val = finalVal)

# ------------------------------------------------------------------------------
# JOIN & CALCULATE DIAGNOSTICS
# ------------------------------------------------------------------------------
joined <- inner_join(new_clean, doug_comp_clean, by = c("Year", "DFA_node"))

series_diagnostics <- joined %>%
  group_by(DFA_node) %>%
  summarise(
    N_Years       = n(),
    Max_Abs_Diff  = max(abs(new_val - orig_val), na.rm = TRUE),
    Mean_Abs_Diff = mean(abs(new_val - orig_val), na.rm = TRUE),
    Correlation   = cor(new_val, orig_val, use = "complete.obs"),
    Exact_Match   = all(abs(new_val - orig_val) < 1e-4),
    Sign_Inverted = all(abs(new_val + orig_val) < 1e-4),
    .groups       = "drop"
  )

cat("\n--- POST-DFA VERIFICATION SUMMARY ---\n")
cat("Total DFA Series Matched:", nrow(series_diagnostics), "\n")
cat("Exact Match Series (Diff < 0.0001):", sum(series_diagnostics$Exact_Match), "\n")
cat("Inverted Sign Series (Rotated Factor):", sum(series_diagnostics$Sign_Inverted), "\n")
cat("Mismatched Series:", sum(!series_diagnostics$Exact_Match & !series_diagnostics$Sign_Inverted), "\n\n")

if (any(!series_diagnostics$Exact_Match & !series_diagnostics$Sign_Inverted)) {
  cat("--- MISMATCHED DFA TRENDS ---\n")
  print(series_diagnostics %>% filter(!Exact_Match & !Sign_Inverted), n = 50)
} else {
  cat("SUCCESS: All post-DFA factor trends match Doug's completeness matrix exactly!\n")
}


#r mostly >0.995----------
#now plot-------------
# ==============================================================================
# Script: 04_plot_dfa_recreation_comparisons.R
# Purpose: Generate side-by-side time series plots comparing recreated DFAs
#          against Doug's original completeness baseline
# Output: Rcode_for_paper/Routput_for_paper/dfa_verification_plots.pdf
# ==============================================================================

library(tidyverse)
library(lubridate)
library(patchwork)

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & DIRECTORIES
# ------------------------------------------------------------------------------
rootdir <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP"

new_dfa_path   <- file.path(rootdir, "copilot/outputs_10/clusDataDFA.rds")
doug_comp_path <- file.path(rootdir, "2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/DAG1A_long/completeness.csv")

out_pdf_path   <- file.path(rootdir, "Rcode_for_paper/Routput_for_paper/dfa_verification_plots.pdf")
dir.create(dirname(out_pdf_path), showWarnings = FALSE, recursive = TRUE)

cat("Loading datasets for visual comparison...\n")

# ------------------------------------------------------------------------------
# STEP 1: LOAD & ALIGN DATASETS
# ------------------------------------------------------------------------------
new_dfa   <- readRDS(new_dfa_path)
doug_comp <- read.csv(doug_comp_path)

# Reshape Doug's completeness matrix into long format
doug_comp_long <- doug_comp %>%
  mutate(Year = year(ymd(date))) %>%
  select(-date) %>%
  pivot_longer(cols = -Year, names_to = "DFA_node", values_to = "orig_val") %>%
  filter(!is.na(orig_val)) %>%
  mutate(DFA_node = str_remove(DFA_node, "^X(?=\\d)"))

# Reshape New DFA output
new_clean <- new_dfa %>%
  mutate(Year = year(date)) %>%
  select(Year, DFA_node = shortName, new_val = finalVal)

# Join datasets
joined_df <- inner_join(new_clean, doug_comp_long, by = c("Year", "DFA_node")) %>%
  pivot_longer(
    cols = c(new_val, orig_val),
    names_to = "Source",
    values_to = "Value"
  ) %>%
  mutate(
    Source = if_else(Source == "new_val", "Recreated DFA", "Doug Baseline")
  )

# Get list of unique paired series
series_list <- sort(unique(joined_df$DFA_node))
cat(sprintf("Generating visual comparison plots for %d paired series...\n", length(series_list)))

# ------------------------------------------------------------------------------
# STEP 2: GENERATE MULTI-PANEL PDF REPORT (6 SERIES PER PAGE)
# ------------------------------------------------------------------------------
pdf(out_pdf_path, width = 11, height = 8.5)

plots_per_page <- 6
total_pages <- ceiling(length(series_list) / plots_per_page)

for (p in seq_len(total_pages)) {
  start_idx <- (p - 1) * plots_per_page + 1
  end_idx   <- min(p * plots_per_page, length(series_list))
  page_series <- series_list[start_idx:end_idx]
  
  plot_list <- list()
  
  for (s in page_series) {
    sub_df <- joined_df %>% filter(DFA_node == s)
    
    # Calculate correlation for plot annotation
    p_wide <- sub_df %>% 
      pivot_wider(names_from = Source, values_from = Value)
    r_val <- round(cor(p_wide$`Recreated DFA`, p_wide$`Doug Baseline`, use = "complete.obs"), 4)
    
    p_obj <- ggplot(sub_df, aes(x = Year, y = Value, color = Source, linetype = Source)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.8) +
      scale_color_manual(values = c("Doug Baseline" = "#1f78b4", "Recreated DFA" = "#e31a1c")) +
      scale_linetype_manual(values = c("Doug Baseline" = "solid", "Recreated DFA" = "dashed")) +
      scale_x_continuous(breaks = seq(1998, 2021, by = 4)) +
      labs(
        title = s,
        subtitle = paste0("r = ", r_val),
        x = "Year",
        y = "Standardized Value",
        color = "Dataset",
        linetype = "Dataset"
      ) +
      theme_minimal(base_size = 9) +
      theme(
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(color = "gray30", size = 8),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    plot_list[[s]] <- p_obj
  }
  
  # Combine plots onto page using patchwork
  combined_page <- wrap_plots(plot_list, ncol = 2, nrow = 3) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  print(combined_page)
}

dev.off()

cat("\nPlotting Complete! Visual comparison PDF exported to:\n", out_pdf_path, "\n")
