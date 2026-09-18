library(dplyr)
library(tidyr)
library(ggplot2)

# Define directories
output_dir <- file.path("copilot/outputs_9/diagnostic_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load datasets
sem_data         <- read.csv("data_Lisa/sem_master_data.csv", row.names = NULL)
sem_legacy <- read.csv("data_Lisa/sem_master_data.csv", row.names = NULL)
sem_new    <- read.csv("copilot/outputs_9/sem_altprey_data.csv", row.names = NULL)

# Standardize year column names
if ("year" %in% names(sem_legacy)) sem_legacy <- rename(sem_legacy, Year = year)
if ("year" %in% names(sem_new))    sem_new    <- rename(sem_new, Year = year)

# Filter 1998-2021 analysis window
sem_legacy_win <- sem_legacy %>% filter(Year >= 1998 & Year <= 2021)
sem_new_win    <- sem_new    %>% filter(Year >= 1998 & Year <= 2021)

# Identify common columns (excluding Year)
common_cols <- intersect(names(sem_legacy_win), names(sem_new_win))
common_cols <- setdiff(common_cols, "Year")

cat(sprintf("Found %d matching time series between sem_data and sem_altprey_data.\n", length(common_cols)))

# ---------------------------------------------------------------------------
# Loop through each matching variable, calculate correlation, and plot
# ---------------------------------------------------------------------------
summary_stats <- data.frame(Variable = character(), Correlation = numeric(), stringsAsFactors = FALSE)

var=common_cols[33]
n<-grep("_2yrLead",common_cols);n #33 34 35 36 57
n<-grep("X15",common_cols);n #33 34 35 36 57

#discrepancies:
#mammals 22:26
#jsoes plankton
#X05_DFA_abudnSardine!!! just a little different, not lagged
#X03_DFA_comMurreDietHerrSard -- TOTALLY DIFFERENT!!!!!
#THAT ONE CANNOT BE USED R=-0.5

#x11_dfa_hARBOUR_P_ws TOtally lack data in the updated file. DONT USE

length(common_cols)

n=1:5
n=6:10
n=11:20
n=21:32
n=37:56
n=57:62

for (var in common_cols[n]) {
  # Build comparison frame
  df_comp <- data.frame(
    Year     = sem_legacy_win$Year,
    Legacy   = sem_legacy_win[[var]],
    Updated  = sem_new_win[[var]]
  ) %>%
    filter(!is.na(Legacy) & !is.na(Updated))
  
  # Calculate Pearson correlation
  corr_val <- if(nrow(df_comp) > 2) cor(df_comp$Legacy, df_comp$Updated) else NA
  summary_stats <- rbind(summary_stats, data.frame(Variable = var, Correlation = round(corr_val, 4)))
  
  # Reshape for ggplot
  df_long <- df_comp %>%
    pivot_longer(cols = c("Legacy", "Updated"), names_to = "Dataset", values_to = "Value")
  
  # Plot overlaid time series
  p <- ggplot(df_long, aes(x = Year, y = Value, color = Dataset, linetype = Dataset)) +
    geom_line(size = 1.1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Legacy" = "firebrick", "Updated" = "dodgerblue3")) +
    scale_linetype_manual(values = c("Legacy" = "solid", "Updated" = "dashed")) +
    labs(
      title = paste("Comparison:", var),
      subtitle = sprintf("Pearson Correlation (1998-2021): r = %.4f", corr_val),
      x = "Year",
      y = "Normalized Indicator Value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  print(p)
  
  # # Save individual plot image
  # file_name <- paste0("compare_", gsub("[^A-Za-z0-9_]", "_", var), ".png")
  # ggsave(file.path(output_dir, file_name), plot = p, width = 8, height = 4.5, dpi = 150)
}

# ---------------------------------------------------------------------------
# Print Correlation Summary
# ---------------------------------------------------------------------------
cat("\n=== TIME SERIES ALIGNMENT SUMMARY (1998-2021) ===\n")
print(summary_stats, row.names = FALSE)

# Export summary table
write.csv(summary_stats, file.path("copilot/outputs_9", "time_series_alignment_correlations.csv"), row.names = FALSE)