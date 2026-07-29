library(dplyr)
library(tidyr)
library(ggplot2)
library(lavaan)

generate_model_time_series_plot_fxn <- function(model_name, m_fit, lookup_df, pretty_names_vec) {
  
  # 1. Identify which variables survived in this specific model run
  active_vars <- lavaan::lavNames(m_fit, type = "ov")
  
  # 2. Match long column names that lavaan used to what exists in the lookup table
  matched_nodes <- lookup_df %>%
    filter(var %in% active_vars | Lisaname %in% active_vars) %>%
    mutate(actual_col = if_else(var %in% active_vars, var, Lisaname)) %>%
    select(SEMnode, actual_col) %>%
    distinct()
  
  cols_to_extract <- matched_nodes$actual_col
  
  # 3. Always use the master dataset guild.dfas1
  raw_data <- guild.dfas1
  
  # Safety check: ensure year column exists
  if (!"year" %in% colnames(raw_data)) {
    stop("Could not find a 'year' column in guild.dfas1. Please adjust the time index column name.")
  }
  
  # 4. Truncate the rows dynamically if it's a "short" model configuration
  # Slices rows to capture only what went into the short model (adjust years if needed)
  if (grepl("short", model_name)) {
    # If you have a separate object holding your exact short boundary dates, you can use that.
    # Otherwise, this safely drops any rows containing NA values in our active variables:
    raw_data <- raw_data %>%
      filter(across(all_of(intersect(cols_to_extract, colnames(raw_data))), ~ !is.na(.)))
  }
  
  cols_present <- intersect(cols_to_extract, colnames(raw_data))
  
  # 5. Filter data columns and pivot to a long layout
  ts_data_long <- raw_data %>%
    select(year, all_of(cols_present)) %>%
    pivot_longer(cols = -year, names_to = "raw_variable", values_to = "value") %>%
    mutate(pretty_label = ifelse(raw_variable %in% names(pretty_names_vec), 
                                 pretty_names_vec[raw_variable], 
                                 raw_variable))
  
  # Determine year range for the subtitle display
  year_range <- range(ts_data_long$year, na.rm = TRUE)
  subtitle_text <- paste0("Time Frame: ", year_range[1], " - ", year_range[2], 
                          " (n = ", length(unique(ts_data_long$year)), " years)")
  
  # 6. Build the facetted time series plot grid
  p_ts <- ggplot(ts_data_long, aes(x = year, y = value)) +
    geom_line(color = "steelblue", size = 1.2) +
    geom_point(color = "darkblue", size = 1.8) +
    facet_wrap(~pretty_label, scales = "free_y", ncol = 2) +
    theme_minimal(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold", color = "black"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey80"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40")
    ) +
    labs(
      title = paste("Raw Data Time Series:", model_name),
      subtitle = subtitle_text,
      x = "Year",
      y = "Value"
    )
  
  return(p_ts)
}
