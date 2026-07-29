library(dplyr)
library(ggplot2)

# ==============================================================================
# 1. TRACE BACKWARDS TO GET THE CORE DFA NAME
# ==============================================================================
target_lisa_name <- "X09_DFA_HakeAge5Plus"

target_guild_info <- var_lookup_NCC_AK %>%
  filter(Lisaname == target_lisa_name)

full_guildname <- target_guild_info$guildname[1]
core_dfa_name  <- gsub("^X|\\.DFA1$|_DFA1$", "", full_guildname) # -> "09.PredFishNCC_b"

# ==============================================================================
# 2. EXTRACT RAW COMPONENT ROWS FROM THE LONG DATAFRAME (guild.data)
# ==============================================================================
# Filter guild.data down to just this DFA's components
raw_ts_long <- guild.data %>%
  filter(DFAguild == core_dfa_name) %>%
  mutate(year=year(date)) %>%
  filter(year>=1998, year<=2021) %>%
  select(year, shortName, finalVal) %>%
  rename(Variable = shortName, Value = finalVal) %>%
  mutate(Source = "Raw Input")

# ==============================================================================
# 3. EXTRACT THE FINAL OUTPUT PRODUCT FROM THE DFA DATAFRAME (guild.dfas1)
# ==============================================================================
# Isolate the final trend column from your wide dfas1 dataset and make it long
final_dfa_long <- guild.dfas1 %>%
  select(year, all_of(target_lisa_name)) %>%
  rename(Value = !!sym(target_lisa_name)) %>%
  mutate(
    Variable = target_lisa_name,
    Source = "Final DFA Output"
  )

# ==============================================================================
# 4. COMBINE AND PLOT
# ==============================================================================
# Stack the raw rows and the final product rows together
combined_plot_data <- bind_rows(raw_ts_long, final_dfa_long) %>%
  # Apply pretty human-readable labels where available
  mutate(Plot_Label = case_when(
    Source == "Final DFA Output" ~ "★ FINAL OUTPUT: Hake/Mackerel DFA",
    Variable %in% names(my_pretty_names) ~ paste(my_pretty_names[Variable]),
    TRUE ~ paste(Variable)
  ))

# Render the vertical stack of timelines
deconstruction_plot <- ggplot(combined_plot_data, aes(x = year, y = Value)) +
  geom_line(aes(color = Source), size = 1.2, show.legend = FALSE) +
  geom_point(aes(color = Source), size = 1.8, show.legend = FALSE) +
  scale_color_manual(values = c("Final DFA Output" = "firebrick", "Raw Input" = "darkgray")) +
  
  facet_wrap(~Plot_Label, scales = "free_y", ncol = 1) + 
  
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey80"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  ) +
  labs(
    title = "DFA Deconstruction: Long-Format Raw Inputs vs. Final Trend",
    x = "Year",
    y = "Value"
  )

print(deconstruction_plot)
