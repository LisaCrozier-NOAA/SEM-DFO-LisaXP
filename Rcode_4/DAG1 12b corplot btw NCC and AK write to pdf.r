library(tidyverse)
library(reshape2)

# 1. Clean and rank your master list of active variables
guild_19yr <- guild.dfas1 %>% filter(year > 2002)

sar_correlations <- guild_19yr %>%
  select(all_of(var_lookup_NCC_AK$Lisaname), X16_SAR) %>%
  cor(use = "complete.obs") %>%
  .[, "X16_SAR"] %>%
  enframe(name = "Lisaname", value = "SAR_Correlation")

ranked_variables <- var_lookup_NCC_AK %>% 
  filter(Lisaname %in% colnames(guild_19yr)) %>% 
  left_join(sar_correlations, by = "Lisaname") %>%
  group_by(SEMnode) %>%
  arrange(desc(SAR_Correlation), .by_groups = TRUE) %>%
  ungroup()

# 2. Extract specific character vectors of the sorted variables per guild
vars_preyNCC <- ranked_variables %>% filter(SEMnode == "PreyNCC") %>% pull(Lisaname)
vars_predNCC <- ranked_variables %>% filter(SEMnode == "PredNCC") %>% pull(Lisaname)
vars_preyAK  <- ranked_variables %>% filter(SEMnode == "PreyAK")  %>% pull(Lisaname)
vars_predAK  <- ranked_variables %>% filter(SEMnode == "PredAK")  %>% pull(Lisaname)

# 3. Generate master melted correlation matrix
matrix_data     <- guild_19yr %>% select(all_of(ranked_variables$Lisaname))
cor_matrix_full <- cor(matrix_data, use = "complete.obs")
melted_master   <- melt(cor_matrix_full) %>%
  rename(Var1 = Var1, Var2 = Var2, Correlation = value)


plot_matrix_page <- function(data_long, x_vars, y_vars, title_text) {
  # Filter master matrix to just the targeted page view
  page_data <- data_long %>%
    filter(Var1 %in% x_vars & Var2 %in% y_vars) %>%
    mutate(
      Var1 = factor(Var1, levels = x_vars),
      Var2 = factor(Var2, levels = y_vars)
    )
  
  p <- ggplot(page_data, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white", size = 0.3) +
    scale_fill_gradient2(
      low = "#0571b0", mid = "white", high = "#ca0020", 
      midpoint = 0, limit = c(-1, 1), name = "Pearson\nCorrelation"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14),
      panel.grid = element_blank()
    ) +
    labs(title = title_text, x = "", y = "") +
    coord_fixed()
  
  print(p)
}


# Open a high-resolution, wide-format PDF canvas
pdf("LisaXP/outputs_4/Guild_Correlation_Atlas.pdf", width = 12, height = 10)

# ------------------------------------------------------------------------------
# PAGE 1: Within-Guild View - The Massive PreyNCC Block
# ------------------------------------------------------------------------------
plot_matrix_page(
  melted_master, x_vars = vars_preyNCC, y_vars = vars_preyNCC,
  title_text = "Page 1: Internal PreyNCC Correlations (Sorted by SAR Correlation)"
)

# ------------------------------------------------------------------------------
# PAGE 2: Within-Guild View - The Massive PredAK Block
# ------------------------------------------------------------------------------
plot_matrix_page(
  melted_master, x_vars = vars_predAK, y_vars = vars_predAK,
  title_text = "Page 2: Internal PredAK Correlations (Sorted by SAR Correlation)"
)

# ------------------------------------------------------------------------------
# PAGE 3: Regional Handshake - PreyNCC vs PredAK (The two giants competing)
# ------------------------------------------------------------------------------
plot_matrix_page(
  melted_master, x_vars = vars_preyNCC, y_vars = vars_predAK,
  title_text = "Page 3: Cross-Regional Teleconnections (PreyNCC vs. PredAK)"
)

# ------------------------------------------------------------------------------
# PAGE 4: Local NCC Dynamics - PreyNCC vs PredNCC
# ------------------------------------------------------------------------------
plot_matrix_page(
  melted_master, x_vars = vars_preyNCC, y_vars = vars_predNCC,
  title_text = "Page 4: Local Northern California Current Dynamics (PreyNCC vs. PredNCC)"
)

# ------------------------------------------------------------------------------
# PAGE 5: The Remaining Groups (PredNCC, PreyAK, PredAK inter-correlations)
# ------------------------------------------------------------------------------
remaining_vars <- c(vars_predNCC, vars_preyAK)
plot_matrix_page(
  melted_master, x_vars = remaining_vars, y_vars = vars_predAK,
  title_text = "Page 5: Regional Predators Matrix (PredNCC & PreyAK vs. PredAK)"
)

# Close the file device and lock in your multi-page PDF
dev.off()
