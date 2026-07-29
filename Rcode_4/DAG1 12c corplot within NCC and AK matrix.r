library(tidyverse)
library(reshape2)

# 1. Slice data to the 19 complete years (post-2002)
guild_19yr <- guild.dfas1 %>% 
  filter(year > 2002)

# 2. Calculate every variable's raw correlation with X16_SAR for unified sorting
sar_correlations <- guild_19yr %>%
  select(all_of(colnames(guild_19yr))) %>% 
  cor(use = "complete.obs") %>%
  .[, "X16_SAR"] %>%
  enframe(name = "Lisaname", value = "SAR_Correlation")

# 3. Create your master ranking map (Preserving ALL columns)
ranked_variables <- guild.dfas1 %>%
  colnames() %>%
  enframe(name = NULL, value = "Lisaname") %>%
  filter(Lisaname != "year") %>% 
  
  # Join your guild lookup table
  left_join(var_lookup_NCC_AK %>% select(Lisaname, SEMnode), by = "Lisaname") %>%
  
  # Keep the original "NA" label for variables not in the lookup table
  mutate(SEMnode = if_else(is.na(SEMnode), "NA", SEMnode)) %>%
  
  # Rank everything within its group by its relationship to SAR
  left_join(sar_correlations, by = "Lisaname") %>%
  group_by(SEMnode) %>%
  arrange(desc(SAR_Correlation), .by_groups = TRUE) %>%
  ungroup() %>%
  mutate(Display_Label = paste0("[", SEMnode, "] ", Lisaname))

# 4. Compute full pairwise correlation matrix
matrix_data <- guild_19yr %>% select(all_of(ranked_variables$Lisaname))
cor_matrix_full <- cor(matrix_data, use = "complete.obs")

# 5. Melt the matrix completely (No filtering out of any rows or columns)
melted_master_uncropped <- melt(cor_matrix_full) %>%
  rename(Var1 = Var1, Var2 = Var2, Correlation = value) %>%
  
  left_join(ranked_variables %>% select(Lisaname, Display_Label, SEMnode), by = c("Var1" = "Lisaname")) %>%
  rename(Label_Var1 = Display_Label, Guild1 = SEMnode) %>%
  
  left_join(ranked_variables %>% select(Lisaname, Display_Label, SEMnode), by = c("Var2" = "Lisaname")) %>%
  rename(Label_Var2 = Display_Label, Guild2 = SEMnode) %>%
  
  # Lock in the factor levels so the ordering stays perfectly aligned
  mutate(
    Label_Var1 = factor(Label_Var1, levels = ranked_variables$Display_Label),
    Label_Var2 = factor(Label_Var2, levels = ranked_variables$Display_Label),
    Guild1 = factor(Guild1, levels = c("PreyNCC", "PredNCC", "PreyAK", "PredAK", "NA")),
    Guild2 = factor(Guild2, levels = c("PreyNCC", "PredNCC", "PreyAK", "PredAK", "NA"))
  )


# 6. Generate the full plot
master_matrix_corplot<-
ggplot(melted_master_uncropped, aes(x = Label_Var1, y = Label_Var2, fill = Correlation)) +
  geom_tile(color = "white", size = 0.1) +
  
  # Symmetrical faceting across all 5 groups (including the original NA block)
  facet_grid(Guild2 ~ Guild1, scales = "free", space = "free") +
  
  scale_fill_gradient2(
    low = "#0571b0", 
    mid = "white", 
    high = "#ca0020", 
    midpoint = 0, 
    limit = c(-1, 1), 
    name = "Pearson\nCorrelation"
  ) +
  theme_bw(base_size = 10) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6), 
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid = element_blank(),
    panel.spacing = unit(0.3, "lines")
  ) +
  labs(
    title = "Complete Symmetrical Ecosystem Correlation Matrix (Post-2002 Window)",
    subtitle = "All 70 variables mapped symmetrically. Salmon metrics sit in the 'NA' row/column panels.",
    x = "",
    y = ""
  )

# 7. Open an ultra-large PDF canvas so all 70x70 tiles are legible
pdf("LisaXP/outputs_4/Master_Matrix.pdf", width = 22, height = 20)
  master_matrix_corplot

dev.off()

ggsave(file="LisaXP/outputs_4/master_matrix_corplot.png",master_matrix_corplot)

#write csv file-----------
# 1. Cast the symmetric matrix into a standard data frame
cor_dataframe <- as.data.frame(cor_matrix_full)

# 2. Prepare the row metadata to match your PDF labels
row_metadata <- ranked_variables %>%
  select(Lisaname, SEMnode, SAR_Correlation) %>%
  rename(Variable_Name = Lisaname, Guild_Group = SEMnode)

# 3. Bind the metadata to the front of the matrix so it's easy to filter in Excel/CSV
final_csv_matrix <- bind_cols(row_metadata, cor_dataframe)

final_csv_matrix[1:5,1:5]

# 4. Write to your working directory
write.csv(
  final_csv_matrix, 
  file = "LisaXP/outputs_4/Master_Matrix.csv", 
  row.names = FALSE
)