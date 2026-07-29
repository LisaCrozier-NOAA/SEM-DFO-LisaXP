library(tidyverse)
library(reshape2)

#Step 1: Calculate SAR Correlations and Rank Your Variables-------

# 1. Slice data to the 19 complete years (post-2002)
guild_19yr <- guild.dfas1 %>% 
  filter(year > 2002)

# 2. Get clean name mappings from your lookup table
name_mapping <- var_lookup_NCC_AK %>% 
  filter(Lisaname %in% colnames(guild_19yr)) %>% 
  select(Lisaname, SEMnode)

# 3. Calculate each variable's raw correlation with X16_SAR
sar_correlations <- guild_19yr %>%
  select(all_of(name_mapping$Lisaname), X16_SAR) %>%
  cor(use = "complete.obs") %>%
  .[, "X16_SAR"] %>%
  enframe(name = "Lisaname", value = "SAR_Correlation")

# 4. Rank variables strictly within their Guild groups based on their SAR correlation
ranked_variables <- name_mapping %>%
  left_join(sar_correlations, by = "Lisaname") %>%
  group_by(SEMnode) %>%
  # Arrange descending (highest positive to lowest/negative correlation with SAR)
  arrange(desc(SAR_Correlation), .by_groups = TRUE) %>%
  ungroup() %>%
  mutate(Display_Label = paste0("[", SEMnode, "] ", Lisaname))

#Step 2: Inject the Empty Spacers into the Axis Ordering---------
# --- Define the ordered X-Axis (NCC Variables) ---
ncc_prey <- ranked_variables %>% filter(SEMnode == "PreyNCC") %>% pull(Display_Label)
ncc_pred <- ranked_variables %>% filter(SEMnode == "PredNCC") %>% pull(Display_Label)

# Place an empty space character right in the middle
ncc_order_with_gap <- c(ncc_prey, " ", ncc_pred)

# --- Define the ordered Y-Axis (BC/AK Variables) ---
ak_prey <- ranked_variables %>% filter(SEMnode == "PreyAK") %>% pull(Display_Label)
ak_pred <- ranked_variables %>% filter(SEMnode == "PredAK") %>% pull(Display_Label)

# Place an empty space character right in the middle
ak_order_with_gap <- c(ak_prey, "  ", ak_pred) # Using 2 spaces to keep the factor levels distinct

#Step 3: Compute Cross-Correlation Matrix and Filter---------
# Calculate the active matrix data
matrix_data <- guild_19yr %>% select(all_of(ranked_variables$Lisaname))
cor_matrix_full <- cor(matrix_data, use = "complete.obs")

# Melt and annotate with our display labels
melted_cor_annotated <- melt(cor_matrix_full) %>%
  rename(Var1 = Var1, Var2 = Var2, Correlation = value) %>%
  left_join(ranked_variables %>% select(Lisaname, Display_Label, SEMnode), by = c("Var1" = "Lisaname")) %>%
  rename(Label_Var1 = Display_Label, Guild1 = SEMnode) %>%
  left_join(ranked_variables %>% select(Lisaname, Display_Label, SEMnode), by = c("Var2" = "Lisaname")) %>%
  rename(Label_Var2 = Display_Label, Guild2 = SEMnode)

# Filter specifically to map NCC (X-axis) against BC/AK (Y-axis)
regional_comparison_grid <- melted_cor_annotated %>%
  filter(grepl("NCC", Guild1) & grepl("AK", Guild2))



#Step 4: Plot the Heatmap with Guild Spacers----------
NCCvAK_corplot<-
  
ggplot(regional_comparison_grid, aes(x = Label_Var1, y = Label_Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  
  # Enforce the strict ranked ordering and preserve the blank gap spaces
  scale_x_discrete(limits = ncc_order_with_gap, drop = FALSE) +
  scale_y_discrete(limits = ak_order_with_gap, drop = FALSE) +
  
  scale_fill_gradient2(
    low = "#0571b0", 
    mid = "white", 
    high = "#ca0020", 
    midpoint = 0, 
    limit = c(-1, 1), 
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal(base_size = 11) + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 13)
  ) +
  labs(
    title = "Regional Teleconnections (NCC vs. BC/AK)",
    subtitle = "Variables ranked within guilds by raw correlation with X16_SAR (Post-2002 Data)",
    x = "NCC Variables (Ordered by SAR Correlation →)",
    y = "BC/AK Variables (Ordered by SAR Correlation →)"
  ) +
  coord_fixed()

#Print NCC v AK corplot-----------
print(NCCvAK_corplot)


#w/in guild corplots-------------
library(tidyverse)
library(reshape2)

# 1. Slice data to the 19 complete years (post-2002)
guild_19yr <- guild.dfas1 %>% 
  filter(year > 2002)

# 2. Extract variable names and calculate their raw correlation with X16_SAR
sar_correlations <- guild_19yr %>%
  select(all_of(var_lookup_NCC_AK$Lisaname), X16_SAR) %>%
  cor(use = "complete.obs") %>%
  .[, "X16_SAR"] %>%
  enframe(name = "Lisaname", value = "SAR_Correlation")

# 3. Create a clean ranking map grouped by SEMnode
ranked_variables <- var_lookup_NCC_AK %>% 
  filter(Lisaname %in% colnames(guild_19yr)) %>% 
  left_join(sar_correlations, by = "Lisaname") %>%
  group_by(SEMnode) %>%
  arrange(desc(SAR_Correlation), .by_groups = TRUE) %>%
  ungroup() %>%
  # Create simple labels so names don't get too long inside the facet panels
  mutate(Display_Label = Lisaname)

# 4. Generate the full pairwise correlation matrix
matrix_data <- guild_19yr %>% select(all_of(ranked_variables$Lisaname))
cor_matrix_full <- cor(matrix_data, use = "complete.obs")


# Melt the full symmetric matrix
melted_all <- melt(cor_matrix_full) %>%
  rename(Var1 = Var1, Var2 = Var2, Correlation = value) %>%
  
  # Join metadata for the X-Axis (Var1)
  left_join(ranked_variables %>% select(Lisaname, Display_Label, SEMnode), by = c("Var1" = "Lisaname")) %>%
  rename(Label_Var1 = Display_Label, Guild1 = SEMnode) %>%
  
  # Join metadata for the Y-Axis (Var2)
  left_join(ranked_variables %>% select(Lisaname, Display_Label, SEMnode), by = c("Var2" = "Lisaname")) %>%
  rename(Label_Var2 = Display_Label, Guild2 = SEMnode) %>%
  
  # Enforce the strict SAR-ranked ordering across the factor levels
  mutate(
    Label_Var1 = factor(Label_Var1, levels = ranked_variables$Display_Label),
    Label_Var2 = factor(Label_Var2, levels = ranked_variables$Display_Label),
    # Enforce a clean regional layout order for the facet blocks
    Guild1 = factor(Guild1, levels = c("PreyNCC", "PredNCC", "PreyAK", "PredAK")),
    Guild2 = factor(Guild2, levels = c("PreyNCC", "PredNCC", "PreyAK", "PredAK"))
  )


withinNCC_withinAK_corplot<-
  
ggplot(melted_all, aes(x = Label_Var1, y = Label_Var2, fill = Correlation)) +
  geom_tile(color = "white", size = 0.2) +
  
  # This creates the clean separation lines between guilds automatically
  facet_grid(Guild2 ~ Guild1, scales = "free", space = "free") +
  
  scale_fill_gradient2(
    low = "#0571b0", 
    mid = "white", 
    high = "#ca0020", 
    midpoint = 0, 
    limit = c(-1, 1), 
    name = "Pearson\nCorrelation"
  ) +
  theme_bw(base_size = 11) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid = element_blank(),
    panel.spacing = unit(0.3, "lines")
  ) +
  labs(
    title = "Symmetric Correlation Space: Within & Between SEMnodes",
    subtitle = "Diagonal blocks show internal guild consistency. Ranked by raw SAR correlation (Top/Left = Highest).",
    x = "Predictor Variables (Sorted by SAR Correlation →)",
    y = "Predictor Variables (Sorted by SAR Correlation →)"
  )


#Print withinNCC_withinAK_corplot----------
print(withinNCC_withinAK_corplot)