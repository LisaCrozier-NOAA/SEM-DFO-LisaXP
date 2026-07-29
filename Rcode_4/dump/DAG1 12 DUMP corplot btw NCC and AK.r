library(tidyverse)
library(reshape2) # For melting the matrix

guild.dfas1<-read.csv("LisaXP/outputs_4/guild.dfa.NCC.AK.csv",row.names = NULL)
var_lookup_NCC_AK<- read.csv("LisaXP/outputs_4/var_lookup_NCC_AK.SAR.csv",row.names = NULL)

head(guild.dfas1)
head(var_lookup_NCC_AK)

# 1. Prepare the 'North' data (assuming it's formatted like your Prey data)
dat_raw_north <- dat_raw %>%
  filter(SEMlatent == "north") %>%
  pivot_wider(id_cols = !SEMlatent, names_from = shortName, values_from = finalVal)

# 2. Compute the correlation matrix
# We remove the first column (the ID/index column) from both
cor_matrix <- cor(dat_raw_prey[,-1], 
                  dat_raw_north[,-1], 
                  use = "pairwise.complete.obs")

# 3. Melt the matrix into a long format for ggplot
melted_cor <- melt(cor_matrix)
colnames(melted_cor) <- c("Prey_Var", "North_Var", "Correlation")

# 4. Plot
NCCvNorth_prey_corplot<-
ggplot(melted_cor, aes(x = Prey_Var, y = North_Var, fill = Correlation)) +
  geom_tile(color = "white") +
  # Customizing the color scale (Red = Positive, Blue = Negative)
  scale_fill_gradient2(low = "#0571b0", 
                       mid = "white", 
                       high = "#ca0020", 
                       midpoint = 0, 
                       limit = c(-1, 1), 
                       space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(title = "Correlation: Prey in NCC vs. BC/AK",
       x = "Prey Variables",
       y = "North Variables") +
  coord_fixed() # Makes the tiles perfectly square
