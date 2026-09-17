
#Prompt for copilot, example code for importance calculation and plot

# The final project for the day is to refit all possible models including the new best index for csl, 
# trying both 2 and 3 variable models.
# Then we need to compute the final variable importance scores for all variables in guild_dfas1_24yr (+csl) using aic weights.
# then plot the variable importance building off code used previously

# Calculate variable importance (see 20may25 Obsidian)-------------
resultsByClus$deltaAIC <- resultsByClus$AIC - min(resultsByClus$AIC, na.rm=T)
resultsByClus$modelWeight <- exp(-0.5*resultsByClus$deltaAIC)/sum(exp(-0.5*resultsByClus$deltaAIC), na.rm=T)

allVars <- resultsByClus %>% select(ends_with("indNames", ignore.case=T)) %>% rowwise() %>%
  mutate(vars=paste(c_across(everything()), collapse=",")) %>%
  ungroup()
resultsByClus$allVars <- allVars$vars

varImportance <- data.frame(var=unique(clusData$shortName))
for(i in 1:nrow(varImportance)) {
  v <- varImportance$var[i]
  thisResults <- resultsByClus
  thisResults$inModel <- grepl(v, thisResults$allVars)
  varImportance[i, "importance"] <- sum(thisResults[thisResults$inModel, "modelWeight"])
}
varImportance <- varImportance %>% arrange(desc(importance))

# Add indicator information if using DFAs
if(useDFA) {
  guilds <- read.csv(file.path(workingDir, "output", "DFA", "guildsDFA.csv"))
  loadings <- read_excel(file.path(workingDir, "output", "DFA", "loadings.xlsx"))
  
  for(i in 1:nrow(varImportance)) {
    if(grepl("smoothed", varImportance[i, "var"])) {
      guildName <- gsub("_smoothed", "", varImportance[i, "var"])
      varImportance[i, "indicators"] <- guilds[guilds$guild==guildName, "shortName"]
      varImportance[i, "SEMnode"] <- guilds[guilds$guild==guildName, "guildSEMnode"]
    } else if(grepl("DFA1", varImportance[i, "var"])) {
      guildName <- gsub("_DFA1", "", varImportance[i, "var"])
      thisLoadings <- loadings[loadings$guild==guildName, ]
      thisLoadings <- thisLoadings %>% arrange(desc(abs(Z_est)))
      varImportance[i, "indicators"] <- paste0(thisLoadings$indicator, collapse=", ")
      varImportance[i, "SEMnode"] <- unique(thisLoadings$guildSEMnode)
    }
  }
}
write.csv(varImportance, file.path(outputDir, "varImportance.csv"), row.names=F)


#Plot variable importance----------

# Define colors that actually match the methodology workflow visuals
cartoon_colors <- c(
  "Foundations"    = "#A3C1AD", # Sage Green (The background of Panel 1 & 3)
  "Data/Assembly"  = "#D9E4EC", # Pale Blue (The paper/background in Panel 2)
  "Functional"     = "#5D8AA8", # Steel Blue (The "Bottom-Up" bin in Panel 3)
  "Simulation"     = "#E67E22", # Burnt Orange (The "Salmon" bin and SEM gears)
  "Dimension/DFA"  = "#7F8C8D", # Slate Gray (The DFA cylinders in Panel 5)
  "Selection/AIC"  = "#F39C12", # Goldenrod (The AIC medal/ribbon in Panel 6)
  "Latent/AK"      = "#2E4053"  # Deep Navy (The deep ocean in Panel 7)
)


#Barchart=================
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggh4x) 

# 1. DAG1C only---------
# Reshape and Refine Data
df_long <- importance_topvar %>%
  mutate(
    Trophic = fct_reorder(Trophic, order),
    Region = fct_reorder(Region, order),
    plot_names = fct_rev(fct_reorder(Lisaname, order))
  ) %>%
  pivot_longer(
    cols = c(DAG1C_long, DAG1C_short),
    #    cols = c(DAG1A_long, DAG1B_long, DAG1A_short, DAG1B_short),
    names_to = c("variable", "type"),
    names_sep = "_",
    values_to = "importance"
  )

# add Missing data labels
missing_labels <- df_long %>%
  group_by(plot_names, type, Trophic, Region) %>%
  summarize(all_missing = all(is.na(importance)), .groups = 'drop') %>%
  filter(all_missing == TRUE)

# 2. Create the Plot
importance_topvar_barchart<-
  ggplot(df_long, aes(x = plot_names, y = importance, fill = variable)) +
  geom_col(position = "dodge") +
  #  scale_fill_manual(values = c("DAG1A" = "#2E4053", "DAG1B" = "#D35400")) +
  scale_fill_manual(values = c("DAG1C" = "#2E4053")) +
  geom_text(data = missing_labels, 
            aes(x = plot_names, y = 0, label = "Missing data"), 
            inherit.aes = FALSE,
            hjust = -0.1, 
            color = "grey40", 
            fontface = "italic",
            size = 3) +
  coord_flip() + 
  facet_nested(Region + Trophic ~ type, 
               scales = "free_y", 
               space = "free_y",
               nest_line = element_line(color = "black"),
               resect = unit(2, "mm")) + 
  # Fix: Defined specific breaks for the vertical grid lines
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1.1), # Ensures space for labels, adjust to 1 if preferred
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal() +
  labs(
    x = "Indicator", 
    y = "Importance",
    fill = NULL 
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # Keeps only the vertical lines at our specific breaks
    panel.grid.minor.x = element_blank(), 
    strip.text.y = element_text(angle = 0),
    panel.spacing.x = unit(2, "lines"), 
    panel.spacing.y = unit(1.5, "lines"),
    strip.background = element_blank()
  )

importance_topvar_barchart

ggsave( "LisaXP/outputs_4/DAG1C_importance_topvar_barchart.png",importance_topvar_barchart)


