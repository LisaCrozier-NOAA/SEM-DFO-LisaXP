
# Table 1 Template (Methods): Indicator Architecture Matrix
# Generates a publication-ready Word/HTML table summarizing all 190 indicators across 7 SEM nodes and 16 guilds.


library(gt)
library(tidyverse)

# Load metadata
indicators_meta <- read.csv("metadata/indicators.csv")

table1_gt <- indicators_meta %>%
  select(SEMlatent, guild, shortName, indicator, source) %>%
  arrange(SEMlatent, guild, shortName) %>%
  gt(groupname_col = "SEMlatent") %>%
  tab_header(
    title = md("**Table 1. Environmental and Ecological Indicators**"),
    subtitle = "Categorized by Structural Equation Model (SEM) nodes and functional guilds"
  ) %>%
  cols_label(
    guild = md("**Guild**"),
    shortName = md("**Code**"),
    indicator = md("**Description**"),
    source = md("**Data Source**")
  ) %>%
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    column_labels.font.weight = "bold"
  )

# Export options
gtsave(table1_gt, "output/tables/Table1_Indicators.html")
 gtsave(table1_gt, "output/tables/Table1_Indicators.docx") # Export directly to Word