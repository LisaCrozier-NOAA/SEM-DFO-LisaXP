```{r}
#| label: fig-indicator_source_guildfile
#| fig-cap: "Supplemental Table S1: Full indicator list."
#| fig-align: center
#| fig-pos: "H"
#| out-width: "90%"
knitr::include_graphics("Routput_for_paper/Supplements/guilds.csv")
```

path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/Routput_for_paper/Supplements/guilds.csv"

guilds<-read.csv(path,stringsAsFactors = TRUE,row.names = NULL, na.strings = c("","NA")) %>%
  filter(!is.na(Lisa_15Oct2026 )) %>%
  mutate(Reference = source_short) %>%
  select(ID,SEMnode,region,subregion,Lisa_15Oct2026,indicator,source_short,Reference,source_raw) %>%
  arrange(Lisa_15Oct2026);head(guilds)

nrow(guilds)
summary(guilds)
guilds %>% filter(is.na(region))
# ID SEMnode region subregion        Lisa_15Oct2026                        indicator source_short Reference                                              source_raw
# 1 226  PredAK   <NA>      <NA> 10.PredMammalNCC_BCAK Killer.whales.NR.BC_2026_adultyr          DFO       DFO Cam Freshwater and Thomas; see 25jun26 e-mail from Lisa

guilds %>% filter(is.na(source_short))
# ID SEMnode region subregion Lisa_15Oct2026   indicator source_short Reference                             source_raw
# 1 184 PredNCC     CC       NCC 09.PredFishNCC     lingcod         <NA>      <NA> pfmc-assessments/nwfscSurvey R package
# 2 330  PreyAK     AK      EGoA     21.Climate SST Jun-Aug         <NA>      <NA>                                   <NA>
#   


library(dplyr)
library(gt)
library(webshot2)

# 1. Define paths and import data
path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/Routput_for_paper/Supplements/guilds.csv"
out_dir <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/Outputs/tables/"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

guilds <- read.csv(path, stringsAsFactors = FALSE, row.names = NULL, na.strings = c("", "NA")) %>%
  filter(!is.na(Lisa_15Oct2026)) %>%
  mutate(Reference = source_short) %>%
  select(ID, SEMnode, region, subregion, Lisa_15Oct2026, indicator, source_short, Reference, source_raw) %>%
  arrange(Lisa_15Oct2026)

# 2. Setup pagination
rows_per_page <- 25
total_rows <- nrow(guilds)
total_pages <- ceiling(total_rows / rows_per_page)

# 3. Loop through slices and save PNGs
for (i in seq_len(total_pages)) {
  start_idx <- ((i - 1) * rows_per_page) + 1
  end_idx <- min(i * rows_per_page, total_rows)
  
  df_page <- guilds[start_idx:end_idx, ]
  
  gt_table <- df_page %>%
    gt() %>%
    # Convert NAs to clean em-dashes for visual clarity
    sub_missing(columns = everything(), missing_text = "—") %>%
    tab_header(
      title = md("**Supplementary Table S1. Guild and SEM Node Definitions**"),
      subtitle = sprintf("Page %d of %d (Rows %d–%d)", i, total_pages, start_idx, end_idx)
    ) %>%
    cols_label(
      ID             = md("**ID**"),
      SEMnode        = md("**SEM Node**"),
      region         = md("**Region**"),
      subregion      = md("**Subregion**"),
      Lisa_15Oct2026 = md("**Guild / Group**"),
      indicator      = md("**Indicator**"),
      source_short   = md("**Source**"),
      Reference      = md("**Reference**"),
      source_raw     = md("**Raw Source Details**")
    ) %>%
    # Left-align text fields, center-align short codes
    cols_align(align = "center", columns = c(ID, SEMnode, region, subregion)) %>%
    cols_align(align = "left", columns = c(Lisa_15Oct2026, indicator, source_short, Reference, source_raw)) %>%
    # Tightened columns 1-4 to give maximum width (290px) to source_raw
    cols_width(
      ID             ~ px(35),
      SEMnode        ~ px(65),
      region         ~ px(55),
      subregion      ~ px(70),
      Lisa_15Oct2026 ~ px(125),
      indicator      ~ px(165),
      source_short   ~ px(95),
      Reference      ~ px(95),
      source_raw     ~ px(290)
    ) %>%
    tab_options(
      table.font.size = px(10.5),
      table.width = px(1000),
      heading.title.font.size = px(13.5),
      heading.subtitle.font.size = px(10.5),
      column_labels.font.weight = "bold",
      column_labels.font.size = px(11),
      data_row.padding = px(2.5)
    )
  
  out_file <- file.path(out_dir, sprintf("obj0_table_S1_page_%02d.png", i))
  gtsave(gt_table, filename = out_file, vwidth = 1100, zoom = 2)
}