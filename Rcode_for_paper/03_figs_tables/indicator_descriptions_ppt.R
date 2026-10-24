# ==============================================================================
# Script: 05_conference_indicator_summary_table_ppt.R
# Purpose: Slide-Optimized Large-Font Presentation Table for SEM Indicator Nodes
# Output: Rcode_for_paper/Routput_for_paper/tables/
# ==============================================================================

library(tidyverse)
library(gt)

# ------------------------------------------------------------------------------
# STEP 1: COMPUTE INDICATOR COUNTS & EXTRACT SUBGUILD CODES
# ------------------------------------------------------------------------------
indicator_table_data <- doug_guilds %>%
  filter(SEMnode %in% c("PreyNCC", "PredNCC", "PreyAK", "PredAK", "Growth", "Abundance")) %>%
  
  group_by(SEMnode) %>%
  summarise(
    # Total initial raw series under pastGuild + new entries
    N_Initial  = sum((!is.na(pastGuild) & pastGuild != "") | (!is.na(Lisa_21oct26) & Lisa_21oct26 != "")),
    
    # Final retained series in Lisa_21oct26
    N_Screened = sum(!is.na(Lisa_21oct26) & Lisa_21oct26 != ""),
    
    # Column 4: Extract unique active subguild codes formatted with HTML line breaks
    Subguild_Codes_HTML = paste(
      unique(Lisa_21oct26[Lisa_21oct26 != "" & !is.na(Lisa_21oct26)]), 
      collapse = "<br>"
    ),
    .groups = "drop"
  ) %>%
  
  mutate(
    # Column 1: Full SEM Node Name (Your Version)
    SEM_Node_Full = case_when(
      SEMnode == "PreyNCC"   ~ "Northern California Current<br><b>NCC Prey</b>",
      SEMnode == "PredNCC"   ~ "Northern California Current<br><b>NCC Predators</b>",
      SEMnode == "PreyAK"    ~ "BC Coast, Gulf of Alaska<br><b>AK Prey</b>",
      SEMnode == "PredAK"    ~ "BC Coast, Gulf of Alaska<br><b>AK Predators</b>",
      SEMnode == "Growth"    ~ "<b>Salmon Growth / Condition</b>",
      SEMnode == "Abundance" ~ "<b>Salmon Early Abundance (CPUE)</b>",
      TRUE                   ~ SEMnode
    ),
    
    # Column 2: Decoded Subguild Descriptions (Your Version)
    Subguild_Names = case_when(
      SEMnode == "PreyNCC"   ~ "Zooplankton (Newport, WA & WCVI surveys), Competitors (Juv Sablefish, Market squid, Jellies), \nKrill, Juv Forage Fish",
      SEMnode == "PredNCC"   ~ "Piscivorous Seabirds, Marine Mammals, Nearshore/Surface-feeding Groundfish",
      SEMnode == "PreyAK"    ~ "Low & Mid-Trophic indices (Ferris et al, EGOA/WGOA), WCVI Krill, Forage Fish, Seabird diets/reprod success, Pink Salmon Competitors",
      SEMnode == "PredAK"    ~ "Predatory Groundfish & Sharks (Dogfish, Sleeper & Salmon Sharks), Marine Mammals",
      SEMnode == "Growth"    ~ "Juvenile Salmon IGF-1 Condition & Stomach Fullness",
      SEMnode == "Abundance" ~ "Juvenile Salmon CPUE Indices",
      TRUE                   ~ ""
    ),
    
    # Column 3: Screened Indicator Counts (Your Version)
    N_Counts = N_Initial# paste0("<b>", N_Initial, " → ", N_Screened, "</b>")
  ) %>%
  
  # Ordering for final layout
  mutate(
    order = case_when(
      SEMnode == "PreyNCC"   ~ 1,
      SEMnode == "PredNCC"   ~ 2,
      SEMnode == "PreyAK"    ~ 3,
      SEMnode == "PredAK"    ~ 4,
      SEMnode == "Growth"    ~ 5,
      SEMnode == "Abundance" ~ 6,
      TRUE                   ~ 99
    )
  ) %>%
  arrange(order)

# ------------------------------------------------------------------------------
# STEP 2: BUILD LARGE-FONT SLIDE GT TABLE
# ------------------------------------------------------------------------------
indicator_summary_gt <- indicator_table_data %>%
  select(SEM_Node_Full, Subguild_Names, N_Counts, Subguild_Codes_HTML) %>%
  gt() %>%
  tab_header(
    title = md("**Ecosystem Indicator Screening & Subguild Consolidation**")
  ) %>%
  cols_label(
    SEM_Node_Full       = md("**SEM Node**"),
    Subguild_Names      = md("**Taxa**"),
    N_Counts            = md("**Time Series (N)**"), #<br>(Raw → Final)
    Subguild_Codes_HTML = md("**Latent Nodes**")
  ) %>%
  fmt_markdown(columns = c(SEM_Node_Full, N_Counts, Subguild_Codes_HTML)) %>%
  cols_align(align = "left", columns = c(SEM_Node_Full, Subguild_Names, Subguild_Codes_HTML)) %>%
  cols_align(align = "center", columns = c(N_Counts)) %>%
  
  # Column width distribution tailored for 16:9 slides
  cols_width(
    SEM_Node_Full       ~ px(210),
    Subguild_Names      ~ px(340),
    N_Counts            ~ px(120),
    Subguild_Codes_HTML ~ px(230)
  ) %>%
  
  # Presentation-Ready Large Styling
  tab_options(
    table.width                = px(900),
    table.font.size            = px(15),         # Large body font
    heading.title.font.size    = px(20),         # Main slide title font
    column_labels.font.size    = px(17),         # Header font
    column_labels.font.weight  = "bold",
    data_row.padding           = px(4),          # Compact padding to fit slide height
    table.border.top.style     = "solid",
    table.border.bottom.style  = "solid"
  )

# Print and Export Slide Image
print(indicator_summary_gt)

gtsave(indicator_summary_gt, file.path(tbl_out_dir, "Table_Indicator_Summary_Large_PPT.html"))
gtsave(
  data     = indicator_summary_gt,
  filename = file.path(tbl_out_dir, "Table_Indicator_Summary_Large_PPT.png"),
  vwidth   = 950,
  vheight  = 550
)

cat("Large-font slide table exported to: ", file.path(tbl_out_dir, "Table_Indicator_Summary_Large_PPT.png"), "\n")