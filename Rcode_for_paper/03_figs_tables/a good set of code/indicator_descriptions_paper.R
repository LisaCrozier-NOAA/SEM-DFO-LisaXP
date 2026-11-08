# ==============================================================================
# Script: 05_conference_indicator_summary_table_ppt.R
# Purpose: Presentation Table Summarizing Indicator Filtering & Subguild Codes
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
  
  # Apply your specific modifications for Columns 1, 2, and 3
  mutate(
    # Column 1: Full SEM Node Name (Your Version)
    SEM_Node_Full = case_when(
      SEMnode == "PreyNCC"   ~ "Northern California Current<br>NCC Prey",
      SEMnode == "PredNCC"   ~ "Northern California Current<br>NCC Predators",
      SEMnode == "PreyAK"    ~ "BC Coast, Gulf of Alaska<br>AK Prey",
      SEMnode == "PredAK"    ~ "BC Coast, Gulf of Alaska<br>AK Predators",
      SEMnode == "Growth"    ~ "Salmon Growth / Condition",
      SEMnode == "Abundance" ~ "Salmon Early Abundance (CPUE)",
      TRUE                   ~ SEMnode
    ),
    
    # Column 2: Decoded Subguild Descriptions (Your Version)
    Subguild_Names = case_when(
      SEMnode == "PreyNCC"   ~ "Zooplankton (Newport, WA & WCVI surveys), Krill, Forage Fish, Competitors (Juv Sablefish, Market squid, Jellies)",
      SEMnode == "PredNCC"   ~ "Piscivorous Seabirds, Nearshore/Surface-feeding Groundfish, Marine Mammals",
      SEMnode == "PreyAK"    ~ "Lower & Mid-Trophic indices (Ferris et al, EGOA/WGOA), WCVI Krill, Forage Fish & Seabirds, Pink Salmon Competitors",
      SEMnode == "PredAK"    ~ "Predatory Groundfish & Sharks (Dogfish, Sleeper & Salmon Sharks), Marine Mammals",
      SEMnode == "Growth"    ~ "Juvenile Salmon IGF-1 Condition & Stomach Fullness",
      SEMnode == "Abundance" ~ "Juvenile Salmon CPUE Indices",
      TRUE                   ~ ""
    ),
    
    # Column 3: Screened Indicator Counts (Your Version)
    N_Counts = paste0(N_Initial, " → ", N_Screened)
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
# STEP 2: BUILD PRESENTATION-READY GT TABLE
# ------------------------------------------------------------------------------
indicator_summary_gt <- indicator_table_data %>%
  select(SEM_Node_Full, Subguild_Names, N_Counts, Subguild_Codes_HTML) %>%
  gt() %>%
  tab_header(
    title = md("**Summary of Ecosystem Indicator Screening & Subguild Consolidation**"),
    subtitle = "Systematic reduction from raw environmental series to Structural Equation Model (SEM) latent nodes"
  ) %>%
  cols_label(
    SEM_Node_Full       = md("**SEM Node**"),
    Subguild_Names      = md("**Consolidated Subguilds**"),
    N_Counts            = md("**Series Count<br>(Raw → Screened)**"),
    Subguild_Codes_HTML = md("**SEM Subguild Codes<br>(Lisa_21oct26)**")
  ) %>%
  fmt_markdown(columns = c(SEM_Node_Full, Subguild_Codes_HTML)) %>%
  cols_align(align = "left", columns = c(SEM_Node_Full, Subguild_Names, Subguild_Codes_HTML)) %>%
  cols_align(align = "center", columns = c(N_Counts)) %>%
  cols_width(
    SEM_Node_Full       ~ px(180),
    Subguild_Names      ~ px(260),
    N_Counts            ~ px(110),
    Subguild_Codes_HTML ~ px(200)
  ) %>%
  tab_options(
    table.width                = px(750),
    table.font.size            = px(11),
    heading.title.font.size    = px(13),
    heading.subtitle.font.size = px(10),
    column_labels.font.weight  = "bold",
    data_row.padding           = px(6)
  )

# Print and Export
print(indicator_summary_gt)

gtsave(indicator_summary_gt, file.path(tbl_out_dir, "Table_Indicator_Screening_Summary.html"))
gtsave(
  data     = indicator_summary_gt,
  filename = file.path(tbl_out_dir, "Table_Indicator_Screening_Summary.png"),
  vwidth   = 800,
  vheight  = 500
)