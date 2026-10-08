# ==============================================================================
# Script: render_trophic_interaction_table.R
# Purpose: Generate a publication-ready gt table summarizing trophic 
#          interactions, spatial drivers, pollock age dynamics, and primary 
#          literature references for Spring Chinook Salmon predators.
# Output: Rcode_for_paper/figures/trophic_interaction_table.png (and .html)
# ==============================================================================

suppressPackageStartupMessages({
  library(gt)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP OUTPUT PATHS
# ------------------------------------------------------------------------------
proj_dir    <- file.path(getwd(), "Rcode_for_paper")
figures_dir <- file.path(proj_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# STEP 1: DEFINE TABLE DATASTRUCTURE (ECOLOGICALLY CORRECTED)
# ------------------------------------------------------------------------------
interaction_data <- data.frame(
  guild = c(
    "Guild 1: Seabirds",
    "Guild 2: Demersal Fish & Mammals (NCC)",
    "Guild 3: Alaska Mammals & Small Piscivores",
    "Guild 4: Alaska Demersal Piscivores",
    "Guild 5: Apex Pelagic Sharks"
  ),
  taxa = c(
    "Common Murre, Cormorants, Gulls, Loons",
    "Pacific Hake, Harbor Seals",
    "Steller Sea Lion, Sablefish, Spiny Dogfish",
    "Arrowtooth Flounder, Pacific Halibut, Pacific Cod",
    "Salmon Shark, Sleeper Shark"
  ),
  target_stage = c(
    "Smolt (Surface outmigration)",
    "Smolt & Adult (Hake on smolts; Seals on adults)",
    "Smolt & Sub-Adult (Inshore to offshore shelf)",
    "Adult (Deep benthic & demersal shelf)",
    "Adult (Epipelagic & mesopelagic offshore)"
  ),
  spatial_driver = c(
    "High HCI contracts coastal upwelling zone; forces birds and smolts into narrow nearshore plume.",
    "Warm SST/HCI drives Hake inshore/northward and concentrates seals near river mouths.",
    "SST shifts coastal forage aggregations; Sablefish and Sea Lions shift bathymetric depth with warming.",
    "Warm GoA SST restricts benthic boundary layer; forces flatfish diel vertical migrations into midwater salmon zones.",
    "SST fronts in GoA/BSAI compress shark foraging depths directly into epipelagic adult salmon return corridors."
  ),
  prey_buffering = c(
    "High prey-switching: Abundant Anchovy, Market Squid, or Herring divert dive effort away from smolts.",
    "Hake track Krill but switch to Anchovy/Squid; Seals switch to Eulachon/Hake.",
    "SSLs & Sablefish track pelagic Age 1–2 Pollock cohorts and Herring in shallow shelf layers.",
    "Piscivorous buffering: Large Arrowtooth & Halibut switch between benthic Age 1+ Pollock cohorts and pelagic Capelin.",
    "Opportunistic apex generalists; target large Adult Pollock (Age 3+), Pink Salmon, and Herring during salmon runs."
  ),
  #working ref list, 1/line for easier review
  references = c(
    "Hostetter et al. (2012);<br> Nur et al. (2011);<br> Phillips et al. (2017, 2021);<br> Tucker et al. (2016);<br> Zamon et al. (2014, 2024)",
    "Emmett & Sampson (2007);<br> Olesiuk (1993);<br> Sabal et al. (2023);<br> Shirk et al. (2023);<br> Steingass (2017);<br> Thomas et al. (2017);<br> Wright et al. (2010)",
    "Bredesen et al. (2004);<br> Coutré et al. (2015);<br> Daly et al. (2024);<br> Marston et al. (2002);<br> Sigler et al. (2009);<br> Sinclair & Zeppelin (2002);<br> Tollit et al. (2015)",
    "Bizzarro et al. (2023);<br> Knoth & Foy (2008);<br> Langan et al. (2024);<br> Wells et al. (2024, 2025);<br> Yang (1993, 2011)",
    "Carlisle et al. (2015);<br> Hulbert et al. (2005);<br> Nagasawa (1998);<br> Riddell et al. (2018);<br> Seitz et al. (2019);<br> Sigler et al. (2006)"
  ),
  #original wrapping refs for final paper
  # references = c(
  #   "Hostetter et al. (2012); Nur et al. (2011); Phillips et al. (2017, 2021); Tucker et al. (2016); Zamon et al. (2014, 2024)",
  #   "Emmett & Sampson (2007); Olesiuk (1993); Sabal et al. (2023); Shirk et al. (2023); Steingass (2017); Thomas et al. (2017); Wright et al. (2010)",
  #   "Bredesen et al. (2004); Coutré et al. (2015); Daly et al. (2024); Marston et al. (2002); Sigler et al. (2009); Sinclair & Zeppelin (2002); Tollit et al. (2015)",
  #   "Bizzarro et al. (2023); Knoth & Foy (2008); Langan et al. (2024); Wells et al. (2024, 2025); Yang (1993, 2011)",
  #   "Carlisle et al. (2015); Hulbert et al. (2005); Nagasawa (1998); Riddell et al. (2018); Seitz et al. (2019); Sigler et al. (2006)"
  # ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# STEP 2: BUILD GT TABLE WITH MANUSCRIPT STYLING
# ------------------------------------------------------------------------------
gt_table <- interaction_data %>%
  gt() %>%
  tab_header(
    title = md("**Predator Guild Interaction Mechanics & Trophic Justification**"),
    subtitle = md("*Ecological mechanisms governing Spring Chinook mortality, temperature-driven spatial overlap, and age-structured prey buffering across migration regions*")
  ) %>%
  cols_label(
    guild          = md("**Functional Guild**"),
    taxa           = md("**Taxonomic Representatives**"),
    target_stage   = md("**Target Life Stage**"),
    spatial_driver = md("**Spatial Overlap Driver (Warm Years)**"),
    prey_buffering = md("**Prey-Switching & Buffering Mechanism**"),
    references     = md("**Primary Scientific References**")
  ) %>%
  cols_align(align = "left", columns = everything()) %>%
  cols_width(
    guild          ~ px(100),
    taxa           ~ px(100),
    target_stage   ~ px(100),
    spatial_driver ~ px(140),
    prey_buffering ~ px(140),
    references     ~ px(140)
    # guild          ~ px(140),
    # taxa           ~ px(150),
    # target_stage   ~ px(140),
    # spatial_driver ~ px(240),
    # prey_buffering ~ px(240),
    # references     ~ px(200)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#0F172A"),
      cell_text(color = "white", weight = "bold", size = px(14))
    ),
    locations = cells_title(groups = "title")
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#1E293B"),
      cell_text(color = "#E2E8F0", style = "italic", size = px(11))
    ),
    locations = cells_title(groups = "subtitle")
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#334155"),
      cell_text(color = "white", weight = "bold", size = px(11))
    ),
    locations = cells_column_labels()
  ) %>%
  opt_row_striping(row_striping = TRUE) %>%
  tab_style(
    style = cell_fill(color = "#F8FAFC"),
    locations = cells_body(rows = seq(1, nrow(interaction_data), 2))
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", color = "#0F172A", size = px(10.5)),
    locations = cells_body(columns = guild)
  ) %>%
  tab_style(
    style = cell_text(size = px(9.5), color = "#1E293B"),
    locations = cells_body(columns = c(taxa, target_stage, spatial_driver, prey_buffering, references))
  ) %>%
  tab_options(
    table.border.top.style = "none",
    table.border.bottom.style = "solid",
    table.border.bottom.width = px(2),
    table.border.bottom.color = "#0F172A",
    column_labels.border.top.color = "#0F172A",
    column_labels.border.bottom.color = "#0F172A",
    table_body.hlines.color = "#E2E8F0",
    table_body.border.bottom.color = "#0F172A",
    data_row.padding = px(8)
  )%>%
  fmt_markdown(columns = references)

# Print in RStudio Viewer
print(gt_table)

# ------------------------------------------------------------------------------
# STEP 3: EXPORT TABLE TO PNG AND HTML
# ------------------------------------------------------------------------------
html_out <- file.path(figures_dir, "trophic_interaction_table.html")
png_out  <- file.path(figures_dir, "trophic_interaction_table.png")

gtsave(gt_table, filename = html_out)

tryCatch({
  gtsave(gt_table, filename = png_out, vwidth = 1200, vheight = 800)
  cat(sprintf("\nTable image successfully saved to: %s\n", png_out))
}, error = function(e) {
  cat(sprintf("\nHTML table saved to: %s\nNote: Install 'webshot2' package to enable PNG rendering via gtsave().\n", html_out))
})