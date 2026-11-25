
# ==============================================================================
# Script: 02_analysis/06f_plot_core_model.R
# ==============================================================================

# ==============================================================================
# Step 3: List Complete (0-NA) Indicators & Plot Backbone 1B via Lisa_sem_graph.v3
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(tidySEM)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. LIST ALL COLUMNS WITH ZERO MISSING DATA (1998-2021)
# ------------------------------------------------------------------------------
complete_cols <- sem_data %>%
  select(-year) %>%
  keep(~ sum(is.na(.x)) == 0) %>%
  names()

cat("Columns with 0 missing years (", length(complete_cols), "total ):\n")
print(complete_cols)

# ------------------------------------------------------------------------------
# 2. FIT BACKBONE 1B WITH fixed.x = TRUE (npar = 12)
# ------------------------------------------------------------------------------
model_1B <- "
  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 X12_DFA_EGOA_DFA_lowerTrophic_smoltyr + 
                                 sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z
"

fit_1B <- sem(
  model_1B,
  data         = sem_data,
  missing      = "ML",
  fixed.x      = TRUE,
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

# ------------------------------------------------------------------------------
# 3. BUILD YOUR 5x5 BACKBONE LAYOUT & PRETTY LABELS
# ------------------------------------------------------------------------------
preyNCC   <- "X02_DFA_SeaNettle_smoltyr"
predNCC   <- "X09_DFA_HakeAge5Plus_smoltyr"
growth    <- "X06_DFA_IGF_mu_smoltyr"
abundance <- "X07_cpue_IntSprJunHW_smoltyr"
preyAK    <- "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"
predAK    <- "sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z"
sar       <- "X16_DFA_sar_UCsprCh_smoltyr"

layout_1B <- matrix(c(
  preyNCC, ""    , "" , ""       , predNCC,
  ""     , growth, "" , abundance, ""     ,
  ""     , ""    , sar, ""       , ""     ,
  ""     , ""    , "" , ""       , ""     ,
  preyAK , ""    , "" , ""       , predAK
), nrow = 5, ncol = 5, byrow = TRUE)

labels_1B <- c(
  "X02_DFA_SeaNettle_smoltyr"                         = "NCC Prey\n(Sea Nettle)",
  "X09_DFA_HakeAge5Plus_smoltyr"                      = "NCC Pred\n(Hake 5+)",
  "X06_DFA_IGF_mu_smoltyr"                            = "Early Marine\nGrowth",
  "X07_cpue_IntSprJunHW_smoltyr"                      = "Smolt\nCPUE",
  "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"             = "AK Prey\n(EGOA Lower)",
  "sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z" = "AK Pred\n(Sablefish Peff)",
  "X16_DFA_sar_UCsprCh_smoltyr"                       = "Chinook\nSAR"
)

# ------------------------------------------------------------------------------
# 4. PLOT BACKBONE 1B USING Lisa_sem_graph.v3
# ------------------------------------------------------------------------------
plot_1B <- Lisa_sem_graph.v3(
  fit               = fit_1B,
  title             = "Backbone Model 1B (Structural Core)",
  layout            = layout_1B,
  node_labels       = labels_1B,
  show_coefficients = TRUE,
  show_r2           = TRUE,
  r2_nodes          = c(growth, abundance, sar),
  width             = 12,
  height            = 8,
  save              = FALSE,
  display           = TRUE
)

