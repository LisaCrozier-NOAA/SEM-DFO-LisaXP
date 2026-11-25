

# ==============================================================================
# Script: 02_analysis/07g_plot_M2fish_NCCprey_NCCpred.R
# Step 1: get raw data from allData.rds
# Step 2: fit SEM
# ==============================================================================


library(tidyverse)
library(lubridate)
library(lavaan)
library(tidySEM)

# ==============================================================================
# 0. DIRECTORY & DATA SETUP
# ==============================================================================
proj_dir      <- file.path(getwd(), "Rcode_for_paper")
master_out    <- file.path(proj_dir, "Outputs")
data_out_dir  <- file.path(master_out, "data")
table_out_dir <- file.path(master_out, "tables")
fig_out_dir   <- file.path(master_out, "figures")
plot_out_dir  <- file.path(master_out, "plots")

data_in_path  <- file.path(data_out_dir, "sem_master_peff_zscored_nccpred_latent.csv")
data_out_path  <- file.path(data_out_dir,"sem_master_peff_zscored_nccpred_latent.csv")
all_data_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/allData.rds"
cw_path       <- file.path(proj_dir, "metadata/crosswalk_names.csv")
summary_fn_path       <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/functions/sem_output_summary_fxn.r"
plot_fn_path       <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/functions/Lisa_sem_graphs.v3_fxn.r"

if (file.exists(summary_fn_path)) source(summary_fn_path)
if (file.exists(plot_fn_path)) source(plot_fn_path)

all_data   <- readRDS(all_data_path)
cw_df      <- read.csv(cw_path, stringsAsFactors = FALSE)
sem_data_merged <- read.csv(data_in_path, stringsAsFactors = FALSE)


model_comparison_table <- read.csv(file.path(data_out_dir,"latent_model_comparison_table.csv"), stringsAsFactors = FALSE,,row.names=NULL)
print(model_comparison_table %>% arrange(Model) %>% filter(Model=="M2fish_preyNCC_predNCC"))          

#note that I dropped all of the akprey models because when I fleshed out all the DFAs, it was no longer significant. 
# And pink wasn't either, if I had to keep NCC prey or growth in the model

# ==============================================================================
# 1. REFIT M2fish_preyNCC_predNCC (fixing X07_cpue ~ predNCC)
# ==============================================================================
model_M2fish_preyNCC_predNCC <- "
  predAK =~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z +
            DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  preyNCC =~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr +
             X02_DFA_SeaNettle_smoltyr +
             X04_marketsquid_GAM_smoltyr +
             X05_DFA_abundSardine_smoltyr

  predNCC =~ X09_latent_HakeAge5Plus +
             X09_latent_pacificmackerel_GAM +
             X09_latent_jackmackerel_GAM +
             X09_latent_hakeFisheryWA  

  DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z ~~ 0.01 * DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z

  X06_DFA_IGF_mu_smoltyr       ~ preyNCC
  X07_cpue_IntSprJunHW_smoltyr ~ predNCC
  X16_DFA_sar_UCsprCh_smoltyr  ~ X07_cpue_IntSprJunHW_smoltyr + 
                                 predAK +
                                 X06_DFA_IGF_mu_smoltyr
"

fit_M2fish_preyNCC_predNCC <- sem(
  model_M2fish_preyNCC_predNCC,
  data         = sem_data_merged,
  std.lv       = TRUE,
  missing      = "ML",
  fixed.x      = TRUE,
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

# ==============================================================================
# 2. LOCKED 9x11 LAYOUT & LABELS (preyAK OMITTED, SPACE PRESERVED)
# ==============================================================================
layout_M2_noPreyAK <- build_master_layout(
  preyNCC_node = "preyNCC",
  predNCC_node = "predNCC",
  preyAK_node  = "",       # Dropped AK Prey (space preserved)
  predAK_node  = "predAK",
  preyNCC_ind  = c(
    "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr",
    "X02_DFA_SeaNettle_smoltyr",
    "X04_marketsquid_GAM_smoltyr",
    "X05_DFA_abundSardine_smoltyr"
  ),
  predNCC_ind  = c(
    "X09_latent_HakeAge5Plus",
    "X09_latent_pacificmackerel_GAM",
    "X09_latent_jackmackerel_GAM",
    "X09_latent_hakeFisheryWA"
  ),
  predAK_ind   = c(
    "sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z",
    "DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z"
  )
)

labels_M2_noPreyAK <- c(
  # Latent & Core Structural Nodes
  "preyNCC"                      = "NCC Prey",
  "predNCC"                      = "NCC Fish\nPredators",
  "predAK"                       = "AK Fish\nPredators",
  "X06_DFA_IGF_mu_smoltyr"       = "Growth",
  "X07_cpue_IntSprJunHW_smoltyr" = "CPUE",
  "X16_DFA_sar_UCsprCh_smoltyr"  = "SAR",
  
  # preyNCC Indicators
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr" = "Plankton\n(Jun)",
  "X02_DFA_SeaNettle_smoltyr"                 = "Sea\nNettle",
  "X04_marketsquid_GAM_smoltyr"               = "Market\nSquid",
  "X05_DFA_abundSardine_smoltyr"              = "Sardine\nHerring",
  
  # predNCC Indicators
  "X09_latent_HakeAge5Plus"        = "Hake\nAge 5+",
  "X09_latent_pacificmackerel_GAM" = "Pacific\nMackerel",
  "X09_latent_jackmackerel_GAM"    = "Jack\nMackerel",
  "X09_latent_hakeFisheryWA"       = "WA Hake\nFishery",
  
  # predAK Indicators
  "sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z"    = "Sablefish\n(3+)",
  "DFA_sleeperSharkBSAI_predAK_smoltyr__peff_TempPrey1_z" = "Sleeper\nShark"
)

# ==============================================================================
# 3. BUILD & RENDER GRAPH (BASE-R SAFE, NO SELF-LOOPS, GROWTH -> SAR TOP ENTRY)
# ==============================================================================
g <- prepare_graph(
  model           = fit_M2fish_preyNCC_predNCC,
  layout          = layout_M2_noPreyAK,
  intercepts      = FALSE,
  rect_width      = 1.8,
  rect_height     = 0.9,
  ellipses_width  = 1.9,
  ellipses_height = 1.1
)

# Apply node labels + R²
g$nodes$label <- g$nodes$name
matches <- match(g$nodes$label, names(labels_M2_noPreyAK))
g$nodes$label[!is.na(matches)] <- labels_M2_noPreyAK[matches[!is.na(matches)]]

rsq     <- inspect(fit_M2fish_preyNCC_predNCC, "r2")
r2_keep <- c("X06_DFA_IGF_mu_smoltyr", "X07_cpue_IntSprJunHW_smoltyr", "X16_DFA_sar_UCsprCh_smoltyr")
rsq_sub <- round(unlist(rsq[names(rsq) %in% r2_keep]), 2)
r2_idx  <- match(names(rsq_sub), g$nodes$name)
g$nodes$label[r2_idx[!is.na(r2_idx)]] <- paste0(
  g$nodes$label[r2_idx[!is.na(r2_idx)]], "\nR² = ", rsq_sub[!is.na(r2_idx)]
)

g$nodes$label_family <- "serif"
g$nodes$fill         <- "grey96"
g$nodes$label_fill   <- "transparent"

# Remove self-variance circles and match standardized estimates + p-values
g$edges <- g$edges[g$edges$from != g$edges$to, ]

pe <- parameterEstimates(fit_M2fish_preyNCC_predNCC, standardized = TRUE, ci = TRUE)
pe_sub <- pe[pe$op %in% c("~", "=~", "~~") & pe$lhs != pe$rhs, ]
pe_sub$from_match <- ifelse(pe_sub$op == "~", pe_sub$rhs, pe_sub$lhs)
pe_sub$to_match   <- ifelse(pe_sub$op == "~", pe_sub$lhs, pe_sub$rhs)

edge_key <- paste(g$edges$from, g$edges$op, g$edges$to)
pe_key1  <- paste(pe_sub$from_match, pe_sub$op, pe_sub$to_match)
pe_key2  <- paste(pe_sub$to_match, pe_sub$op, pe_sub$from_match)

m_idx <- match(edge_key, pe_key1)
m_idx[is.na(m_idx)] <- match(edge_key[is.na(m_idx)], pe_key2)

g$edges$std.all      <- pe_sub$std.all[m_idx]
g$edges$pvalue       <- pe_sub$pvalue[m_idx]
g$edges$linetype     <- ifelse(!is.na(g$edges$pvalue) & g$edges$pvalue < 0.05, "solid", "dashed")
g$edges$label        <- ifelse(is.na(g$edges$std.all), "", sprintf("%.2f", g$edges$std.all))
g$edges$label_family <- "serif"

# Route factor loadings vertically & force Growth -> SAR into top of SAR
node_y     <- setNames(g$nodes$y, g$nodes$name)
is_loading <- g$edges$op == "=~"
from_y     <- node_y[g$edges$from]
to_y       <- node_y[g$edges$to]

g$edges$connect_from[is_loading & to_y > from_y] <- "top"
g$edges$connect_to[is_loading & to_y > from_y]   <- "bottom"
g$edges$connect_from[is_loading & to_y < from_y] <- "bottom"
g$edges$connect_to[is_loading & to_y < from_y]   <- "top"

is_growth_to_sar <- g$edges$op == "~" &
  g$edges$from == "X06_DFA_IGF_mu_smoltyr" &
  g$edges$to   == "X16_DFA_sar_UCsprCh_smoltyr"

g$edges$connect_from[is_growth_to_sar] <- "bottom"
g$edges$connect_to[is_growth_to_sar]   <- "top"

# Render & Save
gg_M2_noPreyAK <- plot(g) +
  labs(title = "M2 Fish/Shark AK Predators + Latent NCC Prey & Predators") +
  theme(
    text       = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

print(gg_M2_noPreyAK)


ggsave(
  filename = file.path(plot_out_dir, "M2fish_preyNCC_predNCC.png"),
  plot     = gg_M2_noPreyAK,
  width    = 16,
  height   = 12,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

# Export companion parameter table
param_table_M2_noPreyAK <- pe_sub %>%
  filter(op %in% c("~", "=~")) %>%
  transmute(
    Path_Type = ifelse(op == "=~", "Factor Loading", "Structural Regression"),
    LHS       = lhs,
    Op        = op,
    RHS       = rhs,
    Std_Est   = round(std.all, 3),
    CI_Lower  = round(ci.lower, 3),
    CI_Upper  = round(ci.upper, 3),
    SE        = round(se, 3),
    P_Value   = round(pvalue, 4)
  )

write_csv(param_table_M2_noPreyAK, file.path(table_out_dir, "M2fish_preyNCC_predNCC_params.csv"))
