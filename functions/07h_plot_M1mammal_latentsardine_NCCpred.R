
# ==============================================================================
# Script: 02_analysis/07h_plot_M1mammal_latentsardine_NCCpred.R
# Step 1: get raw data from allData.rds
# Step 2: fit SEM
#dumped this one because after I expanded the EGOA DFA into a true latent variable, it was not significant, even w/o growth.
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
print(model_comparison_table %>% arrange(Model))          



# ==============================================================================
# 1. UPDATED 9x13 MASTER LAYOUT BUILDER (ACCOMMODATES 7 preyNCC INDICATORS)
# ==============================================================================
build_master_layout_13 <- function(
    preyNCC_node   = "preyNCC",
    predNCC_node   = "predNCC",
    preyAK_node    = "preyAK",
    predAK_node    = "predAK",
    growth_node    = "X06_DFA_IGF_mu_smoltyr",
    abundance_node = "X07_cpue_IntSprJunHW_smoltyr",
    sar_node       = "X16_DFA_sar_UCsprCh_smoltyr",
    preyNCC_ind    = character(0),
    predNCC_ind    = character(0),
    preyAK_ind     = character(0),
    predAK_ind     = character(0)
) {
  mat <- matrix("", nrow = 9, ncol = 13)
  
  # Locked 7-node structural core
  mat[3, 4]  <- preyNCC_node
  mat[3, 10] <- predNCC_node
  mat[4, 6]  <- growth_node
  mat[4, 8]  <- abundance_node
  mat[5, 7]  <- sar_node
  mat[7, 4]  <- preyAK_node
  mat[7, 10] <- predAK_node
  
  # Outer wings
  if (length(preyNCC_ind) > 0) mat[1, seq_along(preyNCC_ind)] <- preyNCC_ind           # Row 1, cols 1..7
  if (length(predNCC_ind) > 0) mat[1, 13 - seq_along(predNCC_ind) + 1] <- predNCC_ind # Row 1, cols 13..10
  if (length(preyAK_ind)  > 0) mat[9, seq_along(preyAK_ind) + 1]  <- preyAK_ind       # Row 9, cols 2..6
  if (length(predAK_ind)  > 0) mat[9, 13 - seq_along(predAK_ind) + 1]  <- predAK_ind  # Row 9, cols 13..9
  
  mat
}

# Grab exact surviving X05_latent_* column names automatically
x05_latent_cols <- grep("^X05_latent_", names(sem_data_merged), value = TRUE)

layout_M1_7prey <- build_master_layout_13(
  preyNCC_node = "preyNCC",
  predNCC_node = "predNCC",
  preyAK_node  = "",
  predAK_node  = "predAK",
  preyNCC_ind  = c(
    "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr",
    "X02_DFA_SeaNettle_smoltyr",
    "X04_marketsquid_GAM_smoltyr",
    x05_latent_cols
  ),
  predNCC_ind  = c(
    "X09_latent_HakeAge5Plus",
    "X09_latent_pacificmackerel_GAM",
    "X09_latent_jackmackerel_GAM",
    "X09_latent_hakeFisheryWA"
  ),
  predAK_ind   = c(
    "DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z",
    "DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z",
    "Killer_whales_NR_BC_adultyr__peff_Temp_z",
    "harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z",
    "ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z"
  )
)

# ==============================================================================
# 2. LABELS (ADD PRETTY NAMES FOR YOUR 4 X05_latent_* VARIABLES BELOW)
# ==============================================================================
# Print full names so you can verify the 4 X05_latent_* labels
print(x05_latent_cols)

labels_M1_7prey <- c(
  # Latent & Core Structural Nodes
  "preyNCC"                      = "NCC Prey",
  "predNCC"                      = "NCC Fish\nPredators",
  "predAK"                       = "Mammal\nPredators",
  "X06_DFA_IGF_mu_smoltyr"       = "Growth",
  "X07_cpue_IntSprJunHW_smoltyr" = "CPUE",
  "X16_DFA_sar_UCsprCh_smoltyr"  = "SAR",
  
  # preyNCC Indicators (update the 4 X05_latent_* keys/labels to match your exact names)
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr" = "Plankton\n(Jun)",
  "X02_DFA_SeaNettle_smoltyr"                 = "Sea\nNettle",
  "X04_marketsquid_GAM_smoltyr"               = "Market\nSquid",
  "X05_latent_northAnch_NCC"                  = "Northern\nAnchovy",
  "X05_latent_sardine_NCC"                    = "Sardine\n(NCC)",
  "X05_latent_herring_GAM"                    = "Herring\n(GAM)",
  "X05_latent_sardine_GAM"                    = "Sardine\n(GAM)",
  
  # predNCC Indicators
  "X09_latent_HakeAge5Plus"        = "Hake\nAge 5+",
  "X09_latent_pacificmackerel_GAM" = "Pacific\nMackerel",
  "X09_latent_jackmackerel_GAM"    = "Jack\nMackerel",
  "X09_latent_hakeFisheryWA"       = "WA Hake\nFishery",
  
  # predAK Indicators
  "DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z"   = "Steller SL\n(WS)",
  "DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z"    = "Sea Lions\n(EMB)",
  "Killer_whales_NR_BC_adultyr__peff_Temp_z"       = "NR Killer\nWhales",
  "harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z" = "Harb. Seal\n(WCVI)",
  "ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z"     = "SEAK SSL"
)

# ==============================================================================
# 3. BUILD & RENDER GRAPH (BASE-R SAFE, NO SELF-LOOPS, GROWTH -> SAR TOP ENTRY)
# ==============================================================================
g <- prepare_graph(
  model           = fit,
  layout          = layout_M1_7prey,
  intercepts      = FALSE,
  rect_width      = 1.6,   # Slightly narrower so 7 adjacent boxes in row 1 don't touch
  rect_height     = 0.9,
  ellipses_width  = 1.9,
  ellipses_height = 1.1
)

print(g$edges[g$edges$op == "~", c("from", "op", "to")])

g$nodes$label <- g$nodes$name
matches <- match(g$nodes$label, names(labels_M1_7prey))
g$nodes$label[!is.na(matches)] <- labels_M1_7prey[matches[!is.na(matches)]]

rsq     <- inspect(fit, "r2")
r2_keep <- c("X06_DFA_IGF_mu_smoltyr", "X07_cpue_IntSprJunHW_smoltyr", "X16_DFA_sar_UCsprCh_smoltyr")
rsq_sub <- round(unlist(rsq[names(rsq) %in% r2_keep]), 2)
r2_idx  <- match(names(rsq_sub), g$nodes$name)
g$nodes$label[r2_idx[!is.na(r2_idx)]] <- paste0(
  g$nodes$label[r2_idx[!is.na(r2_idx)]], "\nR² = ", rsq_sub[!is.na(r2_idx)]
)

g$nodes$label_family <- "serif"
g$nodes$fill         <- "grey96"
g$nodes$label_fill   <- "transparent"

g$edges <- g$edges[g$edges$from != g$edges$to, ]

pe <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
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

# Stagger the 7 preyNCC loading labels slightly so the fan of 7 numbers doesn't touch
# if (any(is_loading)) {
#   g$edges$label_location[is_loading] <- rep(c(0.45, 0.62), length.out = sum(is_loading))
# }

gg <- plot(g) +
  labs(title = "M1 Mammals + Fully Decomposed NCC Prey & Predators") +
  theme(
    text       = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

print(gg)
ggsave(
  filename = file.path(plot_out_dir, "M1mammals_preyNCC7_predNCC.png"),
  plot     = gg,
  width    = 18,
  height   = 12,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)