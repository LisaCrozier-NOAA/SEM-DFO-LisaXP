
# ==============================================================================
# Script: 02_analysis/07g_plot_M1mammal_NCCprey_NCCpred.R
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
print(model_comparison_table %>% arrange(Model))          


#   Model                                     npar    df Chisq p_Chisq   AIC   BIC   CFI  SRMR Growth_R2 CPUE_R2 SAR_R2 Notes                                                                        
# 1 M1mammal_noPreyAK                           25    35 126.        0  398.  428. 0.663 0.365     0.451   0.749  0.823 5-mammal predAK (var floor 0.01), no preyAK, fixed.x=TRUE                    
# 2 M1mammals_preyNCC                           38    64 178.        0  506.  550. 0.758 0.159     0.898   0.749  0.701 5-mammal predAK (var floor 0.01), latent NCC, growth, no preyAK, fixed.x=TRUE
# 3 M1mammals_preyNCC_predNCC                   52   100 251.        0  781.  842. 0.713 0.135     0.895   0.674  0.642 Mammals,latent NCCpred                                                       

#for the mammal models, adding preyNCC helps CFI, SRMR and Growth_R2, but it hurts SAR_R2; adding predNCC hurts R2 (cpue and sar)
#but in fact, predNCC is the assumed process under all of the models, so it is disingenious to choose a higher R2 bec you have a DFA instead of the raw variables

#2. One Small Tweaks Needed in Lisa_sem_graph.v3
#In your current Lisa_sem_graph.v3, tidySEM::prepare_graph() can draw curved covariance arrows (~~) among exogenous indicators when fixed.x = TRUE—which creates a tangled web of dashed lines across the top and bottom indicator rows. Adding variance_diameter = 0 (or filtering g$edges to keep only ~ and =~) keeps the diagram clean while still preserving your linetype styling:
# Inside Lisa_sem_graph.v3, right after prepare_graph(...):
#g$edges <- g$edges[g$edges$op %in% c("~", "=~"), ]

# ==============================================================================
# 1. LOCKED 9x11 MASTER LAYOUT BUILDER -------
# ==============================================================================
build_master_layout <- function(
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
  mat <- matrix("", nrow = 9, ncol = 11)
  
  # Locked 7-node structural core
  mat[3, 3] <- preyNCC_node
  mat[3, 9] <- predNCC_node
  mat[4, 5] <- growth_node
  mat[4, 7] <- abundance_node
  mat[5, 6] <- sar_node
  mat[7, 3] <- preyAK_node
  mat[7, 9] <- predAK_node
  
  # Outer wings for expanded latent indicators (up to 5 per corner)
  if (length(preyNCC_ind) > 0) mat[1, seq_along(preyNCC_ind) + 1] <- preyNCC_ind # Row 1, cols 2..6
  if (length(predNCC_ind) > 0) mat[1, 11 - seq_along(predNCC_ind) + 1] <- predNCC_ind # Row 1, cols 11..7
  if (length(preyAK_ind)  > 0) mat[9, seq_along(preyAK_ind) + 1]  <- preyAK_ind  # Row 9, cols 2..6
  if (length(predAK_ind)  > 0) mat[9, 11 - seq_along(predAK_ind) + 1]  <- predAK_ind  # Row 9, cols 11..7
  
  mat
}

# ==============================================================================
# 2. M1mammals_preyNCC_predNCC ---------
# ==============================================================================
layout_M1_full <- build_master_layout(
  preyNCC_node = "preyNCC",
  predNCC_node = "predNCC",
  preyAK_node  = "", # Absent in M1mammals_preyNCC_predNCC (space preserved)
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
    "DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z",
    "DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z",
    "Killer_whales_NR_BC_adultyr__peff_Temp_z",
    "harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z",
    "ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z"
  )
)

labels_M1_full <- c(
  # Latent & Core Structural Nodes
  "preyNCC"                      = "NCC Prey",
  "predNCC"                      = "NCC Fish\nPredators",
  "predAK"                       = "Mammal \nPredators",
  "X06_DFA_IGF_mu_smoltyr"       = "Growth",
  "X07_cpue_IntSprJunHW_smoltyr" = "CPUE",
  "X16_DFA_sar_UCsprCh_smoltyr"  = "SAR",
  
  # preyNCC Indicators
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr" = "Plankton\n(Jun)",
  "X02_DFA_SeaNettle_smoltyr"                 = "Sea\nNettle",
  "X04_marketsquid_GAM_smoltyr"               = "Market\nSquid",
  "X05_DFA_abundSardine_smoltyr"              = "Sardine\nAbund",
  
  # predNCC Indicators
  "X09_latent_HakeAge5Plus"        = "Hake\nAge 5+",
  "X09_latent_pacificmackerel_GAM" = "Pacific\nMackerel",
  "X09_latent_jackmackerel_GAM"    = "Jack\nMackerel",
  "X09_latent_hakeFisheryWA"       = "WA Hake\nFishery",
  
  # predAK Indicators
  "DFA_Steller_s_l_WS_adultyr__peff_TempPrey1_z"   = "Steller SL\n(WS)",
  "DFA_AllSeaLionsEMB_adultyr__peff_AltPrey1_z"    = "Sea Lions \n(EMB)",
  "Killer_whales_NR_BC_adultyr__peff_Temp_z"       = "NR Killer\nWhales",
  "harbour_seal_DFO_WCVI_smoltyr__peff_AltPrey1_z" = "Harb. Seal\n(WCVI)",
  "ssl_seak_pup_pred_smoltyr__peff_AltPrey2_z"     = "SEAK SSL"
)

# #DUMP------
# plot_M1_full <- Lisa_sem_graph.v3(
#   fit               = fit_M1mammals_preyNCC_predNCC,
#   title             = "M1 Mammals + Latent NCC Prey & Predators",
#   layout            = layout_M1_full,
#   node_labels       = labels_M1_full,
#   rect_width        = 1.8,   # Default is 1.2 (widens observed variable squares)
#   rect_height       = 0.9,   # Default is 0.8
#   ellipses_width    = 1.7,   # Default is 1.0 (widens latent variable ovals)
#   ellipses_height   = 1.1,   # Default is 1.0
#   show_coefficients = TRUE,
#   show_r2           = TRUE,
#   r2_nodes          = c("X06_DFA_IGF_mu_smoltyr", "X07_cpue_IntSprJunHW_smoltyr", "X16_DFA_sar_UCsprCh_smoltyr"),
#   width             = 16,
#   height            = 10,
#   save              = TRUE,
#   savename          = file.path(plot_out_dir, "M1mammals_preyNCC_predNCC.png"),
#   display           = TRUE
# )
# #OR-----------
# # 1. Build the graph object without displaying/saving yet
# g <- Lisa_sem_graph.v3(
#   fit               = fit_M1mammals_preyNCC_predNCC,
#   title             = "M1 Mammals + Latent NCC Prey & Predators",
#   layout            = layout_M1_full,
#   node_labels       = labels_M1_full,
#   rect_width        = 1.8,
#   rect_height       = 0.9,
#   ellipses_width    = 1.7,
#   ellipses_height   = 1.1,
#   show_coefficients = TRUE,
#   show_r2           = TRUE,
#   r2_nodes          = c("X06_DFA_IGF_mu_smoltyr", "X07_cpue_IntSprJunHW_smoltyr", "X16_DFA_sar_UCsprCh_smoltyr"),
#   save              = FALSE,
#   display           = FALSE
# )
# 
# # 2. Customize arrow routing & stagger label positions ONLY for this plot
# node_y     <- setNames(g$nodes$y, g$nodes$name)
# is_loading <- g$edges$op == "=~"
# from_y     <- node_y[g$edges$from]
# to_y       <- node_y[g$edges$to]
# 
# # Top-row indicators: leave TOP of oval -> enter BOTTOM of box
# g$edges$connect_from[is_loading & to_y > from_y] <- "top"
# g$edges$connect_to[is_loading & to_y > from_y]   <- "bottom"
# 
# # Bottom-row indicators: leave BOTTOM of oval -> enter TOP of box
# g$edges$connect_from[is_loading & to_y < from_y] <- "bottom"
# g$edges$connect_to[is_loading & to_y < from_y]   <- "top"
# 
# # Stagger loading label positions along the arrows so they don't overlap
# if (any(is_loading)) {
#   g$edges$label_location[is_loading] <- rep(c(0.35, 0.75), length.out = sum(is_loading))
# }
# 
# # 3. Render and save
# gg <- plot(g) +
#   labs(title = "M1 Mammals + Latent NCC Prey & Predators") +
#   theme(
#     text       = element_text(family = "serif"),
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
#   )
# 
# print(gg)
# ggsave(
#   filename = file.path(plot_out_dir, "M1mammals_preyNCC_predNCC.png"),
#   plot     = gg,
#   width    = 16,
#   height   = 10,
#   units    = "in",
#   dpi      = 300,
#   bg       = "white"
# )
# 
# #OR==========
# # 1. Get the base graph directly from prepare_graph so NO regression lines are dropped
# g <- prepare_graph(
#   model           = fit_M1mammals_preyNCC_predNCC,
#   layout          = layout_M1_full,
#   intercepts      = FALSE,
#   rect_width      = 1.8,
#   rect_height     = 0.9,
#   ellipses_width  = 1.7,
#   ellipses_height = 1.1
# )
# 
# # 2. Apply your pretty node labels + R² on endogenous core nodes
# g$nodes$label <- g$nodes$name
# matches <- match(g$nodes$label, names(labels_M1_full))
# g$nodes$label[!is.na(matches)] <- labels_M1_full[matches[!is.na(matches)]]
# 
# rsq <- inspect(fit_M1mammals_preyNCC_predNCC, "r2")
# r2_keep <- c("X06_DFA_IGF_mu_smoltyr", "X07_cpue_IntSprJunHW_smoltyr", "X16_DFA_sar_UCsprCh_smoltyr")
# rsq_df <- data.frame(name = names(rsq), r2 = round(unlist(rsq), 2), stringsAsFactors = FALSE) %>%
#   filter(name %in% r2_keep)
# 
# g$nodes <- left_join(g$nodes, rsq_df, by = "name") %>%
#   mutate(label = ifelse(is.na(r2), label, paste0(label, "\nR² = ", r2)))
# 
# g$nodes$label_family <- "serif"
# g$nodes$fill         <- "grey96"
# g$nodes$label_fill   <- "transparent"
# 
# # 3. Attach standardized coefficients ONLY (no CIs) and significance linetypes
# pe <- parameterEstimates(fit_M1mammals_preyNCC_predNCC, standardized = TRUE, ci = TRUE) %>%
#   filter(op %in% c("~", "=~", "~~"), lhs != rhs) %>%
#   mutate(
#     # In lavaan: for "~", from = rhs, to = lhs; for "=~" and "~~", from = lhs, to = rhs
#     from_match = ifelse(op == "~", rhs, lhs),
#     to_match   = ifelse(op == "~", lhs, rhs)
#   )
# 
# # 3b. Remove empty self-variance circles (from == to), keep inter-variable covariances (from != to)
# g$edges <- g$edges %>%
#   filter(from != to) %>%
#   left_join(
#     pe %>% select(from_match, to_match, op, std.all, pvalue),
#     by = c("from" = "from_match", "to" = "to_match", "op" = "op")
#   ) %>%
#   mutate(
#     linetype     = ifelse(!is.na(pvalue) & pvalue < 0.05, "solid", "dashed"),
#     label        = ifelse(is.na(std.all), "", sprintf("%.2f", std.all)),
#     label_family = "serif"
#   )
# 
# # 4. Route factor loadings (op == "=~") cleanly to top/bottom & stagger labels
# node_y     <- setNames(g$nodes$y, g$nodes$name)
# is_loading <- g$edges$op == "=~"
# from_y     <- node_y[g$edges$from]
# to_y       <- node_y[g$edges$to]
# 
# g$edges$connect_from[is_loading & to_y > from_y] <- "top"
# g$edges$connect_to[is_loading & to_y > from_y]   <- "bottom"
# g$edges$connect_from[is_loading & to_y < from_y] <- "bottom"
# g$edges$connect_to[is_loading & to_y < from_y]   <- "top"
# 
# # if (any(is_loading)) {
# #   g$edges$label_location[is_loading] <- rep(c(0.40, 0.55, 0.70), length.out = sum(is_loading))
# # }
# 
# 
# # Force Growth -> SAR to enter the TOP of SAR (consistent with CPUE -> SAR)
# is_growth_to_sar <- g$edges$op == "~" &
#   g$edges$from == "X06_DFA_IGF_mu_smoltyr" &
#   g$edges$to   == "X16_DFA_sar_UCsprCh_smoltyr"
# 
# g$edges$connect_from[is_growth_to_sar] <- "bottom"
# g$edges$connect_to[is_growth_to_sar]   <- "top"
# 
# if (any(is_loading)) {
#   g$edges$label_location[is_loading] <- rep(c(0.40, 0.55, 0.70), length.out = sum(is_loading))
# }
# 
# # 5. Plot & Save
# gg <- plot(g) +
#   labs(title = "M1 Mammals + Latent NCC Prey & Predators") +
#   theme(
#     text       = element_text(family = "serif"),
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
#   )
# 
# print(gg)
# ggsave(
#   filename = file.path(plot_out_dir, "M1mammals_preyNCC_predNCC.png"),
#   plot     = gg,
#   width    = 16,
#   height   = 12,
#   units    = "in",
#   dpi      = 300,
#   bg       = "white"
# )



#FINAL---------
# 1. Get the base graph directly from prepare_graph
g <- prepare_graph(
  model           = fit_M1mammals_preyNCC_predNCC,
  layout          = layout_M1_full,
  intercepts      = FALSE,
  rect_width      = 1.8,
  rect_height     = 0.9,
  ellipses_width  = 1.7,
  ellipses_height = 1.1
)

# 2. Apply pretty node labels + R² using base R (preserves tidySEM node structure)
g$nodes$label <- g$nodes$name
matches <- match(g$nodes$label, names(labels_M1_full))
g$nodes$label[!is.na(matches)] <- labels_M1_full[matches[!is.na(matches)]]

rsq <- inspect(fit_M1mammals_preyNCC_predNCC, "r2")
r2_keep <- c("X06_DFA_IGF_mu_smoltyr", "X07_cpue_IntSprJunHW_smoltyr", "X16_DFA_sar_UCsprCh_smoltyr")
rsq_sub <- round(unlist(rsq[names(rsq) %in% r2_keep]), 2)
r2_idx  <- match(names(rsq_sub), g$nodes$name)
g$nodes$label[r2_idx[!is.na(r2_idx)]] <- paste0(
  g$nodes$label[r2_idx[!is.na(r2_idx)]], "\nR² = ", rsq_sub[!is.na(r2_idx)]
)

g$nodes$label_family <- "serif"
g$nodes$fill         <- "grey96"
g$nodes$label_fill   <- "transparent"

# 3. Remove empty self-variance circles (from == to) using base R subsetting
g$edges <- g$edges[g$edges$from != g$edges$to, ]

# Match standardized estimates & p-values without reordering g$edges
pe <- parameterEstimates(fit_M1mammals_preyNCC_predNCC, standardized = TRUE)
pe <- pe[pe$op %in% c("~", "=~", "~~") & pe$lhs != pe$rhs, ]
pe$from_match <- ifelse(pe$op == "~", pe$rhs, pe$lhs)
pe$to_match   <- ifelse(pe$op == "~", pe$lhs, pe$rhs)

# Match both directions for symmetric covariances ("~~")
edge_key <- paste(g$edges$from, g$edges$op, g$edges$to)
pe_key1  <- paste(pe$from_match, pe$op, pe$to_match)
pe_key2  <- paste(pe$to_match, pe$op, pe$from_match)

m_idx <- match(edge_key, pe_key1)
m_idx[is.na(m_idx)] <- match(edge_key[is.na(m_idx)], pe_key2)

g$edges$std.all      <- pe$std.all[m_idx]
g$edges$pvalue       <- pe$pvalue[m_idx]
g$edges$linetype     <- ifelse(!is.na(g$edges$pvalue) & g$edges$pvalue < 0.05, "solid", "dashed")
g$edges$label        <- ifelse(is.na(g$edges$std.all), "", sprintf("%.2f", g$edges$std.all))
g$edges$label_family <- "serif"

# 4. Route factor loadings (op == "=~") cleanly to top/bottom
node_y     <- setNames(g$nodes$y, g$nodes$name)
is_loading <- g$edges$op == "=~"
from_y     <- node_y[g$edges$from]
to_y       <- node_y[g$edges$to]

g$edges$connect_from[is_loading & to_y > from_y] <- "top"
g$edges$connect_to[is_loading & to_y > from_y]   <- "bottom"
g$edges$connect_from[is_loading & to_y < from_y] <- "bottom"
g$edges$connect_to[is_loading & to_y < from_y]   <- "top"

# Force Growth -> SAR to enter the TOP of SAR (consistent with CPUE -> SAR)
is_growth_to_sar <- g$edges$op == "~" &
  g$edges$from == "X06_DFA_IGF_mu_smoltyr" &
  g$edges$to   == "X16_DFA_sar_UCsprCh_smoltyr"

g$edges$connect_from[is_growth_to_sar] <- "bottom"
g$edges$connect_to[is_growth_to_sar]   <- "top"

# 5. Plot & Save
gg <- plot(g) +
  labs(title = "M1 Mammals + Latent NCC Prey & Predators") +
  theme(
    text       = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

print(gg)
ggsave(
  filename = file.path(plot_out_dir, "M1mammals_preyNCC_predNCC.png"),
  plot     = gg,
  width    = 16,
  height   = 12,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)


#EXPORT PARAMETER CI (SINCE WE REMOVED THEM FROM GRAPH)------
# 6. Export companion parameter table (with 95% CIs)
param_table_M1 <- pe %>%
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

write_csv(param_table_M1, file.path(table_out_dir, "M1mammals_preyNCC_predNCC_params.csv"))
