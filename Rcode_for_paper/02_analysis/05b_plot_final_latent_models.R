# ==============================================================================
# Script: 02_analysis/05b_plot_final_latent_models.R
# Purpose: Generate Full-Page Path Diagrams with Custom Box Sizing & White Background
# Outputs:
#   - Outputs/figures/obj4_fig1_final_pruned_models_multipage.pdf
#   - Outputs/figures/obj4_fig1A_sablefish_model.png
#   - Outputs/figures/obj4_fig1B_sleepershark_model.png
#   - Outputs/figures/obj4_fig1C_coastal_baseline_model.png
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(tidySEM)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & LOAD MODELS
# ------------------------------------------------------------------------------
proj_dir    <- file.path(getwd(), "Rcode_for_paper")
master_out  <- file.path(proj_dir, "Outputs")
data_dir    <- file.path(master_out, "data")
fig_out_dir <- file.path(master_out, "figures")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

rdata_path <- file.path(data_dir, "obj4_latent_model_objects.RData")
if (!file.exists(rdata_path)) stop("Model file not found: ", rdata_path)
load(rdata_path)

# ------------------------------------------------------------------------------
# STEP 1: EXPANDED 9-COLUMN LAYOUT (ELIMINATES COLLISION)
# ------------------------------------------------------------------------------
# Adding spacer columns gives wide boxes room to breathe without colliding with arrows
layout_master <- get_layout(
  "sardine_GAM",         "", "",                  "", "",                             "", "",                            "", "",
  "sardine_NCC",         "", "Latent_05_PreyNCC", "", "X06_DFA_IGF_mu_smoltyr",       "", "",                            "", "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr",
  "northAnch_NCC",       "", "",                  "", "",                             "", "",                            "", "",
  "herring_GAM",         "", "",                  "", "",                             "", "X16_DFA_sar_UCsprCh_smoltyr", "", "",
  "",                    "", "",                  "", "",                             "", "",                            "", "",
  "HakeAge5Plus",        "", "",                  "", "",                             "", "",                            "", "P_eff",
  "hakeFisheryWA",       "", "Latent_09_PredNCC", "", "X07_cpue_IntSprJunHW_smoltyr", "", "",                            "", "",
  "jackmackerel_GAM",    "", "",                  "", "",                             "", "",                            "", "",
  "pacificmackerel_GAM", "", "",                  "", "",                             "", "",                            "", "",
  "iphc20bigSkate",      "", "",                  "", "",                             "", "",                            "", "",
  rows = 10
)

# ------------------------------------------------------------------------------
# STEP 2: STANDARDIZED FULL-PAGE PLOT ENGINE (LARGE BOX GEOMETRY)
# ------------------------------------------------------------------------------
build_clean_diagram <- function(model_entry, custom_title = NULL, custom_subtitle = NULL) {
  meta     <- model_entry$metadata
  fit_obj  <- model_entry$Stage3_Final_Pruned_Regression$fit
  aic_val  <- round(as.numeric(fitMeasures(fit_obj, "aic")), 2)
  r2_val   <- round(as.numeric(inspect(fit_obj, "r2")["X16_DFA_sar_UCsprCh_smoltyr"]), 2)
  
  pred_col <- meta$predator
  clean_p  <- case_when(
    grepl("sablefish", pred_col)   ~ "Sablefish (Age 3+)",
    grepl("sleeperShark", pred_col)~ "Sleeper Shark",
    grepl("ssl_seak_pup", pred_col)~ "Steller Sea Lion (Pup)",
    grepl("Killer_whale", pred_col)~ "Killer Whale (NR)",
    grepl("pollock", pred_col)     ~ "Pollock (SSB)",
    TRUE ~ pred_col
  )
  
  plot_title <- if (!is.null(custom_title)) {
    custom_title
  } else {
    sprintf("Structural Model: %s [%s]", clean_p, meta$variant)
  }
  
  plot_sub <- if (!is.null(custom_subtitle)) {
    custom_subtitle
  } else {
    sprintf("AIC = %.2f  |  SAR R² = %.2f  |  Retained Pathways (p < 0.05)", aic_val, r2_val)
  }
  
  node_labels_map <- c(
    "X06_DFA_IGF_mu_smoltyr"                = "Salmon Growth\n[IGF-1]",
    "X07_cpue_IntSprJunHW_smoltyr"          = "Salmon CPUE\n[June Abund]",
    "X16_DFA_sar_UCsprCh_smoltyr"           = "Salmon SAR\n[Adult Return]",
    "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr" = "Prey (AK)\n[EGOA Lower]",
    "P_eff"                                 = paste0("P_eff\n[", clean_p, "]"),
    "Latent_05_PreyNCC"                     = "Prey NCC\n(Latent)",
    "Latent_09_PredNCC"                     = "Pred NCC\n(Latent)",
    "northAnch_NCC"                         = "N. Anchovy",
    "sardine_NCC"                           = "Sardine NCC",
    "herring_GAM"                           = "Herring",
    "sardine_GAM"                           = "Sardine GAM",
    "HakeAge5Plus"                          = "Hake Age 5+",
    "hakeFisheryWA"                         = "Hake Catch",
    "jackmackerel_GAM"                      = "Jack Mackerel",
    "pacificmackerel_GAM"                   = "Pac. Mackerel"
  )
  
  # Large box parameters passed directly to prepare_graph
  g <- prepare_graph(
    model           = fit_obj,
    layout          = layout_master,
    intercepts      = FALSE,
    rect_width      = 1.25,
    rect_height     = 0.52,
    ellipses_width  = 1.25,
    ellipses_height = 0.52
  )
  
  # Assign clean labels
  g$nodes$label <- g$nodes$name
  matches <- match(g$nodes$label, names(node_labels_map))
  g$nodes$label[!is.na(matches)] <- node_labels_map[matches[!is.na(matches)]]
  
  # Append R-squared on endogenous performance nodes
  rsq <- tryCatch(inspect(fit_obj, "r2"), error = function(e) NULL)
  if (!is.null(rsq)) {
    rsq_df <- data.frame(name = names(rsq), r2 = round(unlist(rsq), 2), stringsAsFactors = FALSE)
    rsq_df <- rsq_df[rsq_df$name %in% c("X06_DFA_IGF_mu_smoltyr", "X07_cpue_IntSprJunHW_smoltyr", "X16_DFA_sar_UCsprCh_smoltyr"), ]
    g$nodes <- merge(g$nodes, rsq_df, by = "name", all.x = TRUE)
    g$nodes$label <- ifelse(is.na(g$nodes$r2), g$nodes$label, paste0(g$nodes$label, "\nR² = ", g$nodes$r2))
  }
  
  # Edge styling and standardized parameter loadings
  pe <- parameterEstimates(fit_obj, standardized = TRUE)
  pe_labeled <- pe[pe$op %in% c("~", "=~", "~~") & pe$lhs != pe$rhs, c("lhs", "op", "rhs", "std.all", "pvalue")]
  
  g$edges <- merge(g$edges, pe_labeled, by.x = c("rhs", "op", "lhs"), by.y = c("rhs", "op", "lhs"), all.x = FALSE)
  g$edges$linetype <- ifelse(!is.na(g$edges$pvalue) & g$edges$pvalue < 0.05, "solid", "dashed")
  g$edges$label    <- as.character(round(g$edges[["std.all"]], 2))
  
  # Node styling: solid white fill with thin black outline
  g$nodes$label_family <- g$edges$label_family <- "sans"
  g$nodes$fill         <- "white"
  g$nodes$color        <- "black"
  g$nodes$label_fill   <- "transparent"
  
  p <- plot(g, text_size = 3.2) +
    labs(title = plot_title, subtitle = plot_sub) +
    theme_void() +
    theme(
      text               = element_text(family = "sans"),
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(b = 4)),
      plot.subtitle      = element_text(hjust = 0.5, size = 10.5, color = "grey30", margin = margin(b = 16)),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.margin        = margin(t = 20, r = 20, b = 20, l = 20)
    )
  
  return(p)
}

# ------------------------------------------------------------------------------
# STEP 3: RENDER THE 3 KEY MODELS
# ------------------------------------------------------------------------------
plot_A <- build_clean_diagram(
  dag1b_candidate_models$Model_2,
  custom_title    = "Figure 1A. Final Structural SEM: Subarctic Groundfish Piscivore (Sablefish Age 3+)",
  custom_subtitle = "Peff modulated by Summer SST and Herring/Pollock alternative forage | AIC = -23.22, SAR R² = 0.73"
)

plot_B <- build_clean_diagram(
  dag1b_candidate_models$Model_4,
  custom_title    = "Figure 1B. Final Structural SEM: Subarctic Apex Shark (Sleeper Shark)",
  custom_subtitle = "Peff modulated by Summer SST and Capelin forage buffering | Retained Top-Down Effect"
)

plot_C <- build_clean_diagram(
  dag1b_candidate_models$Model_1,
  custom_title    = "Figure 1C. Final Pruned Baseline Structural Architecture (Coastal-Only Survival)",
  custom_subtitle = "Pruned structural architecture resulting when Steller Sea Lion, Killer Whale, and Pollock Peff drop out (p >= 0.05)"
)

# ------------------------------------------------------------------------------
# STEP 4: EXPORT MULTI-PAGE PDF & INDIVIDUAL WHITE-BACKED PNGs
# ------------------------------------------------------------------------------
pdf_multipage <- file.path(fig_out_dir, "obj4_fig1_final_pruned_models_multipage.pdf")

cairo_pdf(pdf_multipage, width = 11, height = 8.5, onefile = TRUE)
print(plot_A)
print(plot_B)
print(plot_C)
dev.off()

ggsave(file.path(fig_out_dir, "obj4_fig1A_sablefish_model.png"),        plot = plot_A, width = 11, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(fig_out_dir, "obj4_fig1B_sleepershark_model.png"),      plot = plot_B, width = 11, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(fig_out_dir, "obj4_fig1C_coastal_baseline_model.png"), plot = plot_C, width = 11, height = 8.5, dpi = 300, bg = "white")

cat(sprintf("\nSuccessfully exported figures:\n  -> %s\n  -> %s\n  -> %s\n  -> %s\n",
            pdf_multipage,
            file.path(fig_out_dir, "obj4_fig1A_sablefish_model.png"),
            file.path(fig_out_dir, "obj4_fig1B_sleepershark_model.png"),
            file.path(fig_out_dir, "obj4_fig1C_coastal_baseline_model.png")))