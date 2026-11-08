# ==============================================================================
# 22b_Plot_to_pdf_latent_models_from_DAG1B_AKPrey_Peff_2x2.r
# ==============================================================================

library(tidyverse)
library(lavaan)
library(tidySEM)
library(ggplot2)

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

# ==============================================================================
# 0. DIRECTORY & DATA SETUP
# ==============================================================================
proj_dir   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
cw_path    <- file.path(proj_dir, "metadata/crosswalk_names.csv")
rdata_path <- file.path(proj_dir, "Routput_for_paper/data/models_DAG1B_latent_stages.RData")
if (!file.exists(rdata_path)) {
  rdata_path <- file.path(proj_dir, "Routput_for_paper/models_DAG1B_latent_stages.RData")
}
pdf_out    <- file.path(proj_dir, "Routput_for_paper/figures/DAG1B_2Stage_Fits_by_Predator_2x2.pdf")

if (!exists("dag1b_candidate_models")) {
  if (file.exists(rdata_path)) {
    load(rdata_path)
    cat("Loaded 'dag1b_candidate_models' from:", rdata_path, "\n")
  } else {
    stop("Saved model RData file not found at: ", rdata_path)
  }
}

cw_df <- read.csv(cw_path, stringsAsFactors = FALSE)

# ==============================================================================
# 1. ORIGINAL 10-ROW MASTER COORDINATE LAYOUT
# ==============================================================================
layout_master <- get_layout(
  "sardine_GAM",         "",                  "",                             "",                                    "",
  "sardine_NCC",         "Latent_05_PreyNCC", "X06_DFA_IGF_mu_smoltyr",       "",                                    "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr",
  "northAnch_NCC",       "",                  "",                             "",                                    "",
  "herring_GAM",         "",                  "",                             "X16_DFA_sar_UCsprCh_smoltyr",         "",
  "",                    "",                  "",                             "",                                    "",
  "HakeAge5Plus",        "",                  "",                             "",                                    "P_eff",
  "hakeFisheryWA",       "Latent_09_PredNCC", "X07_cpue_IntSprJunHW_smoltyr", "",                                    "",
  "jackmackerel_GAM",    "",                  "",                             "",                                    "",
  "pacificmackerel_GAM", "",                  "",                             "",                                    "",
  "iphc20bigSkate",      "",                  "",                             "",                                    "",
  rows = 10
)

get_clean_label <- function(raw_var) {
  row_match <- cw_df %>% filter(DFAname == raw_var | dfa_cols == raw_var | LisaName == raw_var)
  clean_name <- if (nrow(row_match) > 0 && !is.na(row_match$LisaName[1]) && row_match$LisaName[1] != "") {
    row_match$LisaName[1]
  } else {
    raw_var
  }
  clean_name <- gsub("_smoltyr|_adultyr", "", clean_name)
  clean_name <- gsub("^X[0-9]{2}_", "", clean_name)
  clean_name
}

# ==============================================================================
# 2. ORIGINAL STANDARDIZED PLOT GENERATOR (PRESERVES OVALS)
# ==============================================================================
build_stage_plot <- function(fit_obj, stage_title, col_header = NULL, pred_label = "P_eff") {
  
  node_labels_map <- c(
    "X06_DFA_IGF_mu_smoltyr"               = "Salmon Growth\n[IGF-1]",
    "X07_cpue_IntSprJunHW_smoltyr"         = "Salmon Abund.\n[CPUE]",
    "X16_DFA_sar_UCsprCh_smoltyr"          = "Salmon SAR\n[UCsprCh]",
    "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"= "Prey (AK)\n[EGOA Lower]",
    "P_eff"                                = paste0("P_eff\n[", pred_label, "]"),
    "Latent_05_PreyNCC"                    = "Prey NCC\n(Latent)",
    "Latent_09_PredNCC"                    = "Pred NCC\n(Latent)",
    "northAnch_NCC"                        = "N. Anchovy",
    "sardine_NCC"                          = "Sardine NCC",
    "herring_GAM"                          = "Herring GAM",
    "sardine_GAM"                          = "Sardine GAM",
    "HakeAge5Plus"                         = "Hake Age 5+",
    "iphc20bigSkate"                       = "Big Skate",
    "hakeFisheryWA"                        = "Hake Fishery",
    "jackmackerel_GAM"                     = "Jack Mackerel",
    "pacificmackerel_GAM"                  = "Pac. Mackerel"
  )
  
  # Original graph parameters
  g <- prepare_graph(
    model          = fit_obj,
    layout         = layout_master,
    intercepts     = FALSE,
    rect_width     = 1.8,
    rect_height    = 0.65,
    ellipses_width = 1.8,
    ellipses_height= 0.65
  )
  
  # Assign clean labels
  g$nodes$label <- g$nodes$name
  matches <- match(g$nodes$label, names(node_labels_map))
  g$nodes$label[!is.na(matches)] <- node_labels_map[matches[!is.na(matches)]]
  
  # Merge parameter estimates for path coefficients
  pe <- parameterEstimates(fit_obj, standardized = TRUE, ci = TRUE)
  req_cols <- c("lhs", "op", "rhs", "std.all", "ci.lower", "ci.upper", "pvalue")
  for (col in req_cols) if (!col %in% names(pe)) pe[[col]] <- NA_real_
  pe_labeled <- pe[pe$op %in% c("~", "=~", "~~") & pe$lhs != pe$rhs, req_cols]
  
  g$edges <- merge(g$edges, pe_labeled, by.x = c("rhs", "op", "lhs"),
                   by.y = c("rhs", "op", "lhs"), all.x = FALSE)
  
  g$edges$linetype <- ifelse(!is.na(g$edges$pvalue) & g$edges$pvalue < 0.05, "solid", "dashed")
  g$edges$label    <- as.character(round(g$edges[["std.all"]], 2))
  
  # Append R-squared on endogenous performance nodes
  rsq <- tryCatch(inspect(fit_obj, "r2"), error = function(e) NULL)
  if (!is.null(rsq)) {
    rsq_df <- data.frame(name = names(rsq), r2 = round(unlist(rsq), 2), stringsAsFactors = FALSE)
    rsq_df <- rsq_df[rsq_df$name %in% c("X06_DFA_IGF_mu_smoltyr", "X07_cpue_IntSprJunHW_smoltyr", "X16_DFA_sar_UCsprCh_smoltyr"), ]
    g$nodes <- merge(g$nodes, rsq_df, by = "name", all.x = TRUE)
    g$nodes$label <- ifelse(is.na(g$nodes$r2), g$nodes$label, paste0(g$nodes$label, "\nR² = ", g$nodes$r2))
  }
  
  g$nodes$label_family <- g$edges$label_family <- "serif"
  g$nodes$fill         <- "grey96"
  g$nodes$label_fill   <- "transparent"
  
  full_title <- if (!is.null(col_header)) paste0(col_header, "\n", stage_title) else stage_title
  
  p <- plot(g) +
    labs(title = full_title) +
    theme(
      text        = element_text(family = "serif"),
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 11, lineheight = 1.1),
      plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
    )
  
  return(p)
}

# ==============================================================================
# 3. BUILD 2-PLOT COLUMN (STAGE 1 & STAGE 3 ONLY)
# ==============================================================================
render_predator_column_2stage <- function(model_entry) {
  meta     <- model_entry$metadata
  mod_id   <- meta$model_id
  pred_col <- meta$predator
  variant  <- meta$variant
  clean_p  <- get_clean_label(pred_col)
  
  col_header <- sprintf("Model #%d: %s\n[%s]", mod_id, clean_p, variant)
  
  # Stage 1: Full Original Measurement Model
  fit_s1 <- model_entry$Stage1_Full_Measurement$fit
  aic_s1 <- round(as.numeric(fitMeasures(fit_s1, "aic")), 2)
  t_s1   <- sprintf("Stage 1: Full Measurement (AIC = %.2f)", aic_s1)
  p1     <- build_stage_plot(fit_s1, stage_title = t_s1, col_header = col_header, pred_label = clean_p)
  
  # Stage 3: Final Pruned Structural Model (Pruned Latents + Pruned Regressions)
  fit_s3 <- model_entry$Stage3_Final_Pruned_Regression$fit
  aic_s3 <- round(as.numeric(fitMeasures(fit_s3, "aic")), 2)
  r2_val <- round(as.numeric(inspect(fit_s3, "r2")["X16_DFA_sar_UCsprCh_smoltyr"]), 2)
  t_s3   <- sprintf("Stage 3: Pruned Structural (AIC = %.2f | SAR R² = %.2f)", aic_s3, r2_val)
  p3     <- build_stage_plot(fit_s3, stage_title = t_s3, col_header = NULL, pred_label = clean_p)
  
  # Return stacked column of 2 plots
  return(p1 / p3)
}

# ==============================================================================
# 4. COMPILE 3-PAGE PDF (2 ROWS X 2 COLUMNS PER PAGE)
# ==============================================================================
fig_dir <- dirname(pdf_out)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

while (!is.null(dev.list())) dev.off()
if (file.exists(pdf_out)) tryCatch(file.remove(pdf_out), error = function(e) NULL)

# Standard Tabloid Landscape (17" wide x 13" tall) gives optimal vertical headroom for 2 rows
pdf(pdf_out, width = 17, height = 13, onefile = TRUE, family = "serif")

cat("Compiling 3-page comparative 2x2 PDF:\n", pdf_out, "\n\n")

# PAGE 1: Models 1 & 3
cat("Rendering Page 1 (Model 1 vs Model 3)...\n")
col_m1 <- render_predator_column_2stage(dag1b_candidate_models$Model_1)
col_m3 <- render_predator_column_2stage(dag1b_candidate_models$Model_3)
page1  <- col_m1 | col_m3
print(page1)

# PAGE 2: Models 4 & 5
cat("Rendering Page 2 (Model 4 vs Model 5)...\n")
col_m4 <- render_predator_column_2stage(dag1b_candidate_models$Model_4)
col_m5 <- render_predator_column_2stage(dag1b_candidate_models$Model_5)
page2  <- col_m4 | col_m5
print(page2)

# PAGE 3: Model 7 + Blank Spacer Column
cat("Rendering Page 3 (Model 7 + Spacer)...\n")
col_m7 <- render_predator_column_2stage(dag1b_candidate_models$Model_7)
page3  <- col_m7 | plot_spacer()
print(page3)

dev.off()
cat("\nProcess complete! PDF successfully written to:\n", pdf_out, "\n")