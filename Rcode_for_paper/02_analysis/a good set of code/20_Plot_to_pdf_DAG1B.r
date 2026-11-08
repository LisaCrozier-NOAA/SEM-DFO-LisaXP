#20_Plot_to_pdf_DAG1B.r
library(tidyverse)
library(lubridate)
library(lavaan)
library(tidySEM)
library(ggplot2)

# ==============================================================================
# 0. DIRECTORY & DATA SETUP
# ==============================================================================
clus_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DFA/clusDataDFA.rds"
cw_path   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/crosswalk_names.csv"
dir_1b    <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1B_long"
pdf_out   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/Routput_for_paper/figures/DAG1B_Long_Acceptable_PreyAK_Signs.pdf"

clusData <- readRDS(clus_path)
cw       <- read.csv(cw_path, stringsAsFactors = FALSE)
res_1b   <- read.csv(file.path(dir_1b, "SEMresultsByClus.csv"), stringsAsFactors = FALSE)
pe_1b    <- read.csv(file.path(dir_1b, "parameterEstimates.csv"), stringsAsFactors = FALSE)

# 24-year long dataset (1998-2021)
clusData_sub <- clusData %>%
  filter(minYear == 1998) %>%
  mutate(year = year(ymd(date))) %>%
  filter(year >= 1998, year <= 2021)

# ==============================================================================
# 1. IDENTIFY BOTH POSITIVE & NEGATIVE PREY(AK) MODELS
# ==============================================================================
min_aic_1b <- min(res_1b$AIC, na.rm = TRUE)

res_1b_fit <- res_1b %>%
  mutate(deltaAIC = round(AIC - min_aic_1b, 2)) %>%
  filter(p_Chi2 >= 0.05, CFI >= 0.90, SRMR <= 0.08)

# Extract slope for PreyAK -> logSAR using 'est' and 'pvalue'
preyak_slopes <- pe_1b %>%
  filter(op == "~" & lhs %in% c("logSAR", "SAR") & rhs == "PreyAK") %>%
  select(modNum, est, pvalue) %>%
  rename(beta_preyak = est, p_preyak = pvalue)

res_1b_annotated <- res_1b_fit %>%
  inner_join(preyak_slopes, by = "modNum") %>%
  arrange(AIC)

# Group A: Acceptable models with Delta AIC <= 4 (negative PreyAK)
models_neg_preyak <- res_1b_annotated %>%
  filter(deltaAIC <= 4.0 & beta_preyak < 0) %>%
  mutate(group_label = "DeltaAIC <= 4 (Negative PreyAK)")

# Group B: Top fitting models with POSITIVE PreyAK
models_pos_preyak <- res_1b_annotated %>%
  filter(beta_preyak > 0) %>%
  slice(1:4) %>%
  mutate(group_label = "Top Models (Positive PreyAK)")

cat("=== ACCEPTABLE MODELS WITH NEGATIVE PREY(AK) ===\n")
print(models_neg_preyak %>% select(modNum, AIC, deltaAIC, PreyAKindNames, beta_preyak, p_preyak))

cat("\n=== TOP FITTING MODELS WITH POSITIVE PREY(AK) ===\n")
print(models_pos_preyak %>% select(modNum, AIC, deltaAIC, PreyAKindNames, beta_preyak, p_preyak))

combined_targets <- bind_rows(models_neg_preyak, models_pos_preyak)

# ==============================================================================
# 2. UTILITY & GRAPHING FUNCTIONS (RETURNING GGPLOT OBJECT)
# ==============================================================================
pull_std_node <- function(cdat, var_node, ind_name) {
  df_node <- cdat %>%
    filter(SEMlatent == var_node, shortName == ind_name) %>%
    select(date, finalVal) %>%
    arrange(date)
  df_node[[var_node]] <- as.vector(scale(df_node$finalVal))
  df_node %>% select(date, all_of(var_node))
}

pull_std_sar <- function(cdat) {
  cdat %>%
    filter(SEMlatent == "SAR") %>%
    select(date, finalVal) %>%
    arrange(date) %>%
    mutate(logSAR = as.vector(scale(finalVal))) %>%
    select(date, logSAR)
}

build_explicit_label <- function(node_title, ind_key, cw_df) {
  row_match <- cw_df %>%
    filter(DFAname == ind_key | dfa_cols == ind_key)
  lisa_val <- if (nrow(row_match) > 0 && !is.na(row_match$LisaName[1]) && row_match$LisaName[1] != "") {
    row_match$LisaName[1]
  } else {
    ind_key
  }
  paste0(node_title, "\n[", lisa_val, "]")
}

Lisa_sem_graph.v3 <- function(fit,
                              title = NULL,
                              layout = NULL,
                              node_labels = NULL, 
                              show_coefficients = TRUE,
                              p_threshold = 0.05,
                              sig_linetype = "solid",
                              nonsig_linetype = "dashed",
                              show_r2 = TRUE,
                              r2_nodes = NULL, 
                              font_family = "serif",
                              fill_c = "grey96",
                              save = FALSE,
                              savename = "sem_plot.png",
                              display = TRUE,
                              width = 14,
                              height = 8,
                              dpi = 300,
                              ...) {
  
  g <- prepare_graph(model = fit, layout = layout, intercepts = FALSE, ...)
  g$nodes$label <- g$nodes$name 
  
  if (!is.null(node_labels)) {
    matches <- match(g$nodes$label, names(node_labels))
    g$nodes$label[!is.na(matches)] <- node_labels[matches[!is.na(matches)]]
  }
  
  result <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
  required_cols <- c("lhs", "op", "rhs", "std.all", "ci.lower", "ci.upper", "pvalue")
  for (col in required_cols) {
    if (!col %in% names(result)) result[[col]] <- NA_real_
  }
  
  result_labeled <- result[result$op %in% c("~", "=~", "~~") & result$lhs != result$rhs, required_cols]
  
  g$edges <- merge(g$edges, result_labeled, by.x = c("rhs", "op", "lhs"),
                   by.y = c("rhs", "op", "lhs"), all.x = FALSE)
  
  g$edges$linetype <- ifelse(!is.na(g$edges$pvalue) & g$edges$pvalue < p_threshold, 
                             sig_linetype, 
                             nonsig_linetype)
  
  if (show_coefficients) {
    if (all(!is.na(g$edges$ci.lower)) && all(!is.na(g$edges$ci.upper))) {
      g$edges$label <- paste0(round(g$edges[["std.all"]], 2), "\n[",
                              round(g$edges[["ci.lower"]], 2), ", ",
                              round(g$edges[["ci.upper"]], 2), "]")
    } else {
      g$edges$label <- as.character(round(g$edges[["std.all"]], 2))
    }
  } else {
    g$edges$label <- ""
  }
  
  if (show_r2) {
    rsq <- tryCatch(inspect(fit, "r2"), error = function(e) NULL)
    if (!is.null(rsq)) {
      rsq_df <- data.frame(name = names(rsq), r2 = round(unlist(rsq), 2), stringsAsFactors = FALSE)
      if (!is.null(r2_nodes)) rsq_df <- rsq_df[rsq_df$name %in% r2_nodes, ]
      g$nodes <- merge(g$nodes, rsq_df, by = "name", all.x = TRUE)
      g$nodes$label <- ifelse(is.na(g$nodes$r2), g$nodes$label, paste0(g$nodes$label, "\nR² = ", g$nodes$r2))
    }
  }
  
  g$nodes$label_family <- g$edges$label_family <- font_family
  g$nodes$fill <- fill_c
  g$nodes$label_fill <- "transparent"
  
  gg <- plot(g) + 
    theme(text = element_text(family = font_family),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
  
  if (!is.null(title)) gg <- gg + labs(title = title)
  if (save) ggsave(filename = savename, plot = gg, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  if (display) print(gg) 
  
  return(gg) # Return the ggplot object so print() creates a real PDF page
}

# ==============================================================================
# 3. COMPILE ALL PAGES INTO A MULTI-PAGE PDF (EXPLICIT PRINTING)
# ==============================================================================
layout_dag1b <- get_layout(
  "PreyNCC", "", "",          "", "",       "", "PreyAK",
  "",        "", "Growth",    "", "",       "", "",
  "",        "", "",          "", "logSAR", "", "",
  "",        "", "Abundance", "", "",       "", "",
  "PredNCC", "", "",          "", "",       "", "PredAK",
  rows = 5
)

syntax_1b <- '
  Growth    ~ PreyNCC
  Abundance ~ PredNCC
  logSAR    ~ Growth + Abundance + PreyAK + PredAK
'

fig_dir <- dirname(pdf_out)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Make sure any dangling graphics devices are closed first
while (!is.null(dev.list())) dev.off()

# Open PDF device
pdf(pdf_out, width = 14, height = 8.5, onefile = TRUE, family = "serif")
cat(sprintf("\nWriting %d model pages to PDF: %s\n", nrow(combined_targets), pdf_out))

for (i in seq_len(nrow(combined_targets))) {
  row_entry <- combined_targets[i, ]
  mod_num   <- row_entry$modNum
  curr_aic  <- round(row_entry$AIC, 2)
  curr_daic <- round(row_entry$deltaAIC, 2)
  b_pak     <- round(row_entry$beta_preyak, 2)
  p_pak     <- round(row_entry$p_preyak, 3)
  grp       <- row_entry$group_label
  
  inds <- c(
    "PreyNCC"   = row_entry$PreyNCCindNames,
    "PredNCC"   = row_entry$PredNCCindNames,
    "Growth"    = row_entry$GrowthIndNames,
    "Abundance" = row_entry$AbundanceIndNames,
    "PreyAK"    = row_entry$PreyAKindNames,
    "PredAK"    = row_entry$PredAKindNames
  )
  
  df_mod <- pull_std_node(clusData_sub, "PreyNCC",   inds["PreyNCC"]) %>%
    full_join(pull_std_node(clusData_sub, "PredNCC",   inds["PredNCC"]),   by = "date") %>%
    full_join(pull_std_node(clusData_sub, "Growth",    inds["Growth"]),    by = "date") %>%
    full_join(pull_std_node(clusData_sub, "Abundance", inds["Abundance"]), by = "date") %>%
    full_join(pull_std_node(clusData_sub, "PreyAK",    inds["PreyAK"]),    by = "date") %>%
    full_join(pull_std_node(clusData_sub, "PredAK",    inds["PredAK"]),    by = "date") %>%
    full_join(pull_std_sar(clusData_sub), by = "date") %>%
    filter(!is.na(logSAR)) %>%
    arrange(date)
  
  fit_obj <- sem(syntax_1b, data = df_mod)
  
  labels_explicit <- c(
    "PreyNCC"   = build_explicit_label("Prey (NCC)", inds["PreyNCC"], cw),
    "PredNCC"   = build_explicit_label("Predators (NCC)", inds["PredNCC"], cw),
    "Growth"    = build_explicit_label("Salmon Growth", inds["Growth"], cw),
    "Abundance" = build_explicit_label("Salmon Abund.", inds["Abundance"], cw),
    "PreyAK"    = build_explicit_label("Prey (AK)", inds["PreyAK"], cw),
    "PredAK"    = build_explicit_label("Predators (AK)", inds["PredAK"], cw),
    "logSAR"    = "Salmon SAR\n[X16_SAR]"
  )
  
  plot_title <- sprintf("Page %d of %d | DAG 1B Long (24yr) - Mod #%d [%s]\nAIC = %.2f (dAIC = %.2f) | Prey(AK) beta = %.2f (p = %.3f)",
                        i, nrow(combined_targets), mod_num, grp, curr_aic, curr_daic, b_pak, p_pak)
  
  # Call plotting function and capture ggplot object
  p_gg <- Lisa_sem_graph.v3(
    fit               = fit_obj,
    title             = plot_title,
    layout            = layout_dag1b,
    node_labels       = labels_explicit,
    r2_nodes          = c("Growth", "Abundance", "logSAR"),
    p_threshold       = 0.05,
    nonsig_linetype   = "dashed",
    font_family       = "serif",
    rect_width        = 2.4,
    rect_height       = 0.9,
    ellipses_width    = 2.4,
    ellipses_height   = 0.9,
    text_size         = 3.2,
    variance_diameter = 0.6,
    curvature         = 10,
    width             = 14,
    height            = 8.5,
    display           = FALSE,
    save              = FALSE
  )
  
  # EXPLICIT PRINT TO PDF DEVICE
  print(p_gg)
  
  cat(sprintf("  -> Added Page %d: Mod #%d | PreyAK = %s (beta = %.2f)\n", 
              i, mod_num, inds["PreyAK"], b_pak))
}

# Close the PDF device properly
dev.off()
cat("\nProcess complete! Multi-page PDF written with valid pages to:\n", pdf_out, "\n")