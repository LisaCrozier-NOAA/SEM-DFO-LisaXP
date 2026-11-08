

#Here is a complete, automated script to identify, print, fit, and plot all acceptable models ($\Delta\text{AIC} \le 4$, $p_{\chi^2} \ge 0.05$, $\text{CFI} \ge 0.90$, $\text{SRMR} \le 0.08$) for both DAG 1A (Long) and DAG 1B (Long).The script iterates through every winning row, dynamically extracts and standardizes the exact indicator data from clusDataDFA.rds, builds the explicit LisaNames via crosswalk_names.csv, and outputs each diagram with dashed lines for non-significant paths ($p \ge 0.05$).
#Found 12 acceptable models for DAG 1A and 4 for DAG 1B


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
dir_1a    <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1A_long"
dir_1b    <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1B_long"

clusData <- readRDS(clus_path)
cw       <- read.csv(cw_path, stringsAsFactors = FALSE)
res_1a   <- read.csv(file.path(dir_1a, "SEMresultsByClus.csv"), stringsAsFactors = FALSE)
res_1b   <- read.csv(file.path(dir_1b, "SEMresultsByClus.csv"), stringsAsFactors = FALSE)

# 24-year long dataset (1998-2021)
clusData_sub <- clusData %>%
  filter(minYear == 1998) %>%
  mutate(year = year(ymd(date))) %>%
  filter(year >= 1998, year <= 2021)

# ==============================================================================
# 1. FILTER ACCEPTABLE MODELS (Delta AIC <= 4 & Adequate Fit)
# ==============================================================================
get_acceptable <- function(res_df, case_name) {
  min_aic <- min(res_df$AIC, na.rm = TRUE)
  res_df %>%
    mutate(
      model_case = case_name,
      deltaAIC   = round(AIC - min_aic, 2)
    ) %>%
    filter(deltaAIC <= 4.0, p_Chi2 >= 0.05, CFI >= 0.90, SRMR <= 0.08) %>%
    arrange(AIC)
}

acc_1a <- get_acceptable(res_1a, "DAG1A_long")
acc_1b <- get_acceptable(res_1b, "DAG1B_long")

cat("====================================================================\n")
cat(sprintf("Found %d acceptable models for DAG 1A and %d for DAG 1B\n", nrow(acc_1a), nrow(acc_1b)))
cat("====================================================================\n")

print(acc_1a %>% select(modNum, AIC, deltaAIC, p_Chi2, CFI, SRMR, RMSEA))
print(acc_1b %>% select(modNum, AIC, deltaAIC, p_Chi2, CFI, SRMR, RMSEA))

# ==============================================================================
# 2. CORE UTILITY FUNCTIONS
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

# Robust Lisa_sem_graph.v3
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
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  if (!is.null(title)) gg <- gg + labs(title = title)
  if (save) ggsave(filename = savename, plot = gg, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  if (display) print(gg) 
  invisible(g)
}

# ==============================================================================
# 3. REUSABLE PLOTTING WRAPPER
# ==============================================================================
fit_and_plot_model <- function(row_entry, dag_type, rank_idx) {
  
  mod_num  <- row_entry$modNum
  curr_aic <- round(row_entry$AIC, 2)
  curr_daic<- round(row_entry$deltaAIC, 2)
  
  inds <- c(
    "PreyNCC"   = row_entry$PreyNCCindNames,
    "PredNCC"   = row_entry$PredNCCindNames,
    "Growth"    = row_entry$GrowthIndNames,
    "Abundance" = row_entry$AbundanceIndNames,
    "PreyAK"    = row_entry$PreyAKindNames,
    "PredAK"    = row_entry$PredAKindNames
  )
  
  # Assemble standardized data
  df_mod <- pull_std_node(clusData_sub, "PreyNCC",   inds["PreyNCC"]) %>%
    full_join(pull_std_node(clusData_sub, "PredNCC",   inds["PredNCC"]),   by = "date") %>%
    full_join(pull_std_node(clusData_sub, "Growth",    inds["Growth"]),    by = "date") %>%
    full_join(pull_std_node(clusData_sub, "Abundance", inds["Abundance"]), by = "date") %>%
    full_join(pull_std_node(clusData_sub, "PreyAK",    inds["PreyAK"]),    by = "date") %>%
    full_join(pull_std_node(clusData_sub, "PredAK",    inds["PredAK"]),    by = "date") %>%
    full_join(pull_std_sar(clusData_sub), by = "date") %>%
    filter(!is.na(logSAR)) %>%
    arrange(date)
  
  # Structural syntax & Layouts
  if (dag_type == "DAG1A") {
    syntax <- '
      Growth    ~ PreyNCC
      Abundance ~ Growth + PredNCC
      logSAR    ~ Abundance + PreyAK + PredAK
    '
    layout <- get_layout(
      "PreyNCC", "", "",          "", "",       "", "PreyAK",
      "",        "", "Growth",    "", "",       "", "",
      "",        "", "Abundance", "", "logSAR", "", "",
      "PredNCC", "", "",          "", "",       "", "PredAK",
      rows = 4
    )
    plt_h <- 8
  } else {
    syntax <- '
      Growth    ~ PreyNCC
      Abundance ~ PredNCC
      logSAR    ~ Growth + Abundance + PreyAK + PredAK
    '
    layout <- get_layout(
      "PreyNCC", "", "",          "", "",       "", "PreyAK",
      "",        "", "Growth",    "", "",       "", "",
      "",        "", "",          "", "logSAR", "", "",
      "",        "", "Abundance", "", "",       "", "",
      "PredNCC", "", "",          "", "",       "", "PredAK",
      rows = 5
    )
    plt_h <- 9
  }
  
  fit_obj <- sem(syntax, data = df_mod)
  
  labels_explicit <- c(
    "PreyNCC"   = build_explicit_label("Prey (NCC)", inds["PreyNCC"], cw),
    "PredNCC"   = build_explicit_label("Predators (NCC)", inds["PredNCC"], cw),
    "Growth"    = build_explicit_label("Salmon Growth", inds["Growth"], cw),
    "Abundance" = build_explicit_label("Salmon Abund.", inds["Abundance"], cw),
    "PreyAK"    = build_explicit_label("Prey (AK)", inds["PreyAK"], cw),
    "PredAK"    = build_explicit_label("Predators (AK)", inds["PredAK"], cw),
    "logSAR"    = "Salmon SAR\n[X16_SAR]"
  )
  
  plot_title <- sprintf("%s Long (24yr) - Rank #%d (Mod #%d, AIC = %.2f, dAIC = %.2f)",
                        dag_type, rank_idx, mod_num, curr_aic, curr_daic)
  
  Lisa_sem_graph.v3(
    fit               = fit_obj,
    title             = plot_title,
    layout            = layout,
    node_labels       = labels_explicit,
    r2_nodes          = c("Growth", "Abundance", "logSAR"),
    p_threshold       = 0.05,
    nonsig_linetype   = "dashed",
    font_family       = "serif",
    rect_width        = 2.4, rect_height = 0.9,
    ellipses_width    = 2.4, ellipses_height = 0.9,
    text_size         = 3.2, variance_diameter = 0.6, curvature = 10,
    width             = 14, height = plt_h,
    display           = TRUE,
    save              = FALSE
  )
}

# ==============================================================================
# 4. RENDER ALL ACCEPTABLE MODELS
# ==============================================================================
# Plot all acceptable DAG 1B models
cat("\nRendering DAG 1B acceptable models...\n")
for (i in seq_len(nrow(acc_1b))) {
  fit_and_plot_model(acc_1b[i, ], "DAG1B", i)
}

# Plot all acceptable DAG 1A models
cat("\nRendering DAG 1A acceptable models...\n")
if (nrow(acc_1a) > 0) {
  for (i in seq_len(nrow(acc_1a))) {
    fit_and_plot_model(acc_1a[i, ], "DAG1A", i)
  }
} else {
  cat("No models met the strict acceptability criteria (p_Chi2 >= 0.05, CFI >= 0.90, SRMR <= 0.08) for DAG 1A.\n")
}
