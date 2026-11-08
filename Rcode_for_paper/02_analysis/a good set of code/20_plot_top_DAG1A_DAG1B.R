library(tidyverse)
library(lubridate)
library(lavaan)
library(tidySEM)
library(ggplot2)

# ==============================================================================
# 0. SETUP PATHS
# ==============================================================================
clus_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DFA/clusDataDFA.rds"
cw_path   <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/crosswalk_names.csv"

clusData <- readRDS(clus_path)
cw       <- read.csv(cw_path, stringsAsFactors = FALSE)

# Filter 24-year long dataset (1998-2021)
clusData_sub <- clusData %>%
  filter(minYear == 1998) %>%
  mutate(year = year(ymd(date))) %>%
  filter(year >= 1998, year <= 2021)

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================
# Pull, standardize, and rename a single indicator series
pull_std_node <- function(cdat, var_node, ind_name) {
  df_node <- cdat %>%
    filter(SEMlatent == var_node, shortName == ind_name) %>%
    select(date, finalVal) %>%
    arrange(date)
  
  df_node[[var_node]] <- as.vector(scale(df_node$finalVal))
  df_node %>% select(date, all_of(var_node))
}

# Pull and standardize logSAR
pull_std_sar <- function(cdat) {
  cdat %>%
    filter(SEMlatent == "SAR") %>%
    select(date, finalVal) %>%
    arrange(date) %>%
    mutate(logSAR = as.vector(scale(finalVal))) %>%
    select(date, logSAR)
}

# Lookup explicit LisaNames from crosswalk
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

# My SEM Summary fallback
if (!exists("my_sem_summary")) {
  my_sem_summary <- function(fit, bootstrap = FALSE) {
    list(estimates = parameterEstimates(fit, standardized = TRUE))
  }
}

# ==============================================================================
# 2. PLOTTING FUNCTION (Lisa_sem_graph.v3)
# ==============================================================================
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
  
  summary_obj <- my_sem_summary(fit, bootstrap = FALSE)
  result <- summary_obj$estimates
  
  result_labeled <- result[result$op %in% c("~", "=~", "~~") & result$lhs != result$rhs, 
                           c("lhs", "op", "rhs", "std.all", "ci.lower", "ci.upper", "pvalue")]
  
  g$edges <- merge(g$edges, result_labeled, by.x = c("rhs", "op", "lhs"),
                   by.y = c("rhs", "op", "lhs"), all.x = FALSE)
  
  g$edges$linetype <- ifelse(!is.na(g$edges$pvalue) & g$edges$pvalue < p_threshold, 
                             sig_linetype, 
                             nonsig_linetype)
  
  if (show_coefficients) {
    g$edges$label <- paste0(round(g$edges[["std.all"]], 2), "\n[",
                            round(g$edges[["ci.lower"]], 2), ", ",
                            round(g$edges[["ci.upper"]], 2), "]")
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
# 3. FIT & PLOT: DAG 1A LONG (24yr)
# ==============================================================================
dag1a_inds <- c(
  "PreyNCC"   = "04.CompFishNCC_smoothed",
  "PredNCC"   = "09.PredFishNCC_b_DFA1",
  "Growth"    = "06.Cond1NCC_DFA1",
  "Abundance" = "07.Cond2NCC_smoothed",
  "PreyAK"    = "13.FishPreyAK_DFA1",
  "PredAK"    = "10.PredMammalNCC_0_smoothed"
)

# Assemble DAG 1A data
df_dag1a <- pull_std_node(clusData_sub, "PreyNCC",   dag1a_inds["PreyNCC"]) %>%
  full_join(pull_std_node(clusData_sub, "PredNCC",   dag1a_inds["PredNCC"]),   by = "date") %>%
  full_join(pull_std_node(clusData_sub, "Growth",    dag1a_inds["Growth"]),    by = "date") %>%
  full_join(pull_std_node(clusData_sub, "Abundance", dag1a_inds["Abundance"]), by = "date") %>%
  full_join(pull_std_node(clusData_sub, "PreyAK",    dag1a_inds["PreyAK"]),    by = "date") %>%
  full_join(pull_std_node(clusData_sub, "PredAK",    dag1a_inds["PredAK"]),    by = "date") %>%
  full_join(pull_std_sar(clusData_sub), by = "date") %>%
  filter(!is.na(logSAR)) %>%
  arrange(date)

# Model syntax DAG 1A
syntax_1a <- '
  Growth    ~ PreyNCC
  Abundance ~ Growth + PredNCC
  logSAR    ~ Abundance + PreyAK + PredAK
'
fit_1a <- sem(syntax_1a, data = df_dag1a)

# Labels & Layout DAG 1A
labels_1a <- c(
  "PreyNCC"   = build_explicit_label("Prey (NCC)", dag1a_inds["PreyNCC"], cw),
  "PredNCC"   = build_explicit_label("Predators (NCC)", dag1a_inds["PredNCC"], cw),
  "Growth"    = build_explicit_label("Salmon Growth", dag1a_inds["Growth"], cw),
  "Abundance" = build_explicit_label("Salmon Abund.", dag1a_inds["Abundance"], cw),
  "PreyAK"    = build_explicit_label("Prey (AK)", dag1a_inds["PreyAK"], cw),
  "PredAK"    = build_explicit_label("Predators (AK)", dag1a_inds["PredAK"], cw),
  "logSAR"    = "Salmon SAR\n[X16_SAR]"
)

layout_1a <- get_layout(
  "PreyNCC", "", "",          "", "",       "", "PreyAK",
  "",        "", "Growth",    "", "",       "", "",
  "",        "", "Abundance", "", "logSAR", "", "",
  "PredNCC", "", "",          "", "",       "", "PredAK",
  rows = 4
)

# Plot DAG 1A
p_1a <- Lisa_sem_graph.v3(
  fit               = fit_1a,
  title             = "DAG 1A Long (24yr) - Top Model",
  layout            = layout_1a,
  node_labels       = labels_1a,
  r2_nodes          = c("Growth", "Abundance", "logSAR"),
  p_threshold       = 0.05,
  nonsig_linetype   = "dashed",
  font_family       = "serif",
  rect_width        = 2.4, rect_height = 0.9,
  ellipses_width    = 2.4, ellipses_height = 0.9,
  text_size         = 3.2, variance_diameter = 0.6, curvature = 10,
  width             = 14, height = 8
)

# ==============================================================================
# 4. FIT & PLOT: DAG 1B LONG (24yr)
# ==============================================================================
dag1b_inds <- c(
  "PreyNCC"   = "04.CompFishNCC_smoothed",
  "PredNCC"   = "09.PredFishNCC_b_DFA1",
  "Growth"    = "06.Cond1NCC_DFA1",
  "Abundance" = "07.Cond2NCC_smoothed",
  "PreyAK"    = "12.ZooPreyAK_0_smoothed",
  "PredAK"    = "15.PredFishAK_2_smoothed"
)

# Assemble DAG 1B data
df_dag1b <- pull_std_node(clusData_sub, "PreyNCC",   dag1b_inds["PreyNCC"]) %>%
  full_join(pull_std_node(clusData_sub, "PredNCC",   dag1b_inds["PredNCC"]),   by = "date") %>%
  full_join(pull_std_node(clusData_sub, "Growth",    dag1b_inds["Growth"]),    by = "date") %>%
  full_join(pull_std_node(clusData_sub, "Abundance", dag1b_inds["Abundance"]), by = "date") %>%
  full_join(pull_std_node(clusData_sub, "PreyAK",    dag1b_inds["PreyAK"]),    by = "date") %>%
  full_join(pull_std_node(clusData_sub, "PredAK",    dag1b_inds["PredAK"]),    by = "date") %>%
  full_join(pull_std_sar(clusData_sub), by = "date") %>%
  filter(!is.na(logSAR)) %>%
  arrange(date)

# Model syntax DAG 1B: Converging Model
# Growth goes directly to SAR; Abundance goes directly to SAR[cite: 1]
syntax_1b <- '
  Growth    ~ PreyNCC
  Abundance ~ PredNCC
  logSAR    ~ Growth + Abundance + PreyAK + PredAK
'
fit_1b <- sem(syntax_1b, data = df_dag1b)

# Print parameter estimates for inspection
cat("\n--- DAG 1B Long Parameter Estimates ---\n")
parameterEstimates(fit_1b, standardized = TRUE) %>%
  filter(op == "~") %>%
  select(lhs, op, rhs, est, se, z, pvalue, std.all) %>%
  print()

# Labels & Layout DAG 1B
labels_1b <- c(
  "PreyNCC"   = build_explicit_label("Prey (NCC)", dag1b_inds["PreyNCC"], cw),
  "PredNCC"   = build_explicit_label("Predators (NCC)", dag1b_inds["PredNCC"], cw),
  "Growth"    = build_explicit_label("Salmon Growth", dag1b_inds["Growth"], cw),
  "Abundance" = build_explicit_label("Salmon Abund.", dag1b_inds["Abundance"], cw),
  "PreyAK"    = build_explicit_label("Prey (AK)", dag1b_inds["PreyAK"], cw),
  "PredAK"    = build_explicit_label("Predators (AK)", dag1b_inds["PredAK"], cw),
  "logSAR"    = "Salmon SAR\n[X16_SAR]"
)

# Symmetrical layout with SAR centered vertically
layout_1b <- get_layout(
  "PreyNCC", "", "",          "", "",       "", "PreyAK",
  "",        "", "Growth",    "", "",       "", "",
  "",        "", "",          "", "logSAR", "", "",
  "",        "", "Abundance", "", "",       "", "",
  "PredNCC", "", "",          "", "",       "", "PredAK",
  rows = 5
)

# Plot DAG 1B
p_1b <- Lisa_sem_graph.v3(
  fit               = fit_1b,
  title             = "DAG 1B Long (24yr) - Top Model",
  layout            = layout_1b,
  node_labels       = labels_1b,
  r2_nodes          = c("Growth", "Abundance", "logSAR"),
  p_threshold       = 0.05,
  nonsig_linetype   = "dashed",
  font_family       = "serif",
  rect_width        = 2.4, rect_height = 0.9,
  ellipses_width    = 2.4, ellipses_height = 0.9,
  text_size         = 3.2, variance_diameter = 0.6, curvature = 10,
  width             = 14, height = 9
)
