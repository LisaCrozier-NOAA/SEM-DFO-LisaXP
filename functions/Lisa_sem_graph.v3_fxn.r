library(tidyverse)
library(tidySEM)
library(ggplot2)

# 1. Ensure helper generates CI columns
my_sem_summary <- function(fit, bootstrap = FALSE) {
  # ci = TRUE guarantees ci.lower and ci.upper exist
  pe <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
  list(estimates = pe)
}

# 2. Robust Lisa_sem_graph.v3
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
  
  # Prepare graph
  g <- prepare_graph(model = fit, layout = layout, intercepts = FALSE, ...)
  g$nodes$label <- g$nodes$name 
  
  if (!is.null(node_labels)) {
    matches <- match(g$nodes$label, names(node_labels))
    g$nodes$label[!is.na(matches)] <- node_labels[matches[!is.na(matches)]]
  }
  
  # Extract estimates safely with CIs
  result <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
  
  # Verify required columns exist; fill with NA if missing
  required_cols <- c("lhs", "op", "rhs", "std.all", "ci.lower", "ci.upper", "pvalue")
  for (col in required_cols) {
    if (!col %in% names(result)) {
      result[[col]] <- NA_real_
    }
  }
  
  result_labeled <- result[result$op %in% c("~", "=~", "~~") & result$lhs != result$rhs, 
                           required_cols]
  
  # Merge with graph edge definitions
  g$edges <- merge(g$edges, result_labeled, by.x = c("rhs", "op", "lhs"),
                   by.y = c("rhs", "op", "lhs"), all.x = FALSE)
  
  # Determine line type from significance
  g$edges$linetype <- ifelse(!is.na(g$edges$pvalue) & g$edges$pvalue < p_threshold, 
                             sig_linetype, 
                             nonsig_linetype)
  
  # Label coefficients
  if (show_coefficients) {
    # If CIs are present, show estimate and interval; otherwise just estimate
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
  
  # R² logic
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