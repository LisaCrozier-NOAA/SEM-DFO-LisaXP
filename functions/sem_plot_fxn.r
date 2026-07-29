Lisa_sem_graph.v2 <- function(fit,
                              title = NULL,
                              layout = NULL,
                              node_labels = NULL, 
                              show_coefficients = TRUE, # NEW: Toggle for numeric labels
                              node_width = 8,
                              node_height = 1,
                              lavaan = TRUE,
                              blavaan = FALSE,
                              standardized = TRUE,
                              bootstrap = FALSE,
                              show_r2 = TRUE,
                              r2_nodes = NULL, 
                              bayesian_ci_level = .95,
                              font_family = "Times New Roman",
                              fill_c = "grey96",
                              save = FALSE,
                              savename = "sem_plot.png",
                              display = TRUE,
                              width = 6,
                              height = 4,
                              dpi = 300,
                              ...) {
  
  # 1. Prepare graph
  g <- prepare_graph(model = fit, layout = layout, intercepts = FALSE, ...)
  g$nodes$label <- g$nodes$name 
  
  if (!is.null(node_labels)) {
    matches <- match(g$nodes$label, names(node_labels))
    g$nodes$label[!is.na(matches)] <- node_labels[matches[!is.na(matches)]]
  }
  
  # 2. Extract results
  if (lavaan) {
    summary_obj <- if (bootstrap) my_sem_summary(fit, bootstrap = TRUE) else my_sem_summary(fit, bootstrap = FALSE)
    result <- summary_obj$estimates
  } else if (blavaan) {
    summary_obj <- my_bsem_summary(fit, digits = 2, ci_level = bayesian_ci_level)
    result <- summary_obj$estimates
  }
  
  est_name <- if (standardized) (if (blavaan) "est.median.std" else "est.std") else (if (blavaan) "est.median" else "est")
  lower_ci_name <- if (standardized) "ci.lower.std" else "ci.lower"
  upper_ci_name <- if (standardized) "ci.upper.std" else "ci.upper"
  
  result_labeled <- result[result$op %in% c("~", "=~", "~~") & result$lhs != result$rhs, 
                           c("lhs", "op", "rhs", est_name, lower_ci_name, upper_ci_name)]
  
  g$edges <- merge(g$edges, result_labeled, by.x = c("rhs", "op", "lhs"),
                   by.y = c("rhs", "op", "lhs"), all.x = FALSE)
  
  # --- MODIFIED LABEL LOGIC ---
  if (show_coefficients) {
    g$edges$label <- paste0(round(g$edges[[est_name]], 2), "\n[",
                            round(g$edges[[lower_ci_name]], 2), ", ",
                            round(g$edges[[upper_ci_name]], 2), "]")
  } else {
    g$edges$label <- "" # Set to empty string to suppress
  }
  # ----------------------------
  
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



Lisa_sem_graph <- function(fit,
                         title = NULL,
                         layout = NULL,
                         node_labels = NULL, # NEW: Named vector for pretty names
                         node_width = 8,
                         node_height = 1,
                         lavaan = TRUE,
                         blavaan = FALSE,
                         standardized = TRUE,
                         bootstrap = FALSE,
                         show_r2 = TRUE,
                         r2_nodes = NULL, 
                         bayesian_ci_level = .95,
                         font_family = "Times New Roman",
                         fill_c = "grey96",
                         save = FALSE,
                         savename = "sem_plot.png",
                         display = TRUE,
                         width = 6,
                         height = 4,
                         dpi = 300,
                         ...) {
  
  # 1. Prepare graph
  g <- prepare_graph(model = fit, layout = layout, intercepts = FALSE, ...)
  
  # 2. Reset labels to node names (Cleans intercepts)
  g$nodes$label <- g$nodes$name 
  
  # --- NEW: SWAP FOR PRETTY NAMES ---
  if (!is.null(node_labels)) {
    # Match the names in the graph to the names in your translation vector
    matches <- match(g$nodes$label, names(node_labels))
    # Replace only the ones that have a match
    g$nodes$label[!is.na(matches)] <- node_labels[matches[!is.na(matches)]]
  }
  # ----------------------------------
  
  # (Summary & Edge logic remains same as before...)
  if (lavaan) {
    summary_obj <- if (bootstrap) my_sem_summary(fit, bootstrap = TRUE) else my_sem_summary(fit, bootstrap = FALSE)
    result <- summary_obj$estimates
  } else if (blavaan) {
    summary_obj <- my_bsem_summary(fit, digits = 2, ci_level = bayesian_ci_level)
    result <- summary_obj$estimates
  }
  
  if (standardized) {
    est_name <- if (blavaan) "est.median.std" else "est.std"
    lower_ci_name <- "ci.lower.std"
    upper_ci_name <- "ci.upper.std"
  } else {
    est_name <- if (blavaan) "est.median" else "est"
    lower_ci_name <- "ci.lower"
    upper_ci_name <- "ci.upper"
  }
  
  result_labeled <- result[result$op %in% c("~", "=~", "~~"),
                           c("lhs", "op", "rhs", est_name, lower_ci_name, upper_ci_name)]
  
  
  # NEW: Remove "Self-loops" (Variances/Residual Variances)
  # We only want to keep rows where the LHS and RHS are different
  result_labeled <- result_labeled[result_labeled$lhs != result_labeled$rhs, ]
  
  
  
  g$edges <- merge(g$edges, result_labeled, by.x = c("rhs", "op", "lhs"),
                   by.y = c("rhs", "op", "lhs"), all.x = FALSE, suffixes = c(".x", ""))
  
  # Fix symmetric covariances
  cov_fix <- result_labeled[result_labeled$op == "~~" & result_labeled$lhs != result_labeled$rhs, ]
  for (i in seq_len(nrow(cov_fix))) {
    pair <- c(cov_fix$lhs[i], cov_fix$rhs[i])
    idx <- which(g$edges$op == "~~" & g$edges$from %in% pair & g$edges$to %in% pair & g$edges$from != g$edges$to)
    if (length(idx) == 1 && is.na(g$edges[[est_name]][idx])) {
      g$edges[[est_name]][idx] <- cov_fix[[est_name]][i]
      g$edges[[lower_ci_name]][idx] <- cov_fix[[lower_ci_name]][i]
      g$edges[[upper_ci_name]][idx] <- cov_fix[[upper_ci_name]][i]
    }
  }
  
  g$edges$label <- paste0(round(g$edges[[est_name]], 2), "\n[",
                          round(g$edges[[lower_ci_name]][idx], 2), ", ", # Fix indexing if needed
                          round(g$edges[[upper_ci_name]][idx], 2), "]")
  
  # Simplified Edge Labeling to ensure it works
  g$edges$label <- paste0(round(g$edges[[est_name]], 2), "\n[",
                          round(g$edges[[lower_ci_name]], 2), ", ",
                          round(g$edges[[upper_ci_name]], 2), "]")
  
  g$edges$label_family <- font_family
  
  # 3. R² logic (Appends to the NEW pretty labels)
  if (show_r2) {
    rsq <- tryCatch(inspect(fit, "r2"), error = function(e) NULL)
    if (!is.null(rsq)) {
      rsq_df <- data.frame(name = names(rsq),
                           r2 = round(unlist(rsq), 2),
                           stringsAsFactors = FALSE)
      if (!is.null(r2_nodes)) rsq_df <- rsq_df[rsq_df$name %in% r2_nodes, ]
      
      g$nodes <- merge(g$nodes, rsq_df, by = "name", all.x = TRUE)
      g$nodes$label <- ifelse(is.na(g$nodes$r2),
                              g$nodes$label,
                              paste0(g$nodes$label, "\nR² = ", g$nodes$r2))
    }
  }
  
  g$nodes$label_family <- font_family
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