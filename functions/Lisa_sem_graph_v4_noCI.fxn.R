# ==============================================================================
# File: functions/Lisa_sem_graph_v4_noCI.fxn.R
# ==============================================================================
library(tidyverse)
library(lavaan)
library(tidySEM)
library(ggplot2)

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
  
  # Outer indicator wings
  if (length(preyNCC_ind) > 0) mat[1, seq_along(preyNCC_ind)]          <- preyNCC_ind # Row 1, cols 1..7
  if (length(predNCC_ind) > 0) mat[1, 13 - seq_along(predNCC_ind) + 1] <- predNCC_ind # Row 1, cols 13..10
  if (length(preyAK_ind)  > 0) mat[9, seq_along(preyAK_ind) + 1]       <- preyAK_ind  # Row 9, cols 2..6
  if (length(predAK_ind)  > 0) mat[9, 13 - seq_along(predAK_ind) + 1]  <- predAK_ind  # Row 9, cols 13..9
  
  mat
}

Lisa_sem_graph_v4_noCI.fxn <- function(
    fit,
    layout,
    node_labels       = NULL,
    title             = NULL,
    r2_nodes          = c("X06_DFA_IGF_mu_smoltyr", "X07_cpue_IntSprJunHW_smoltyr", "X16_DFA_sar_UCsprCh_smoltyr"),
    growth_node       = "X06_DFA_IGF_mu_smoltyr",
    sar_node          = "X16_DFA_sar_UCsprCh_smoltyr",
    rect_width        = 1.7,
    rect_height       = 0.9,
    ellipses_width    = 2.4,
    ellipses_height   = 1.1,
    p_threshold       = 0.05,
    p_threshold2       = 0.1,
    font_family       = "serif",
    fill_c            = "grey96",
    save              = FALSE,
    savename          = "sem_plot.png",
    export_table_path = NULL,
    display           = TRUE,
    width             = 18,
    height            = 12,
    dpi               = 300
) {
  # 1. Build base graph using do.call so numeric dimensions evaluate cleanly------------
  g <- do.call(
    tidySEM::prepare_graph,
    list(
      model           = fit,
      layout          = layout,
      intercepts      = FALSE,
      rect_width      = rect_width,
      rect_height     = rect_height,
      ellipses_width  = ellipses_width,
      ellipses_height = ellipses_height
    )
  )
  
  # 2. Apply node labels + R² using base R (preserves tidySEM node attributes)-------------
  g$nodes$label <- g$nodes$name
  if (!is.null(node_labels)) {
    matches <- match(g$nodes$label, names(node_labels))
    g$nodes$label[!is.na(matches)] <- node_labels[matches[!is.na(matches)]]
  }
  
  if (!is.null(r2_nodes) && length(r2_nodes) > 0) {
    rsq <- tryCatch(inspect(fit, "r2"), error = function(e) NULL)
    if (!is.null(rsq)) {
      rsq_vals <- unlist(rsq[names(rsq) %in% r2_nodes])
      rsq_sub  <- setNames(sprintf("%.2f", rsq_vals), names(rsq_vals))
      r2_idx  <- match(names(rsq_sub), g$nodes$name)
      valid   <- !is.na(r2_idx)
      g$nodes$label[r2_idx[valid]] <- paste0(
        g$nodes$label[r2_idx[valid]], "\nR² = ", rsq_sub[valid]
      )
    }
  }
  
  g$nodes$label_family <- font_family
  g$nodes$fill         <- fill_c
  g$nodes$label_fill   <- "transparent"
  
  # 3. Remove empty self-variance circles (from == to) & match standardized estimates------------
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
  g$edges$linetype <- ifelse(
    !is.na(g$edges$pvalue) & g$edges$pvalue < p_threshold,  "solid",
    ifelse(!is.na(g$edges$pvalue) & g$edges$pvalue < p_threshold2, "dashed", "dotted")
  )
  g$edges$label        <- ifelse(is.na(g$edges$std.all), "", sprintf("%.2f", g$edges$std.all))
  g$edges$label_family <- font_family
  
  # 4. Route factor loadings vertically & force Growth -> SAR into top of SAR
  node_y     <- setNames(g$nodes$y, g$nodes$name)
  is_loading <- g$edges$op == "=~"
  from_y     <- node_y[g$edges$from]
  to_y       <- node_y[g$edges$to]
  
  g$edges$connect_from[is_loading & to_y > from_y] <- "top"
  g$edges$connect_to[is_loading & to_y > from_y]   <- "bottom"
  g$edges$connect_from[is_loading & to_y < from_y] <- "bottom"
  g$edges$connect_to[is_loading & to_y < from_y]   <- "top"
  
  is_growth_to_sar <- g$edges$op == "~" &
    g$edges$from == growth_node &
    g$edges$to   == sar_node
  
  if (any(is_growth_to_sar)) {
    g$edges$connect_from[is_growth_to_sar] <- "bottom"
    g$edges$connect_to[is_growth_to_sar]   <- "top"
  }
  
  # 5. Render ggplot object
  gg <- plot(g) +
    theme(
      text       = element_text(family = font_family),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  if (!is.null(title)) gg <- gg + labs(title = title)
  
  if (save) {
    ggsave(
      filename = savename,
      plot     = gg,
      width    = width,
      height   = height,
      units    = "in",
      dpi      = dpi,
      bg       = "white"
    )
  }
  
  # 6. Optional companion parameter table export
  if (!is.null(export_table_path)) {
    param_tbl <- pe_sub %>%
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
    write_csv(param_tbl, export_table_path)
  }
  
  if (display) print(gg)
  invisible(gg)
}