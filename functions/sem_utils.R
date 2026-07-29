


# load packages:
# create vector of packages needed
packages <- c("lavaan",
              "blavaan",
              "dplyr",
              "tidySEM",
              "semhelpinghands",
              "DT",
              "stringr",
              "ggplot2",
              "extrafont")

# loop through required packages,
# check if they are installed,
# install them if not,
# and load them
for (p in packages) {
  if (p %in% rownames(installed.packages())) {
    library(p, character.only=TRUE)
  } else {
    install.packages(p)
    library(p,character.only = TRUE)
  }
}



# summary
my_sem_summary <- function(fit, 
                           bootstrap = TRUE, 
                           digits = 3) {

  # --- Unstandardized estimates (will be percentile if bootstrap, SE based if not)---
  unstd <- parameterEstimates(fit, standardized = FALSE) %>%
    filter(op %in% c("~", "=~", "~~")) %>%
    select(lhs, op, rhs, est, se, ci.lower, ci.upper, pvalue)
  
  # --- Standardized estimates with bootstrap CIs (percentile) ---
  if (bootstrap) {
    std <- semhelpinghands::standardizedSolution_boot_ci(fit) %>%
      filter(op %in% c("~", "=~", "~~")) %>%
      select(lhs, op, rhs, est.std, boot.se, boot.ci.lower, boot.ci.upper, pvalue) %>%
      rename(est.std = est.std,
             se.std = boot.se,
             ci.lower.std = boot.ci.lower,
             ci.upper.std = boot.ci.upper,
             pvalue.std = pvalue)
  } else {
    # --- Standardized estimates with default standardized CIs (assume normality) ---
    std <- standardizedSolution(fit) %>%
      filter(op %in% c("~", "=~", "~~")) %>%
      select(lhs, op, rhs, est.std, se, ci.lower, ci.upper, pvalue) %>%
      rename(est.std = est.std,
             se.std = se,
             ci.lower.std = ci.lower,
             ci.upper.std = ci.upper,
             pvalue.std = pvalue)
  }
  
  # --- Merge estimates ---
  merged <- full_join(unstd, std, by = c("lhs", "op", "rhs"))
  
  # Round numeric columns
  num_cols <- sapply(merged, is.numeric)
  merged[num_cols] <- lapply(merged[num_cols], round, digits = digits)
  
  # --- R^2 values ---
  r2_tbl <- tryCatch({
    r2 <- inspect(fit, "r2")
    tibble(variable = names(r2), r2 = round(as.numeric(r2), digits))
  }, error = function(e) NULL)
  
  # --- Fit indices ---
  fitmeasures_tbl <- fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", "cfi","tli", "srmr"))
  fit_indices <- tibble::tibble(
    Measure = c("Chi-squared","df","Chi-squared p-value", "RMSEA", "CFI","TLI", "SRMR"),
    Value = round(as.numeric(fitmeasures_tbl), digits)
  )
  
  return(list(
    estimates = merged %>%
      mutate(
        op_sort = case_when(
          op == "=~" ~ 1,
          op == "~"  ~ 2,
          op == "~~" ~ 3,
          op == ":=" ~ 4,
          TRUE       ~ 5
        )
      ) %>%
      arrange(op_sort, str_to_lower(lhs), str_to_lower(rhs)) %>%
      select(-op_sort),
    r_squared = r2_tbl,
    fit_indices = fit_indices
  ))
}

my_bsem_summary <- function(fit,
                            digits = 3,
                            ci_level = 0.95) {
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(posterior)
  library(blavaan)
  
  alpha <- (1 - ci_level) / 2
  lower_prob <- alpha
  upper_prob <- 1 - alpha
  
  # ---- Extract Raw Posterior Draws ----
  stan_draws <- as.data.frame(fit@external$mcmcout)
  colnames(stan_draws) <- make.names(colnames(stan_draws))  # R-safe names
  
  # ---- Get Parameter Metadata ----
  param_info <- as.data.frame(fit@external$origpt) %>%
    select(lhs, op, rhs, pxnames, label) %>%
    filter(!is.na(pxnames)) %>%
    distinct(pxnames, .keep_all = TRUE)
  
  param_info_unstd <- param_info %>% filter(op != ":=")
  monitored_names <- make.names(param_info_unstd$pxnames)
  monitored_draws <- stan_draws %>% select(any_of(monitored_names))
  
  # ---- Convergence Diagnostics ----
  rhat_vals <- blavInspect(fit, "rhat")
  neff_vals <- blavInspect(fit, "neff")
  
  diagnostics_tbl <- tibble(
    parameter = names(rhat_vals),
    rhat = round(rhat_vals, digits),
    neff = round(neff_vals, digits)
  )
  
  # ---- Compute Stats for Monitored Parameters ----
  summary_stats <- monitored_draws %>%
    summarise(across(
      everything(),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        ci.lower = ~quantile(.x, probs = lower_prob, na.rm = TRUE),
        ci.upper = ~quantile(.x, probs = upper_prob, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ))

  ci_long <- summary_stats %>%
    pivot_longer(
      everything(),
      names_to = c("pxnames", "stat"),
      names_pattern = "(.+?)_(mean|median|sd|ci\\.lower|ci\\.upper)$",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value)

  unstd_monitored <- ci_long %>%
    mutate(pxnames = make.names(pxnames)) %>%
    left_join(param_info_unstd %>% mutate(pxnames = make.names(pxnames)), by = "pxnames") %>%
    select(lhs, op, rhs, est.mean = mean, est.median = median, sd, ci.lower, ci.upper)
  
  # ---- Compute Stats for := Parameters ----
  defined_params <- parameterEstimates(fit, standardized = FALSE) %>%
    filter(op == ":=") %>%
    select(lhs, rhs)
  
  label_map <- param_info %>%
    filter(!is.na(label), label != "", op != ":=") %>%
    select(label, pxnames) %>%
    mutate(pxnames = make.names(pxnames))
  
  if (nrow(defined_params) > 0) {
    defined_draws <- lapply(1:nrow(defined_params), function(i) {
      expr_str <- defined_params$rhs[i]
      for (j in seq_len(nrow(label_map))) {
        expr_str <- gsub(paste0("\\b", label_map$label[j], "\\b"),
                         paste0("`", label_map$pxnames[j], "`"), expr_str)
      }
      expr <- parse(text = expr_str)[[1]]
      eval(expr, envir = stan_draws)
    })
    
    names(defined_draws) <- defined_params$lhs
    defined_draws_df <- as.data.frame(defined_draws)
    
    defined_summary <- defined_draws_df %>%
      summarise(across(
        everything(),
        list(
          mean = ~mean(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE),
          sd = ~sd(.x, na.rm = TRUE),
          ci.lower = ~quantile(.x, probs = lower_prob, na.rm = TRUE),
          ci.upper = ~quantile(.x, probs = upper_prob, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )) %>%
      pivot_longer(
        everything(),
        names_to = c("lhs", "stat"),
        names_pattern = "(.+?)_(mean|median|sd|ci\\.lower|ci\\.upper)$",
        values_to = "value"
      ) %>%
      pivot_wider(names_from = stat, values_from = value) %>%
      mutate(
        op = ":=",
        rhs = defined_params$rhs[match(lhs, defined_params$lhs)]
      ) %>%
      select(lhs, op, rhs, est.mean = mean, est.median = median, sd, ci.lower, ci.upper)
  } else {
    defined_summary <- NULL
  }
  
  # ---- Combine Unstandardized Results ----
  unstandardized_results <- bind_rows(unstd_monitored, defined_summary)
  
  # ---- Append fixed parameters that were not monitored ----
  # Get all parameters (including fixed ones)
  pe_all <- parameterEstimates(fit, standardized = FALSE) %>%
    filter(op %in% c("=~", "~", "~~")) %>%
    mutate(px_id = paste(lhs, op, rhs, sep = "|"))
  
  # Identify which were monitored (i.e., appeared in posterior draws)
  monitored_ids <- unstandardized_results %>%
    mutate(px_id = paste(lhs, op, rhs, sep = "|")) %>%
    pull(px_id)
  
  # Get fixed parameters that were not in posterior draws
  fixed_missing <- pe_all %>%
    filter(!(px_id %in% monitored_ids)) %>%
    select(lhs, op, rhs, est) %>%
    transmute(
      lhs,
      op,
      rhs,
      est.mean   = est,
      est.median = est,
      sd         = 0,
      ci.lower   = est,
      ci.upper   = est
    )
  
  # Add fixed parameters to unstandardized results
  unstandardized_results <- bind_rows(unstandardized_results, fixed_missing)
  
  # ---- Standardized Posterior Summary ----
  standardized_draws <- as.data.frame(standardizedPosterior(fit))
  
  summary_std <- lapply(standardized_draws, function(x) {
    ci_vals <- quantile(x, probs = c(lower_prob, upper_prob), na.rm = TRUE)
    names(ci_vals) <- c("ci.lower.std", "ci.upper.std")
    c(
      est.mean.std = mean(x, na.rm = TRUE),
      est.median.std = median(x, na.rm = TRUE),
      sd.std = sd(x, na.rm = TRUE),
      ci_vals
    )
  })
  
  summary_std_df <- as.data.frame(do.call(rbind, summary_std))
  summary_std_df$parameter <- rownames(summary_std_df)
  
  parsed <- str_match(summary_std_df$parameter, "^(.+?)(=~|~~|~|:=)(.+)$")
  
  standardized_results <- summary_std_df %>%
    mutate(
      lhs = parsed[, 2],
      op = parsed[, 3],
      rhs = parsed[, 4]
    ) %>%
    select(lhs, op, rhs, est.mean.std, est.median.std, sd.std, ci.lower.std, ci.upper.std)
  
  # ---- Merge All Results ----
  estimates <- full_join(unstandardized_results, standardized_results,
                         by = c("lhs", "op", "rhs")) %>%
    mutate(
      op_sort = case_when(
        op == "=~" ~ 1,
        op == "~"  ~ 2,
        op == "~~" ~ 3,
        op == ":=" ~ 4,
        TRUE       ~ 5
      )
    ) %>%
    arrange(op_sort, str_to_lower(lhs), str_to_lower(rhs)) %>%
    select(-op_sort) %>%
    mutate(across(where(is.numeric), ~round(.x, digits)))
  
  # ---- Bayesian R² ----
  r2_vals <- blavInspect(fit, "r2")
  r2_tbl <- tibble(variable = names(r2_vals),
                   r2 = round(as.numeric(r2_vals), digits))
  
  # ---- Fit Indices ----
  fm <- suppressWarnings(fitMeasures(fit))
  fit_inds <- tibble(Measure = names(fm),
                     Value = round(as.numeric(fm), digits))
  
  # ---- Return all ----
  return(list(
    estimates = estimates,
    r_squared = r2_tbl,
    fit_indices = fit_inds,
    convergence = diagnostics_tbl
  ))
}


my_sem_graph <- function(fit,
                         title= NULL,     # NEW: Add title argument,
                         layout = NULL,
                         node_width = 8,   # New default
                         node_height = 1,  # New default
                         lavaan = TRUE,
                         blavaan = FALSE,
                         standardized = TRUE,
                         bootstrap = FALSE,
                         show_r2 = TRUE,
                         bayesian_ci_level = .95,
                         font_family = "Times New Roman",
                         fill_c="grey96",
                         save = FALSE,
                         savename = "sem_plot.png",
                         display = TRUE,
                         width = 6,
                         height = 4,
                         dpi = 300,
                         ...) {
  
  g <- prepare_graph(model = fit, layout = layout, ...)

  if (lavaan) {
    summary <- if (bootstrap) {
      my_sem_summary(fit, bootstrap = TRUE)
    } else {
      my_sem_summary(fit, bootstrap = FALSE)
    }
    result <- summary$estimates
  } else if (blavaan) {
    summary <- my_bsem_summary(fit, digits = 2, ci_level = bayesian_ci_level)
    result <- summary$estimates
  }

  # Select estimate columns
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
                           c("lhs", "op", "rhs", est_name,
                             lower_ci_name, upper_ci_name)]

  # Merge labeled results into g$edges
  g$edges <- merge(g$edges,
                   result_labeled,
                   by.x = c("rhs", "op", "lhs"),
                   by.y = c("rhs", "op", "lhs"),
                   all.x = TRUE,
                   suffixes = c(".x", ""))

  # Fix symmetric covariances
  cov_fix <- result_labeled[result_labeled$op == "~~" & result_labeled$lhs != result_labeled$rhs, ]
  for (i in seq_len(nrow(cov_fix))) {
    pair <- c(cov_fix$lhs[i], cov_fix$rhs[i])
    idx <- which(
      g$edges$op == "~~" &
        g$edges$from %in% pair &
        g$edges$to %in% pair &
        g$edges$from != g$edges$to
    )
    if (length(idx) == 1 && is.na(g$edges[[est_name]][idx])) {
      g$edges[[est_name]][idx] <- cov_fix[[est_name]][i]
      g$edges[[lower_ci_name]][idx] <- cov_fix[[lower_ci_name]][i]
      g$edges[[upper_ci_name]][idx] <- cov_fix[[upper_ci_name]][i]
    }
  }

  # Label with estimate + CI
  g$edges$label <- paste0(round(g$edges[[est_name]], 2), "\n[",
                          round(g$edges[[lower_ci_name]], 2), ", ",
                          round(g$edges[[upper_ci_name]], 2), "]")

  # Add font to edges
  g$edges$label_family <- font_family

  # Add R² values
  if (show_r2) {
    rsq <- tryCatch(inspect(fit, "r2"), error = function(e) NULL)
    if (!is.null(rsq)) {
      rsq_df <- data.frame(name = names(rsq),
                           r2 = round(unlist(rsq), 2),
                           stringsAsFactors = FALSE)
      g$nodes <- merge(g$nodes, rsq_df, by = "name", all.x = TRUE)
      g$nodes$label <- ifelse(is.na(g$nodes$r2),
                              g$nodes$label,
                              paste0(g$nodes$label, "\nR² = ", g$nodes$r2))
    }
  }

  # Add font to nodes
  g$nodes$label_family <- font_family
  
  # node fill color
  g$nodes$fill <- fill_c
  g$nodes$label_fill <- "transparent"

 
  # Finalizing the plot object
  gg <- plot(g) + 
      theme(
        text = element_text(family = font_family),
        # NEW: Center the title (hjust = 0.5)
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14) 
      )
  # NEW: Add title if provided
  if (!is.null(title)) {
    gg <- gg + labs(title = title)
  }
  
  # Save and/or display
  if (save) {
    # Use 'gg' which now contains the title
    ggsave(filename = savename, plot = gg,
           width = width, height = height, units = "in", 
           dpi = dpi, bg = "white", create.dir = TRUE)
  }
  
  if (display) {
    print(gg) # Print 'gg' instead of re-plotting g
  } 
 
  invisible(g)
}

my_sem_graph_unfit <- function(model_syntax,
                               layout,
                               save = TRUE,
                               display = TRUE,
                               savename = "sem_plot.png",
                               width = 6,
                               height = 4,
                               dpi = 300,
                               text_size=6,
                               font_family = "Times New Roman",
                               fill_c = "transparent",
                               ...) {
  requireNamespace("lavaan")
  requireNamespace("tidySEM")
  requireNamespace("ggplot2")
  requireNamespace("dplyr")
  
  #Parse lavaan syntax
  pt <- lavaan::lavaanify(model_syntax)
  
  # Extract edges
  
  # Regression paths:
  regression_edges <- pt %>%
    dplyr::filter(op == "~") %>%
    dplyr::transmute(
      from = rhs,
      to = lhs,
      label = ifelse(is.na(label), "", label),
      arrow = "last",
      linetype = "solid",
      curvature = NA,
      label_fill = "white",
      label_family = font_family,
      label_location = 0.4
    )
  
  # Latent measurement paths
  measurement_edges <- pt %>%
    dplyr::filter(op == "=~") %>%
    dplyr::transmute(
      from = lhs,  # latent variable
      to = rhs,    # indicator
      label = ifelse(is.na(label), "", label),
      arrow = "last",
      linetype = "solid",
      curvature = NA,
      label_fill = "white",
      label_family = font_family,
      label_location = 0.5
    )
  
  # Variances
  variance_edges <- pt %>%
    dplyr::filter(op == "~~", lhs == rhs, !is.na(label) & nzchar(label)) %>%
    dplyr::transmute(
      from = lhs,
      to = rhs,
      label = label,
      arrow = "both",
      linetype = "solid",
      curvature = NA,
      label_fill = "white",
      label_family = font_family,
      label_location = NA
    )
  
  # Covariances
  covariance_edges <- pt %>%
    dplyr::filter(op == "~~", lhs != rhs) %>%
    dplyr::transmute(
      from = lhs,
      to = rhs,
      label = ifelse(is.na(label), "", label),  # ← keep label if present
      arrow = "none",
      linetype = "dashed",
      curvature = 60,
      label_fill = "white",
      label_family = font_family,
      label_location = 0.5
    )
  
  # Combine edge components only if they are non-empty and valid
  edges <- list(
    regression_edges,
    measurement_edges,
    variance_edges,
    covariance_edges
  ) |>
    purrr::discard(~ nrow(.x) == 0) |>  # remove empty data frames
    purrr::map(~ dplyr::mutate(.x, label = as.character(label))) |>  # ensure label is character
    dplyr::bind_rows() 
  
  # Identify latent variables
  latent_vars <- pt %>%
    dplyr::filter(op == "=~") %>%
    dplyr::pull(lhs) %>%
    unique()
  
  all_vars <- unique(c(edges$from, edges$to))
  
  # Define node shapes and labels
  nodes <- data.frame(
    name = all_vars,
    shape = ifelse(all_vars %in% latent_vars, "oval", "rect"),
    label_family = font_family,
    label = "",  # blank labels; parsed separately below
    fill = fill_c,
    label_fill = "transparent",
    stringsAsFactors = FALSE
  )
  
  # Build tidySEM graph
  g <- tidySEM::prepare_graph(edges = edges, nodes = nodes, layout = layout, ...)
  
  # Node label parsing 
  node_labels <- g$nodes %>%
    dplyr::mutate(
      label = gsub("([a-zA-Z]+)([0-9]+)", "\\1[\\2]", name)
    )
  
  # Render original tidySEM plot
  p <- plot(g)
  
  # Extract tidySEM's automatic edge label positions
  built <- ggplot2::ggplot_build(p)

  edge_label_layer <- purrr::keep(
    built$data,
    ~ "label" %in% names(.x) && any(nzchar(.x$label)) && any(.x$label != " ")
  )
  
  if (is.null(edge_label_layer)) {
    warning("Could not find edge label layer.")
    df_edges_parsed <- NULL
  } else {
    df_edges_parsed <- edge_label_layer %>%
      dplyr::bind_rows() %>%
      dplyr::filter(nzchar(label)) %>%
      dplyr::mutate(label = gsub("([a-z]+)([0-9A-Za-z_]+)", "\\1[\\2]", label)) %>%
      dplyr::select(x, y, label)
  }
  
  # Force blank node/edge labels so we can overlay parsed ones
  g$nodes$label <- "   "
  g$edges$label <- "    "
  
  # create plot
  p <- plot(g)
  
  # Overlay parsed labels
  final_plot <- p +
    ggplot2::geom_text(data = node_labels,
                       ggplot2::aes(x = x, y = y, label = label),
                       parse = TRUE,
                       size = text_size) +
    ggplot2::geom_text(data = df_edges_parsed,
                       ggplot2::aes(x = x, y = y, label = label),
                       parse = TRUE,
                       size = text_size)
  
  # Save figure
  if (save) {
    grDevices::png(filename = savename,
                   width = width,
                   height = height,
                   units = "in",
                   res = dpi,
                   bg = "white")
    grid::grid.draw(final_plot)
    grDevices::dev.off()
  }
  
  # Display plot
  if (display) {
    print(final_plot)
  }
  
  invisible(final_plot)
}


my_scatterplot_gg <- function(x, y, xlab, ylab, xlim, ylim,
                              col = "darkblue", bg = "deepskyblue",
                              savename = 'plot.png',
                              save = TRUE,
                              display = TRUE,
                              width = 3,
                              height = 3,
                              line=FALSE,
                              equal_aspect = TRUE,
                              dpi=400) {

  df <- data.frame(x = x, y = y)

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(shape = 21, color = col, stroke=1,fill = bg, size = 2) +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_family = "Times New Roman",
                  base_size = 12) +
    theme(
      panel.grid.major = element_line(color = alpha("black", .3), linetype = "11", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks = element_line(linewidth=.5,color = "black"),
      axis.text = element_text(color = "black", family = "Times New Roman"),
      axis.title = element_text(color = "black",family = "Times New Roman"),
      plot.margin = margin(t = 15, r = 15, b = 5, l = 5)
    )

  if (line) {
    p <- p+geom_line(color=bg)
  }
  
  if (equal_aspect) {
    p <- p+theme(aspect.ratio=1)
  }
  
  if (save) {
    ggsave(filename = savename,
           plot = p,
           width = width,
           height = height,
           units = "in",
           dpi = dpi,
           bg = "white",
           create.dir = TRUE)
  }
  if (display) {
    print(p)
  }
}




