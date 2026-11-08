library(tidyverse)
library(lavaan)

# Paths to Step 1 Doug sweeps
sweep_files <- list(
  "DAG1A_long"  = "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1A_long/SEMresultsByClus.csv",
  "DAG1B_long"  = "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1B_long/SEMresultsByClus.csv",
  "DAG1A_short" = "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1A_short/SEMresultsByClus.csv",
  "DAG1B_short" = "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1B_short/SEMresultsByClus.csv"
)

# Extract winning indicator combinations for Step 1
top_step1_models <- imap_dfr(sweep_files, function(fpath, label) {
  if (!file.exists(fpath)) return(tibble())
  read.csv(fpath, stringsAsFactors = FALSE) %>%
    slice(1) %>%
    mutate(model_label = label) %>%
    dplyr::select(model_label, modNum, AIC, p_Chi2, CFI, SRMR, RMSEA,
                  PreyNCCindNames, PredNCCindNames, GrowthIndNames, 
                  AbundanceIndNames, PreyAKindNames, PredAKindNames)
})

print(top_step1_models)


crosswalk_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/crosswalk_names.csv"
crosswalk_df   <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

# Build a lookup vector: raw name -> clean name
make_label_lookup <- function(cw_df) {
  lookup <- c(
    # Structural Latent Nodes
    "PreyNCC"    = "Prey\n(NCC)",
    "PredNCC"    = "Predators\n(NCC)",
    "Growth"     = "Salmon\nGrowth",
    "Abundance"  = "Salmon\nAbundance",
    "CPUE"       = "Salmon\nAbundance",
    "PreyAK"     = "Prey\n(AK)",
    "PredAK"     = "Predators\n(AK)",
    "SAR"        = "Salmon\nSAR",
    "logSAR"     = "Salmon\nSAR",
    "P_eff_test" = "P_eff (Mod.)",
    "EMBssl"     = "Steller\nSea Lion"
  )
  
  for (i in seq_len(nrow(cw_df))) {
    clean_lbl <- ifelse(!is.na(cw_df$LisaName[i]) && cw_df$LisaName[i] != "", 
                        cw_df$LisaName[i], 
                        cw_df$indicator1[i])
    
    # Store exact DFA names and X-prefixed column names
    if (!is.na(cw_df$DFAname[i]))  lookup[cw_df$DFAname[i]]  <- clean_lbl
    if (!is.na(cw_df$dfa_cols[i])) lookup[cw_df$dfa_cols[i]] <- clean_lbl
    if (!is.na(cw_df$indicator1[i])) lookup[cw_df$indicator1[i]] <- clean_lbl
  }
  return(lookup)
}

node_label_map <- make_label_lookup(crosswalk_df)



# Retrieve model object and plot with automated labels
plot_winning_step1_model <- function(dag_name, run_dir, top_df, label_map) {
  row_match <- top_df %>% filter(model_label == dag_name)
  mod_num   <- row_match$modNum[1]
  
  # Load saved lavaan fit
  mod_path <- file.path(run_dir, "doPar", paste0("modList_", mod_num, ".rds"))
  fit_obj  <- readRDS(mod_path)
  
  # Plot
  Lisa_sem_graph.v2(
    fit               = fit_obj,
    title             = paste("Top Model:", dag_name, "(AIC =", round(row_match$AIC[1], 2), ")"),
    node_labels       = label_map,
    show_coefficients = TRUE,
    standardized      = TRUE,
    show_r2           = TRUE
  )
}



# 1. Helper to extract standardized structural regressions (~ only)
extract_structural_paths <- function(fit, model_name, time_period = "24yr") {
  parameterEstimates(fit, standardized = TRUE) %>%
    filter(op == "~") %>%
    mutate(
      model       = model_name,
      period      = time_period,
      path_string = paste(rhs, "->", lhs)
    ) %>%
    dplyr::select(model, period, path_string, lhs, rhs, 
                  est.std = std.all, ci.lower = ci.lower, ci.upper = ci.upper, pvalue)
}

# 2. Bind across fitted models (example combining DAG 1A, 1B, 1D, and Goal 4)
# all_paths_df <- bind_rows(
#   extract_structural_paths(fit_dag1a, "DAG 1A"),
#   extract_structural_paths(fit_dag1b, "DAG 1B"),
#   extract_structural_paths(fit_dag1d, "DAG 1D"),
#   extract_structural_paths(fit_goal4, "Goal 4 (Ocean Pred + Sea Lion)")
# )

# 3. Forest Plot
plot_comparative_forest <- function(all_paths_df) {
  ggplot(all_paths_df, aes(x = est.std, y = path_string, color = model, shape = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(xmin = ci.lower, xmax = ci.upper), 
                    position = position_dodge(width = 0.6), size = 0.6, linewidth = 0.8) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Comparison of Shared Structural Coefficients Across Models",
      subtitle = "Standardized betas with 95% confidence intervals",
      x = "Standardized Path Coefficient (std.all)",
      y = "Structural Path",
      color = "Model Topology",
      shape = "Model Topology"
    ) +
    theme_bw(base_family = "Times New Roman") +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 14)
    )
}


latent_dag1a_syntax <- '
  # ====================================================================
  # MEASUREMENT MODEL (=~)
  # ====================================================================
  # Latent Predator NCC (e.g., Hake DFA decomposed into raw indicators)
  PredNCC =~ hake_acoustic_biomass + hake_age5plus + other_pred_ind
  
  # Latent Early Growth (if measured by Condition Factor + Length Increment)
  # Growth =~ Cond1_ind1 + Cond1_ind2

  # ====================================================================
  # STRUCTURAL MODEL (~) - DAG 1A Horizontal Topology
  # ====================================================================
  Growth    ~ PreyNCC_smoothed
  Abundance ~ Growth + PredNCC
  SAR       ~ Abundance + PreyAK_DFA1 + PredAK_smoothed

  # ====================================================================
  # RESIDUAL COVARIANCES / DISTURBANCES (~~)
  # ====================================================================
  # Optional: allow non-causal correlations between upstream inputs
  PreyNCC_smoothed ~~ PredNCC
'

# Fit with standard errors via ML / MLR
latent_fit_dag1a <- sem(
  model      = latent_dag1a_syntax,
  data       = raw_df_master,
  std.lv     = TRUE,          # Set latent metric via variance = 1
  missing    = "ML"
)

summary(latent_fit_dag1a, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

top_step1_models %>% 
  dplyr::select(model_label, AIC, PreyNCCindNames, PredNCCindNames, GrowthIndNames, AbundanceIndNames, PreyAKindNames, PredAKindNames)
