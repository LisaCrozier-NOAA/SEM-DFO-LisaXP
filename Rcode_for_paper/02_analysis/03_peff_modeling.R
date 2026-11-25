# ==============================================================================
# Script: 02_analysis/03_peff_modeling.R
# Purpose: Objective 3 - Data scaling, lookup table generation, full Peff sweep,
#          export of Peff time series, and consolidated comparison table.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(gt)
})

# ------------------------------------------------------------------------------
# STEP 0: PATHS & MASTER DATA LOAD
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Outputs")
data_out_dir <- file.path(master_out, "data")
tbl_out_dir  <- file.path(master_out, "tables")

dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)

raw_csv_path <- file.path(proj_dir, "metadata", "sem_altprey_data_1998_2021_10212026.csv")
if (!file.exists(raw_csv_path)) {
  stop("Missing master raw altprey dataset at: ", raw_csv_path)
}

raw_df_master <- read.csv(raw_csv_path, stringsAsFactors = FALSE, row.names = 1)
raw_df_master <- raw_df_master[, names(raw_df_master) != "" & !is.na(names(raw_df_master))]

# Crosswalk Lookup
crosswalk_path <- file.path(proj_dir, "metadata", "master_name_crosswalk.csv")
lookup_dict <- character(0)
if (file.exists(crosswalk_path)) {
  cw <- read.csv(crosswalk_path, stringsAsFactors = FALSE)
  lookup_dict <- c(setNames(cw$LisaName, cw$dfa_cols), setNames(cw$LisaName, cw$DFAname))
}

# Baseline exclusions enforced across the entire paper
baseline_exclusions <- c("X15_sablefishAK_biomass4plus", "X15_pollockGOA_SSB_lag1")

# ------------------------------------------------------------------------------
# STEP 1: HELPER FUNCTIONS FOR SCALING & PREY POOLS
# ------------------------------------------------------------------------------
scale_custom_bounds <- function(x, a = 0.2, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

get_prey_pool <- function(guild_num, stage_suffix = "smoltyr", pred_col = "", all_cols) {
  sfx <- if (nzchar(stage_suffix)) paste0("_", stage_suffix) else "_smoltyr"
  
  m_squid     <- paste0("X04_marketsquid_GAM", sfx)
  anchovy     <- paste0("X05_anchovy_GAM", sfx)
  sardine     <- paste0("X05_DFA_abundSardine", sfx)
  herring_ncc <- paste0("X05_herring_NCC", sfx)
  eulachon    <- paste0("X05_eulachon_during_chinook", sfx)
  hake        <- paste0("X09_DFA_HakeAge5Plus", sfx)
  euphausiids <- paste0("biomassEuphShelfSum", sfx)
  
  pollock     <- paste0("X13_pollockGOA_age1and2", sfx)
  herring_ak  <- paste0("sitkaHerring_EGoA", sfx)
  capelin     <- paste0("capelin_WGoA", sfx)
  krill_ak    <- paste0("X12_egoa_krill", sfx)
  pink_ak     <- paste0("X14_pinkSalmonNorthAmerica", sfx)
  
  prey_list <- switch(
    as.character(guild_num),
    "1" = c(m_squid, anchovy, sardine, herring_ncc, eulachon),
    "2" = c(euphausiids, m_squid, anchovy, sardine, herring_ncc, eulachon, hake),
    "3" = c(pollock, herring_ak, capelin, krill_ak),
    "4" = c(pollock, herring_ak, capelin),
    "5" = c(pink_ak, pollock, herring_ak, capelin),
    character(0)
  )
  
  prey_list <- setdiff(prey_list, pred_col)
  intersect(prey_list, all_cols)
}

get_spatial_overlap_col <- function(region, stage_suffix = "smoltyr", all_cols) {
  sfx <- if (nzchar(stage_suffix)) paste0("_", stage_suffix) else "_smoltyr"
  col <- if (region == "NCC") paste0("SST_NCC_aprjun_DFObccm", sfx) else paste0("X21_sst_egoa_junjulaug", sfx)
  if (col %in% all_cols) col else NA_character_
}

# ------------------------------------------------------------------------------
# STEP 2: PIPELINE FUNCTION RUNNING COMPLETE SWEEP PER TIME-SERIES SET
# ------------------------------------------------------------------------------
run_peff_pipeline_for_years <- function(years_included = "24yr") {
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf(">>> Running Objective 3 Peff Pipeline for: %s\n", years_included))
  cat(sprintf("======================================================\n"))
  
  # 1. Filter rows & columns
  min_n <- if (years_included == "24yr") 24 else 19
  valid_cols <- names(raw_df_master)[colSums(!is.na(raw_df_master)) >= min_n]
  
  df_raw <- raw_df_master %>% 
    select(all_of(valid_cols)) %>% 
    filter(complete.cases(.))
  
  # 2. Scale dataframe: [0.0, 1.0] bounds
  df_scaled <- df_raw %>%
    mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))
  
  all_cols <- names(df_scaled)
  
  # 3. Build Lookup Table (Filtering baseline exclusions directly)
  predator_cols <- sort(all_cols[grepl("^(X08_|X09_|X10_|X11_|X15_)", all_cols)])
  predator_cols <- predator_cols[!sapply(predator_cols, function(x) any(sapply(baseline_exclusions, grepl, x)))]
  
  lookup_rows <- list()
  for (pred in predator_cols) {
    stage_suffix <- if (grepl("_adultyr$", pred)) "adultyr" else "smoltyr"
    region <- if (grepl("^(X08_|X09_|X10_)|Harbour_p", pred, ignore.case = TRUE)) "NCC" else "AK"
    
    guild_num <- case_when(
      region == "NCC" & grepl("^(X08_|Harbour_p)", pred, ignore.case = TRUE) ~ 1,
      region == "NCC" ~ 2,
      region == "AK" & grepl("(ssl|sablefish|spinyDogfish)", pred, ignore.case = TRUE) ~ 3,
      region == "AK" & grepl("(halibut|PacificCod|Arrowtooth)", pred, ignore.case = TRUE) ~ 4,
      region == "AK" & grepl("shark", pred, ignore.case = TRUE) ~ 5,
      TRUE ~ 3
    )
    
    prey_candidates <- get_prey_pool(guild_num, stage_suffix, pred, all_cols)
    spatial_col     <- get_spatial_overlap_col(region, stage_suffix, all_cols)
    all_cands       <- c(spatial_col, prey_candidates)
    all_cands       <- all_cands[!is.na(all_cands) & all_cands != ""]
    
    lookup_rows[[length(lookup_rows) + 1]] <- tibble(
      region              = region,
      guild_num           = guild_num,
      pred_data_col       = pred,
      stage_suffix        = stage_suffix,
      spatial_overlap_col = spatial_col,
      all_candidates_str  = paste(all_cands, collapse = ", ")
    )
  }
  
  guild_pred_lookup <- bind_rows(lookup_rows)
  write_csv(guild_pred_lookup, file.path(data_out_dir, paste0("obj3_predator_altprey_lookup_", years_included, ".csv")))
  
  # 4. Target SEM Specifications
  var_growth   <- "X06_DFA_IGF_mu_smoltyr"
  var_prey_ncc <- if ("X05_DFA_abundSardine_smoltyr" %in% all_cols) "X05_DFA_abundSardine_smoltyr" else "X04_marketsquid_GAM_smoltyr"
  var_cpue     <- "X07_cpue_IntSprJunHW_smoltyr"
  var_hake     <- "X09_DFA_HakeAge5Plus_smoltyr"
  var_sar      <- "X16_DFA_sar_UCsprCh_smoltyr"
  var_prey_ak  <- "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"
  
  syntax_template <- paste0(
    var_growth, " ~ ", var_prey_ncc, "\n",
    var_cpue,   " ~ ", var_hake, "\n",
    var_sar,    " ~ ", var_growth, " + ", var_cpue, " + ", var_prey_ak, " + TARGET_PRED"
  )
  
  # 5. Fit Function
  fit_variant <- function(p_vec, var_label, prey_label) {
    if (is.na(sd(p_vec)) || sd(p_vec) < 1e-6) return(NULL)
    d_fit <- df_scaled
    d_fit$TARGET_PRED <- p_vec
    
    fit <- tryCatch(
      sem(syntax_template, data = d_fit, missing = "ML", warn = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
    
    pe  <- parameterEstimates(fit, standardized = TRUE)
    rsq <- inspect(fit, "r2")
    p_est <- pe %>% filter(lhs == var_sar, rhs == "TARGET_PRED", op == "~")
    if (nrow(p_est) == 0) return(NULL)
    
    list(
      variant = var_label,
      prey    = prey_label,
      aic     = as.numeric(fitMeasures(fit, "aic")),
      r2      = as.numeric(rsq[[var_sar]]),
      b_pred  = p_est$std.all[1],
      p_pred  = p_est$pvalue[1],
      p_vector = p_vec
    )
  }
  
  # 6. Run Sweeps & Collect Peff Time Series
  pred_candidates <- guild_pred_lookup %>%
    filter(region == "AK" | (region == "NCC" & stage_suffix == "adultyr") | grepl("DFO", pred_data_col)) %>%
    pull(pred_data_col) %>%
    unique()
  pred_candidates <- intersect(pred_candidates, all_cols)
  pred_candidates <- pred_candidates[!sapply(pred_candidates, function(x) any(sapply(baseline_exclusions, grepl, x)))]
  
  summary_records <- list()
  peff_time_series_list <- list()
  
  for (pred_col in pred_candidates) {
    p_row <- guild_pred_lookup %>% filter(pred_data_col == pred_col) %>% slice(1)
    P_norm <- df_scaled[[pred_col]]
    
    is_salmon_specialist <- grepl("Killer_whale", pred_col, ignore.case = TRUE)
    spatial_col <- trimws(as.character(p_row$spatial_overlap_col[1]))
    has_temp    <- (!is.na(spatial_col) && spatial_col != "" && spatial_col %in% all_cols)
    P_temp      <- if (has_temp) P_norm * scale_custom_bounds(df_raw[[spatial_col]], a = 0.2, b = 1.0) else NULL
    
    if (is_salmon_specialist) {
      prey_cands <- character(0)
    } else {
      all_cands  <- unlist(strsplit(p_row$all_candidates_str[1], ", "))
      prey_cands <- setdiff(all_cands[!is.na(all_cands) & all_cands != ""], c(spatial_col, pred_col))
      prey_cands <- intersect(prey_cands, all_cols)
    }
    
    # Run Baseline
    base_res <- fit_variant(P_norm, "Nominal Biomass (Base)", "None")
    if (is.null(base_res)) next
    
    all_runs <- list(base_res)
    
    # Temp Only
    if (has_temp && !is.null(P_temp)) {
      t_res <- fit_variant(P_temp, "Temp Overlap Only", "None")
      if (!is.null(t_res)) all_runs[[length(all_runs) + 1]] <- t_res
    }
    
    # Alt Prey Combinations (Up to 2)
    if (length(prey_cands) > 0) {
      max_k <- min(length(prey_cands), 2)
      for (k in 1:max_k) {
        combos <- combn(prey_cands, k, simplify = FALSE)
        for (cb in combos) {
          raw_f <- if (length(cb) == 1) df_scaled[[cb[1]]] else rowSums(df_scaled[, cb, drop = FALSE], na.rm = TRUE)
          gate_alt <- 1.0 - scale_custom_bounds(raw_f, a = 0.2, b = 1.0)
          
          ap_res <- fit_variant(P_norm * gate_alt, paste0("AltPrey Buffered (", k, ")"), paste(cb, collapse = " + "))
          if (!is.null(ap_res)) all_runs[[length(all_runs) + 1]] <- ap_res
          
          if (has_temp && !is.null(P_temp)) {
            jt_res <- fit_variant(P_temp * gate_alt, paste0("Joint Temp + AltPrey (", k, ")"), paste(cb, collapse = " + "))
            if (!is.null(jt_res)) all_runs[[length(all_runs) + 1]] <- jt_res
          }
        }
      }
    }
    
    df_runs <- bind_rows(lapply(all_runs, function(x) {
      x$p_vector <- NULL
      as_tibble(x)
    })) %>%
      mutate(delta_aic = round(aic - base_res$aic, 2))
    
    best_abs <- df_runs %>% arrange(aic) %>% slice(1)
    best_neg_df <- df_runs %>% filter(b_pred < 0) %>% arrange(aic)
    best_neg <- if (nrow(best_neg_df) > 0) best_neg_df %>% slice(1) else NULL
    
    selected_model <- best_abs
    if (!is.null(best_neg)) {
      if (best_abs$b_pred >= 0 || (best_neg$aic - best_abs$aic) <= 2.0) {
        selected_model <- best_neg
      }
    }
    
    best_run_idx <- which(sapply(all_runs, function(x) x$variant == selected_model$variant && x$prey == selected_model$prey))[1]
    best_p_vec   <- all_runs[[best_run_idx]]$p_vector
    
    clean_pname <- gsub("^(X08_|X09_|X10_|X11_|X15_)", "", pred_col)
    peff_time_series_list[[paste0("Peff_", clean_pname)]] <- round(best_p_vec, 4)
    peff_time_series_list[[paste0("Pbase_", clean_pname)]] <- round(P_norm, 4)
    
    clean_prey_display <- selected_model$prey
    if (clean_prey_display != "None") {
      clean_prey_display <- gsub("X[0-9]{2}_", "", clean_prey_display)
      clean_prey_display <- gsub("(_smoltyr|_adultyr)", "", clean_prey_display)
    }
    
    summary_records[[length(summary_records) + 1]] <- tibble(
      years_set      = years_included,
      predator       = pred_col,
      base_beta      = round(base_res$b_pred, 3),
      base_p         = round(base_res$p_pred, 3),
      base_aic       = round(base_res$aic, 2),
      mod_variant    = selected_model$variant,
      mod_prey       = clean_prey_display,
      mod_beta       = round(selected_model$b_pred, 3),
      mod_p          = round(selected_model$p_pred, 3),
      mod_aic        = round(selected_model$aic, 2),
      delta_aic      = round(selected_model$aic - base_res$aic, 2),
      sign_restored  = if_else(base_res$b_pred >= 0 & selected_model$b_pred < 0, "Restored (-)", "Consistent")
    )
  }
  
  # Export Time Series Matrix
  time_vec <- if ("Year" %in% names(df_raw)) df_raw$Year else (1998:(1998 + nrow(df_raw) - 1))
  peff_ts_df <- bind_cols(tibble(Year = time_vec), as_tibble(peff_time_series_list))
  write_csv(peff_ts_df, file.path(data_out_dir, paste0("obj3_data_peff_timeseries_", years_included, ".csv")))
  
  sweep_summary_df <- bind_rows(summary_records)
  write_csv(sweep_summary_df, file.path(data_out_dir, paste0("obj3_data_peff_sweep_all_models_", years_included, ".csv")))
  
  cat(sprintf("[%s] Successfully generated lookup, sweep results, and Peff time series!\n", years_included))
  return(sweep_summary_df)
}

# Run both 24yr and 19yr sweeps
res_24 <- run_peff_pipeline_for_years("24yr")
res_19 <- run_peff_pipeline_for_years("19yr")

# ------------------------------------------------------------------------------
# STEP 3: CONSTRUCT UNIFIED PUBLICATION TABLE (TIGHTENED FOR PORTRAIT PDF)
# ------------------------------------------------------------------------------
cat("\nGenerating Publication Table for Objective 3...\n")

# Helper function to truncate and shorten verbose predator names
shorten_pred_label <- function(x) {
  x <- gsub("^(X08_|X09_|X10_|X11_|X15_)", "", x)
  x <- gsub("DFA_", "", x)
  x <- gsub("_predAK", "", x)
  x <- gsub("Biomass", "", x)
  x <- ifelse(nchar(x) > 24, paste0(substr(x, 1, 22), ".."), x)
  return(x)
}

# Combine and enforce Panel ordering (24yr long series on top)
combined_summary <- bind_rows(res_24, res_19) %>%
  mutate(
    raw_mod_variant = mod_variant,  # Keep unshortened mechanism for conditional styling logic
    Panel = factor(
      if_else(years_set == "24yr", 
              "Panel A: Long / 24-Year Series (1998–2021)", 
              "Panel B: Short / 19-Year Series (2003–2021)"),
      levels = c(
        "Panel A: Long / 24-Year Series (1998–2021)", 
        "Panel B: Short / 19-Year Series (2003–2021)"
      )
    ),
    Predator_Clean = case_when(
      predator %in% names(lookup_dict) ~ shorten_pred_label(lookup_dict[predator]),
      TRUE ~ shorten_pred_label(predator)
    ),
    # Compact mechanism labels for presentation
    mod_variant = gsub("Joint Temp \\+ AltPrey", "Joint Temp+Prey", mod_variant),
    mod_variant = gsub("AltPrey Buffered", "AltPrey Only", mod_variant),
    mod_variant = gsub("Temp Overlap Only", "Temp Only", mod_variant),
    mod_variant = gsub("Nominal Biomass \\(Base\\)", "Base (Nominal)", mod_variant),
    
    Base_Display = paste0(sprintf("%.2f", base_beta), " (", sprintf("%.2f", base_p), ")"),
    Mod_Display  = paste0(sprintf("%.2f", mod_beta), " (", sprintf("%.2f", mod_p), ")")
  ) %>%
  arrange(Panel, mod_aic)

table_obj3_gt <- combined_summary %>%
  select(
    Panel, Predator_Clean, 
    Base_Display, base_aic, 
    mod_variant, mod_prey, 
    Mod_Display, mod_aic, delta_aic,
    # Hidden columns for evaluation in styling expressions
    raw_mod_variant, base_beta, base_p, mod_beta, mod_p
  ) %>%
  gt(groupname_col = "Panel") %>%
  tab_header(
    title = md("**Table 3. Effective Predation (P<sub>eff</sub>) vs. Nominal Biomass**"),
    subtitle = "Evaluating whether temperature overlap and alternative forage buffering restore expected top-down control."
  ) %>%
  tab_spanner(label = "Nominal (Base)", columns = c(Base_Display, base_aic)) %>%
  tab_spanner(label = "Effective Predation (Peff)", columns = c(mod_variant, mod_prey, Mod_Display, mod_aic, delta_aic)) %>%
  cols_label(
    Predator_Clean = md("**Predator**"),
    Base_Display   = md("**Base &beta; (p)**"),
    base_aic       = md("**Base AIC**"),
    mod_variant    = md("**Mechanism**"),
    mod_prey       = md("**Alt. Prey Pool**"),
    Mod_Display    = md("**P<sub>eff</sub> &beta; (p)**"),
    mod_aic        = md("**P<sub>eff</sub> AIC**"),
    delta_aic      = md("**&Delta;AIC**")
  ) %>%
  cols_hide(columns = c(raw_mod_variant, base_beta, base_p, mod_beta, mod_p)) %>%
  cols_align(align = "left", columns = c(Predator_Clean, mod_variant, mod_prey)) %>%
  cols_align(align = "center", columns = c(Base_Display, base_aic, Mod_Display, mod_aic, delta_aic)) %>%
  cols_width(
    Predator_Clean ~ px(125),
    Base_Display   ~ px(75),
    base_aic       ~ px(55),
    mod_variant    ~ px(100),
    mod_prey       ~ px(135),
    Mod_Display    ~ px(75),
    mod_aic        ~ px(55),
    delta_aic      ~ px(48)
  ) %>%
  # ----------------------------------------------------------------------------
# CONDITIONAL CELL HIGHLIGHTING
# ----------------------------------------------------------------------------
# 1. Base beta green: Nominal Biomass selected AND base p <= 0.05 AND base beta < 0
tab_style(
  style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#1B5E20", weight = "bold")),
  locations = cells_body(
    columns = Base_Display,
    rows = (raw_mod_variant == "Nominal Biomass (Base)" | mod_variant == "Base (Nominal)") & 
      base_p <= 0.05 & 
      base_beta < 0
  )
) %>%
  # 2. Peff beta green: Mechanism is NOT Nominal Base AND mod p <= 0.05 AND mod beta < 0
  tab_style(
    style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#1B5E20", weight = "bold")),
    locations = cells_body(
      columns = Mod_Display,
      rows = (raw_mod_variant != "Nominal Biomass (Base)" & mod_variant != "Base (Nominal)") & 
        mod_p <= 0.05 & 
        mod_beta < 0
    )
  ) %>%
  # 3. Delta AIC green: delta_aic is negative
  tab_style(
    style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#1B5E20", weight = "bold")),
    locations = cells_body(
      columns = delta_aic,
      rows = delta_aic < 0
    )
  ) %>%
  tab_source_note(
    source_note = md("P<sub>eff</sub> attenuates predation via temperature gating (P &times; SST) and/or alternative prey buffering (P &times; [1 &minus; AltPrey]). Shaded green cells highlight statistically significant negative top-down impacts (&beta; < 0, p &le; 0.05) or model improvements (&Delta;AIC < 0).")
  ) %>%
  tab_options(
    table.font.size            = px(8.5),
    heading.title.font.size    = px(10.5),
    heading.subtitle.font.size = px(8.5),
    column_labels.font.size    = px(8.5),
    column_labels.font.weight  = "bold",
    row_group.font.size        = px(9.0),
    row_group.font.weight      = "bold",
    row_group.padding          = px(2),
    data_row.padding           = px(1.5),
    source_notes.font.size     = px(7.5)
  )

# Output PDF and PNG with compact viewport dimensions
gtsave(table_obj3_gt, file.path(tbl_out_dir, "obj3_table1_peff_model_comparisons.pdf"))
gtsave(table_obj3_gt, file.path(tbl_out_dir, "obj3_table1_peff_model_comparisons.png"), vwidth = 680, vheight = 920)

cat("\nTightened Table 3 successfully exported to Outputs/tables/!\n")



#Quarto prep-------
# ------------------------------------------------------------------------------
# STEP 4: EXPORT CLEAN DISPLAY CSVs FOR QUARTO (24yr Main + 19yr Supplement)
# ------------------------------------------------------------------------------
export_peff_display_csv <- function(sweep_df, out_csv_path) {
  sweep_df %>%
    mutate(
      raw_mod_variant = mod_variant,
      Predator_Clean = case_when(
        predator %in% names(lookup_dict) ~ shorten_pred_label(lookup_dict[predator]),
        TRUE ~ shorten_pred_label(predator)
      ),
      # Replace underscores with spaces so LaTeX never needs escaping
      Predator_Clean = gsub("_", " ", Predator_Clean),
      mod_prey_clean = gsub("_", " ", mod_prey),
      
      # Compact mechanism labels
      mod_variant = gsub("Joint Temp \\+ AltPrey", "Joint Temp+Prey", mod_variant),
      mod_variant = gsub("AltPrey Buffered", "AltPrey Only", mod_variant),
      mod_variant = gsub("Temp Overlap Only", "Temp Only", mod_variant),
      mod_variant = gsub("Nominal Biomass \\(Base\\)", "Base (Nominal)", mod_variant),
      
      base_str = paste0(sprintf("%.2f", base_beta), " (", sprintf("%.2f", base_p), ")"),
      mod_str  = paste0(sprintf("%.2f", mod_beta), " (", sprintf("%.2f", mod_p), ")"),
      daic_str = sprintf("%.2f", delta_aic),
      
      # Green = significant negative top-down or delta_aic < 0; Red = positive & p < 0.20
      Base_Display = case_when(
        (raw_mod_variant == "Nominal Biomass (Base)" | mod_variant == "Base (Nominal)") &
          base_p <= 0.05 & base_beta < 0 ~ paste0("\\textcolor[HTML]{1B5E20}{\\textbf{", base_str, "}}"),
        base_beta > 0 & base_p < 0.20   ~ paste0("\\textcolor[HTML]{B71C1C}{\\textbf{", base_str, "}}"),
        TRUE                            ~ base_str
      ),
      Mod_Display = if_else(
        (raw_mod_variant != "Nominal Biomass (Base)" & mod_variant != "Base (Nominal)") &
          mod_p <= 0.05 & mod_beta < 0,
        paste0("\\textcolor[HTML]{1B5E20}{\\textbf{", mod_str, "}}"),
        mod_str
      ),
      Delta_Display = if_else(
        delta_aic < 0,
        paste0("\\textcolor[HTML]{1B5E20}{\\textbf{", daic_str, "}}"),
        daic_str
      ),
      Base_AIC_Display = sprintf("%.2f", base_aic),
      Mod_AIC_Display  = sprintf("%.2f", mod_aic)
    ) %>%
    arrange(mod_aic) %>%
    select(
      Predator_Clean, Base_Display, Base_AIC_Display,
      mod_variant, mod_prey_clean, Mod_Display, Mod_AIC_Display, Delta_Display
    ) %>%
    write_csv(out_csv_path)
}

export_peff_display_csv(res_24, file.path(data_out_dir, "obj3_table_display_24yr.csv"))
export_peff_display_csv(res_19, file.path(data_out_dir, "obj3_table_display_19yr.csv"))