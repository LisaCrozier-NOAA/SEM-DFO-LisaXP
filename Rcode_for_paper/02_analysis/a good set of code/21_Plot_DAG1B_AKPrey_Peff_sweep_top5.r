

#21_Plot_DAG1B_AKPrey_Peff_sweep_top5.r

#Plot best Peff
library(tidyverse)
library(lavaan)
library(tidySEM)
library(ggplot2)

# ==============================================================================
# 0. DIRECTORY & DATA SETUP
# ==============================================================================
proj_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
data_path    <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
lookup_path  <- file.path(proj_dir, "Routput_for_paper/data/guild_predator_lookup_table_24yr.csv")
cw_path      <- file.path(proj_dir, "metadata/crosswalk_names.csv")

raw_df_master <- read.csv(data_path, stringsAsFactors = FALSE)
raw_df_master <- raw_df_master[, names(raw_df_master) != "" & !is.na(names(raw_df_master))]
df_24         <- raw_df_master %>% filter(complete.cases(.))

lookup_df <- read_csv(lookup_path, show_col_types = FALSE)
cw_df     <- read.csv(cw_path, stringsAsFactors = FALSE)

scale_custom_bounds <- function(x, a = 0.2, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

# Base min-max scaled dataframe
df_base <- df_24 %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

# Fixed structural anchors
var_growth   <- "X06_DFA_IGF_mu_smoltyr"
var_prey_ncc <- "X05_DFA_abundSardine_smoltyr"
var_cpue     <- "X07_cpue_IntSprJunHW_smoltyr"
var_hake     <- "X09_DFA_HakeAge5Plus_smoltyr"
var_sar      <- "X16_DFA_sar_UCsprCh_smoltyr"
var_prey_ak  <- "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"

# ==============================================================================
# 1. DEFINE THE 5 SELECTED CANDIDATE MODELS
# ==============================================================================
selected_models <- tribble(
  ~model_id, ~predator,                                 ~variant,                    ~alt_prey,
  1,         "X15_pollockGOA_SSB_lag1_smoltyr",         "AltPrey Only (2)",          "capelin_WGoA_smoltyr + X12_egoa_krill_smoltyr",
  3,         "X11_ssl_seak_pup_pred_smoltyr",           "AltPrey Only (2)",          "capelin_WGoA_smoltyr + X12_egoa_krill_smoltyr",
  4,         "X15_sablefishAK_biomass3plus_smoltyr",    "Joint Temp + AltPrey (2)",  "sitkaHerring_EGoA_smoltyr + capelin_WGoA_smoltyr",
  5,         "X15_DFA_sleeperSharkBSAI_predAK_smoltyr", "Raw Baseline",              "None",
  7,         "X10_Killer_whales_NR_BC_adultyr",         "Temp Only",                 "None"
)

# ==============================================================================
# 2. LABEL AND LAYOUT HELPERS
# ==============================================================================
get_clean_name <- function(raw_var) {
  match_row <- cw_df %>% filter(DFAname == raw_var | dfa_cols == raw_var | LisaName == raw_var)
  clean_name <- if (nrow(match_row) > 0 && !is.na(match_row$LisaName[1]) && match_row$LisaName[1] != "") {
    match_row$LisaName[1]
  } else {
    raw_var
  }
  clean_name <- gsub("_smoltyr|_adultyr", "", clean_name)
  clean_name <- gsub("^X[0-9]{2}_", "", clean_name)
  clean_name
}

layout_dag1b <- get_layout(
  "PreyNCC", "", "",          "", "",       "", "PreyAK",
  "",        "", "Growth",    "", "",       "", "",
  "",        "", "",          "", "logSAR", "", "",
  "",        "", "Abundance", "", "",       "", "",
  "PredNCC", "", "",          "", "",       "", "P_eff",
  rows = 5
)

# ==============================================================================
# 3. FIT AND PLOT LOOP
# ==============================================================================
plot_list <- list()

for (i in seq_len(nrow(selected_models))) {
  row_m    <- selected_models[i, ]
  mod_idx  <- row_m$model_id
  pred_col <- row_m$predator
  variant  <- row_m$variant
  prey_str <- row_m$alt_prey
  
  # A. Isolate data frame
  d_fit <- df_base
  p_vec <- d_fit[[pred_col]]
  
  # B. Lookup and apply spatial temperature overlap
  l_row       <- lookup_df %>% filter(pred_data_col == pred_col) %>% slice(1)
  spatial_col <- trimws(as.character(l_row$spatial_overlap_col[1]))
  has_temp    <- (grepl("Temp|Joint", variant) && !is.na(spatial_col) && spatial_col %in% names(df_24))
  
  if (has_temp) {
    p_vec <- p_vec * scale_custom_bounds(df_24[[spatial_col]], a = 0.2, b = 1.0)
  }
  
  # C. Apply alternate prey attenuation gate
  if (prey_str != "None" && !is.na(prey_str)) {
    cands       <- trimws(unlist(strsplit(prey_str, "\\+")))
    valid_cands <- intersect(cands, names(d_fit))
    if (length(valid_cands) > 0) {
      raw_forage <- if (length(valid_cands) == 1) d_fit[[valid_cands[1]]] else rowSums(d_fit[, valid_cands, drop = FALSE], na.rm = TRUE)
      gate_alt   <- 1.0 - 1.0 * scale_custom_bounds(raw_forage, a = 0.2, b = 1.0)
      p_vec      <- p_vec * gate_alt
    }
  }
  
  # Attach P_eff into the isolated dataset
  d_fit$P_eff <- p_vec
  
  # Map columns to structural names for Lisa_sem_graph.v3
  d_fit$Growth    <- d_fit[[var_growth]]
  d_fit$PreyNCC   <- d_fit[[var_prey_ncc]]
  d_fit$Abundance <- d_fit[[var_cpue]]
  d_fit$PredNCC   <- d_fit[[var_hake]]
  d_fit$logSAR    <- d_fit[[var_sar]]
  d_fit$PreyAK    <- d_fit[[var_prey_ak]]
  
  syntax_dag1b <- '
    Growth    ~ PreyNCC
    Abundance ~ PredNCC
    logSAR    ~ Growth + Abundance + PreyAK + P_eff
  '
  
  fit_obj <- sem(syntax_dag1b, data = d_fit, missing = "ML", warn = FALSE)
  
  # Print parameter checks to confirm exact match with sweep table
  aic_val <- round(as.numeric(fitMeasures(fit_obj, "aic")), 2)
  pe      <- parameterEstimates(fit_obj, standardized = TRUE)
  
  p_peff  <- pe %>% filter(lhs == "logSAR", rhs == "P_eff")
  p_grw   <- pe %>% filter(lhs == "logSAR", rhs == "Growth")
  p_pak   <- pe %>% filter(lhs == "logSAR", rhs == "PreyAK")
  
  cat(sprintf("\n--- Model #%d: %s [%s] ---\n", mod_idx, get_clean_name(pred_col), variant))
  cat(sprintf("AIC: %.2f | P_eff beta: %.3f (p = %.4f) | Growth beta: %.3f (p = %.4f) | PreyAK beta: %.3f (p = %.4f)\n",
              aic_val, p_peff$std.all[1], p_peff$pvalue[1], p_grw$std.all[1], p_grw$pvalue[1], p_pak$std.all[1], p_pak$pvalue[1]))
  
  # D. Node label formatting
  pred_clean <- get_clean_name(pred_col)
  
  temp_tag <- if (has_temp) paste0("\nTemp: ", get_clean_name(spatial_col)) else ""
  alt_tag  <- if (prey_str != "None") {
    p_items <- trimws(unlist(strsplit(prey_str, "\\+")))
    paste0("\nAlt: ", paste(sapply(p_items, get_clean_name), collapse = " + "))
  } else {
    ""
  }
  
  labels_map <- c(
    "PreyNCC"   = paste0("Prey (NCC)\n[", get_clean_name(var_prey_ncc), "]"),
    "PredNCC"   = paste0("Predators (NCC)\n[", get_clean_name(var_hake), "]"),
    "Growth"    = paste0("Salmon Growth\n[", get_clean_name(var_growth), "]"),
    "Abundance" = paste0("Salmon Abund.\n[", get_clean_name(var_cpue), "]"),
    "PreyAK"    = paste0("Prey (AK)\n[", get_clean_name(var_prey_ak), "]"),
    "logSAR"    = "Salmon SAR\n[UCsprCh]",
    "P_eff"     = paste0("P_eff: ", pred_clean, temp_tag, alt_tag)
  )
  
  # E. Generate plot
  plot_title <- sprintf("Candidate Model #%d: %s (%s, AIC = %.2f)", 
                        mod_idx, pred_clean, variant, aic_val)
  
  p_graph <- Lisa_sem_graph.v3(
    fit               = fit_obj,
    title             = plot_title,
    layout            = layout_dag1b,
    node_labels       = labels_map,
    r2_nodes          = c("Growth", "Abundance", "logSAR"),
    p_threshold       = 0.05,
    nonsig_linetype   = "dashed",
    font_family       = "serif",
    rect_width        = 2.8,
    rect_height       = 0.95,
    ellipses_width    = 2.8,
    ellipses_height   = 0.95,
    text_size         = 3.2,
    variance_diameter = 0.6,
    curvature         = 10,
    width             = 15,
    height            = 9,
    display           = TRUE,
    save              = FALSE
  )
  
  plot_list[[as.character(mod_idx)]] <- p_graph
}