# ==============================================================================
# Script: 02_analysis/05a_fit_latent_models_DAG1B.R
# Purpose: Fit Latent DAG 1B models across all 3 stages for 5 top candidates:
#          Stage 1: Full Measurement Model
#          Stage 2: Pruned Measurement Model (p < 0.05 loadings)
#          Stage 3: Pruned Structural Regressions (p < 0.05 backward elimination)
# Outputs:
#   - Outputs/data/obj4_latent_model_objects.RData
#   - Outputs/data/obj4_data_latent_model_parameters.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: DIRECTORIES & DATA SETUP
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Outputs")
data_out_dir <- file.path(master_out, "data")
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

all_data_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/allData.rds"
cw_path       <- file.path(proj_dir, "metadata", "master_name_crosswalk.csv")
alt_data_path <- file.path(proj_dir, "metadata", "sem_altprey_data_1998_2021_10212026.csv")
lookup_path   <- file.path(data_out_dir, "obj3_predator_altprey_lookup_24yr.csv")
save_rdata    <- file.path(data_out_dir, "obj4_latent_model_objects.RData")

if (!file.exists(all_data_path)) stop("Missing allData.rds at: ", all_data_path)
if (!file.exists(alt_data_path)) stop("Missing altprey raw data at: ", alt_data_path)

all_data   <- readRDS(all_data_path)
cw_df      <- read.csv(cw_path, stringsAsFactors = FALSE)
df_24_raw  <- read.csv(alt_data_path, stringsAsFactors = FALSE, row.names = 1) %>% filter(complete.cases(.))
lookup_df  <- read.csv(lookup_path, stringsAsFactors = FALSE)

scale_custom_bounds <- function(x, a = 0.2, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

df_base <- df_24_raw %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

yr_col <- grep("year", names(df_base), value = TRUE, ignore.case = TRUE)
df_base$year <- if (length(yr_col) > 0) df_base[[yr_col[1]]] else 1998:2021

var_growth  <- "X06_DFA_IGF_mu_smoltyr"
var_cpue    <- "X07_cpue_IntSprJunHW_smoltyr"
var_sar     <- "X16_DFA_sar_UCsprCh_smoltyr"
var_prey_ak <- "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"

# ------------------------------------------------------------------------------
# STEP 1: PARSE & EXTRACT 100% COMPLETE MEASUREMENT CONSTITUENTS
# ------------------------------------------------------------------------------
parse_cw_tokens <- function(node_key, cw) {
  row_match <- cw %>% filter(LisaName == node_key | DFAname == node_key | dfa_cols == node_key)
  raw_tokens <- trimws(unlist(strsplit(row_match$rankedIndicators[1], "\\s*-\\s*|\\s*,\\s*")))
  raw_tokens <- raw_tokens[raw_tokens != "" & !is.na(raw_tokens)]
  unique(gsub("_[0-9]{4}$", "", raw_tokens))
}

raw_x05_tokens <- parse_cw_tokens("X05_DFA_abundSardine", cw_df)
raw_x09_tokens <- parse_cw_tokens("X09_DFA_HakeAge5Plus", cw_df)

extract_complete_constituents <- function(all_df, token_list) {
  df_wide <- all_df %>%
    mutate(
      year            = year(ymd(date)),
      shortName_clean = gsub("_[0-9]{4}$", "", shortName)
    ) %>%
    filter(
      year >= 1998, year <= 2021,
      yearType == "smolt",
      shortName_clean %in% token_list
    ) %>%
    select(year, shortName = shortName_clean, finalVal) %>%
    distinct(year, shortName, .keep_all = TRUE) %>%
    pivot_wider(names_from = shortName, values_from = finalVal) %>%
    arrange(year)
  
  cand_cols     <- setdiff(names(df_wide), "year")
  col_nas       <- sapply(cand_cols, function(col) sum(is.na(df_wide[[col]])))
  complete_cols <- names(col_nas[col_nas == 0])
  df_wide %>% select(year, all_of(complete_cols))
}

df_x05_comp <- extract_complete_constituents(all_data, raw_x05_tokens)
df_x09_comp <- extract_complete_constituents(all_data, raw_x09_tokens)

df_x05_std <- df_x05_comp %>% mutate(across(-year, ~ scale_custom_bounds(.x, 0.0, 1.0)))
df_x09_std <- df_x09_comp %>% mutate(across(-year, ~ scale_custom_bounds(.x, 0.0, 1.0)))

avail_x05 <- setdiff(names(df_x05_std), "year")
avail_x09 <- setdiff(names(df_x09_std), "year")

df_master <- df_base %>%
  left_join(df_x05_std, by = "year") %>%
  left_join(df_x09_std, by = "year")

# ------------------------------------------------------------------------------
# STEP 2: TOP 5 CANDIDATES (INCLUDING SLEEPER SHARK)
# ------------------------------------------------------------------------------
selected_models <- tribble(
  ~model_id, ~predator,                                  ~variant,                   ~alt_prey,
  1,         "X11_ssl_seak_pup_pred_smoltyr",            "AltPrey Only (2)",         "sitkaHerring_EGoA_smoltyr + X12_egoa_krill_smoltyr",
  2,         "X15_sablefishAK_biomass3plus_smoltyr",      "Joint Temp + AltPrey (2)", "pollockGOA_age1and2_smoltyr + sitkaHerring_EGoA_smoltyr",
  3,         "X10_Killer_whales_NR_BC_adultyr",           "Temp Only",                "None",
  4,         "X15_DFA_sleeperSharkBSAI_predAK_smoltyr",   "Joint Temp + AltPrey (1)", "capelin_WGoA_smoltyr",
  5,         "X15_pollockGOA_SSB_smoltyr",                "Joint Temp + AltPrey (2)", "sitkaHerring_EGoA_smoltyr + X12_egoa_krill_smoltyr"
)
# ------------------------------------------------------------------------------
# STEP 3: SYNTAX BUILDER & PRUNING ENGINES
# ------------------------------------------------------------------------------
build_full_syntax <- function(x05_vec, x09_vec, reg_growth = NULL, reg_cpue = NULL, reg_sar = NULL) {
  meas_blocks <- character(0)
  
  if (length(x05_vec) >= 2) {
    meas_blocks <- c(meas_blocks, paste0("Latent_05_PreyNCC =~ ", paste(x05_vec, collapse = " + ")))
    default_rhs_growth <- "Latent_05_PreyNCC"
  } else if (length(x05_vec) == 1) {
    default_rhs_growth <- x05_vec[1]
  } else {
    default_rhs_growth <- "X05_DFA_abundSardine_smoltyr"
  }
  
  if (length(x09_vec) >= 2) {
    meas_blocks <- c(meas_blocks, paste0("Latent_09_PredNCC =~ ", paste(x09_vec, collapse = " + ")))
    default_rhs_cpue <- "Latent_09_PredNCC"
  } else if (length(x09_vec) == 1) {
    default_rhs_cpue <- x09_vec[1]
  } else {
    default_rhs_cpue <- "X09_DFA_HakeAge5Plus_smoltyr"
  }
  
  rhs_growth <- if (is.null(reg_growth)) default_rhs_growth else reg_growth
  rhs_cpue   <- if (is.null(reg_cpue))   default_rhs_cpue   else reg_cpue
  rhs_sar    <- if (is.null(reg_sar))    c(var_growth, var_cpue, var_prey_ak, "P_eff") else reg_sar
  
  struct_lines <- character(0)
  if (length(rhs_growth) > 0) struct_lines <- c(struct_lines, paste0(var_growth, " ~ ", paste(rhs_growth, collapse = " + ")))
  if (length(rhs_cpue) > 0)   struct_lines <- c(struct_lines, paste0(var_cpue,   " ~ ", paste(rhs_cpue,   collapse = " + ")))
  if (length(rhs_sar) > 0)    struct_lines <- c(struct_lines, paste0(var_sar,    " ~ ", paste(rhs_sar,    collapse = " + ")))
  
  paste(c(meas_blocks, struct_lines), collapse = "\n")
}

prune_measurement_loadings <- function(df_data, init_x05, init_x09, p_thresh = 0.05) {
  cur_x05 <- init_x05
  cur_x09 <- init_x09
  repeat {
    syntax <- build_full_syntax(cur_x05, cur_x09)
    fit    <- tryCatch(sem(syntax, data = df_data, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
    if (is.null(fit) || !lavInspect(fit, "converged")) break
    
    pe <- parameterEstimates(fit)
    ns_loadings <- pe %>% filter(op == "=~" & !is.na(pvalue) & pvalue >= p_thresh) %>% arrange(desc(pvalue))
    if (nrow(ns_loadings) == 0) break
    
    worst_ind <- ns_loadings$rhs[1]
    worst_lat <- ns_loadings$lhs[1]
    if (worst_lat == "Latent_05_PreyNCC") cur_x05 <- setdiff(cur_x05, worst_ind)
    if (worst_lat == "Latent_09_PredNCC") cur_x09 <- setdiff(cur_x09, worst_ind)
  }
  final_syntax <- build_full_syntax(cur_x05, cur_x09)
  final_fit    <- sem(final_syntax, data = df_data, std.lv = TRUE, missing = "ML", warn = FALSE)
  list(fit = final_fit, syntax = final_syntax, x05 = cur_x05, x09 = cur_x09)
}

prune_structural_regressions <- function(df_data, x05_vec, x09_vec, p_thresh = 0.05) {
  rhs_growth <- if (length(x05_vec) >= 2) "Latent_05_PreyNCC" else x05_vec[1]
  rhs_cpue   <- if (length(x09_vec) >= 2) "Latent_09_PredNCC" else x09_vec[1]
  rhs_sar    <- c(var_growth, var_cpue, var_prey_ak, "P_eff")
  
  iteration <- 1
  repeat {
    syntax <- build_full_syntax(x05_vec, x09_vec, reg_growth = rhs_growth, reg_cpue = rhs_cpue, reg_sar = rhs_sar)
    fit    <- tryCatch(sem(syntax, data = df_data, std.lv = TRUE, missing = "ML", warn = FALSE), error = function(e) NULL)
    if (is.null(fit) || !lavInspect(fit, "converged")) break
    
    pe <- parameterEstimates(fit)
    ns_regs <- pe %>% filter(op == "~" & !is.na(pvalue) & pvalue >= p_thresh) %>% arrange(desc(pvalue))
    if (nrow(ns_regs) == 0) break
    
    worst_lhs <- ns_regs$lhs[1]
    worst_rhs <- ns_regs$rhs[1]
    if (worst_lhs == var_sar)    rhs_sar    <- setdiff(rhs_sar, worst_rhs)
    if (worst_lhs == var_growth) rhs_growth <- setdiff(rhs_growth, worst_rhs)
    if (worst_lhs == var_cpue)   rhs_cpue   <- setdiff(rhs_cpue, worst_rhs)
    iteration <- iteration + 1
  }
  final_syntax <- build_full_syntax(x05_vec, x09_vec, reg_growth = rhs_growth, reg_cpue = rhs_cpue, reg_sar = rhs_sar)
  final_fit    <- sem(final_syntax, data = df_data, std.lv = TRUE, missing = "ML", warn = FALSE)
  list(fit = final_fit, syntax = final_syntax, rhs_sar = rhs_sar, rhs_growth = rhs_growth, rhs_cpue = rhs_cpue)
}

# ------------------------------------------------------------------------------
# STEP 4: EXECUTE ALL 3 STAGES ACROSS 5 MODELS
# ------------------------------------------------------------------------------
dag1b_candidate_models <- list()
param_summary_list     <- list()

for (i in seq_len(nrow(selected_models))) {
  row_m    <- selected_models[i, ]
  mod_id   <- row_m$model_id
  pred_col <- row_m$predator
  variant  <- row_m$variant
  prey_str <- row_m$alt_prey
  m_key    <- paste0("Model_", mod_id)
  
  cat(sprintf("\n>>> Fitting %s: %s [%s]...\n", m_key, pred_col, variant))
  
  d_fit <- df_master
  # Validate predator column exists
  if (!pred_col %in% names(d_fit)) {
    # Check if adding or removing 'X' prefix resolves it
    alt_name <- if (grepl("^X", pred_col)) sub("^X", "", pred_col) else paste0("X", pred_col)
    if (alt_name %in% names(d_fit)) {
      pred_col <- alt_name
    } else {
      stop(sprintf("Column '%s' (or '%s') not found in df_master!", pred_col, alt_name))
    }
  }
  
  p_vec <- d_fit[[pred_col]]

  # Temperature overlap
  l_row       <- lookup_df %>% filter(pred_data_col == pred_col) %>% slice(1)
  spatial_col <- trimws(as.character(l_row$spatial_overlap_col[1]))
  has_temp    <- (grepl("Temp|Joint", variant) && !is.na(spatial_col) && spatial_col %in% names(df_24_raw))
  if (has_temp) p_vec <- p_vec * scale_custom_bounds(df_24_raw[[spatial_col]], a = 0.2, b = 1.0)
  
  # Alternate prey gate
  if (prey_str != "None" && !is.na(prey_str)) {
    cands       <- trimws(unlist(strsplit(prey_str, "\\+")))
    valid_cands <- intersect(cands, names(d_fit))
    if (length(valid_cands) > 0) {
      raw_forage <- if (length(valid_cands) == 1) d_fit[[valid_cands[1]]] else rowSums(d_fit[, valid_cands, drop = FALSE], na.rm = TRUE)
      gate_alt   <- 1.0 - scale_custom_bounds(raw_forage, a = 0.2, b = 1.0)
      p_vec      <- p_vec * gate_alt
    }
  }
  d_fit$P_eff <- p_vec
  
  # Stage 1: Full Measurement
  syntax_s1 <- build_full_syntax(avail_x05, avail_x09)
  fit_s1    <- sem(syntax_s1, data = d_fit, std.lv = TRUE, missing = "ML", warn = FALSE)
  
  # Stage 2: Pruned Measurement
  res_s2 <- prune_measurement_loadings(d_fit, avail_x05, avail_x09, p_thresh = 0.05)
  
  # Stage 3: Structural Pruned
  res_s3 <- prune_structural_regressions(d_fit, res_s2$x05, res_s2$x09, p_thresh = 0.05)
  
  dag1b_candidate_models[[m_key]] <- list(
    metadata                       = as.list(row_m),
    Stage1_Full_Measurement        = list(fit = fit_s1, syntax = syntax_s1),
    Stage2_Pruned_Measurement      = list(fit = res_s2$fit, syntax = res_s2$syntax, x05 = res_s2$x05, x09 = res_s2$x09),
    Stage3_Final_Pruned_Regression = list(fit = res_s3$fit, syntax = res_s3$syntax, sar_paths = res_s3$rhs_sar)
  )
  
  pe_s3 <- parameterEstimates(res_s3$fit, standardized = TRUE) %>%
    filter(op %in% c("=~", "~")) %>%
    mutate(model_key = m_key, predator = pred_col, variant = variant)
  param_summary_list[[i]] <- pe_s3
}

# ------------------------------------------------------------------------------
# STEP 5: SAVE ARTIFACTS
# ------------------------------------------------------------------------------
save(dag1b_candidate_models, df_master, file = save_rdata)
write_csv(bind_rows(param_summary_list), file.path(data_out_dir, "obj4_data_latent_model_parameters.csv"))

cat(sprintf("\nSuccessfully saved full latent model objects to:\n -> %s\n", save_rdata))