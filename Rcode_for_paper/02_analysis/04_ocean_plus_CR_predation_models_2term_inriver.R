# ==============================================================================
# Script: 02_analysis/04_ocean_plus_CR_predation_models_2term_inriver.R
# Purpose: Objective 4 - Life-Stage Arena Bottlenecks
#          Single: SAR ~ Growth + CPUE
#          Dual:   SAR ~ Growth + CPUE + PreyAK + P_ocean_eff
#          Triple: SAR ~ Growth + CPUE + PreyAK + P_ocean_eff + 
#                        InRiver_SeaLions + F_known_removals
# Outputs:
#   - Outputs/data/obj4_data_ocean_plus_CR_predation_2term_inriver.csv
#   - Outputs/tables/obj4_table1_ocean_plus_CR_predation_2term_inriver.[pdf|png]
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(gt)
})

# ------------------------------------------------------------------------------
# STEP 0: SETUP PATHS & DATA
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Outputs")
data_out_dir <- file.path(master_out, "data")
tbl_out_dir  <- file.path(master_out, "tables")

dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_out_dir, showWarnings = FALSE, recursive = TRUE)

meta_file <- file.path(proj_dir, "metadata", "sem_data_with_harvest_bird_pinn_predation.csv")
if (!file.exists(meta_file)) {
  stop("Missing merged dataset at: ", meta_file)
}

raw_df_master <- read.csv(meta_file, stringsAsFactors = FALSE)
raw_df_master <- raw_df_master[, names(raw_df_master) != "" & !is.na(names(raw_df_master))]

# Standardize continuous columns to Z-scores (mean = 0, sd = 1)
df_base <- raw_df_master %>%
  mutate(across(where(is.numeric) & !matches("Year|year|smoltyr|return_yr|id"), ~ as.numeric(scale(.x))))

# Structural anchors
target_growth        <- "X06_DFA_IGF_mu_smoltyr"
target_prey_ncc      <- if ("X05_DFA_abundSardine_smoltyr" %in% names(df_base)) "X05_DFA_abundSardine_smoltyr" else "X04_marketsquid_GAM_smoltyr"
target_cpue          <- "X07_cpue_IntSprJunHW_smoltyr"
hake_col             <- "X09_DFA_HakeAge5Plus_smoltyr"
target_prey_ak       <- "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"
target_sar           <- "X16_DFA_sar_UCsprCh_smoltyr"

# In-River 2-Term Predictors
inriver_sealions_col <- "X10_DFA_AllSeaLionsEMB_adultyr"
inriver_mort_col     <- "z_f_inriver_total"  # bird + tailrace pinniped + Z1-5 harvest

# Load Objective 3 deliverables
peff_ts_path    <- file.path(data_out_dir, "obj3_data_peff_timeseries_24yr.csv")
peff_sweep_path <- file.path(data_out_dir, "obj3_data_peff_sweep_all_models_24yr.csv")

if (!file.exists(peff_ts_path) || !file.exists(peff_sweep_path)) {
  stop("Missing Objective 3 files in Outputs/data. Run peff pipeline first.")
}

peff_ts_df    <- read_csv(peff_ts_path, show_col_types = FALSE)
peff_sweep_df <- read_csv(peff_sweep_path, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# STEP 1: FIT BASELINE SINGLE-ARENA MODEL
# ------------------------------------------------------------------------------
cat("Fitting Single-Arena baseline model (Z-scale)...\n")

sem_single <- paste0(
  target_growth, " ~ ", target_prey_ncc, "\n",
  target_cpue,   " ~ ", hake_col, "\n",
  target_sar,    " ~ ", target_growth, " + ", target_cpue
)

fit_single <- sem(sem_single, data = df_base,  fixed.x = TRUE, warn = FALSE)
#fit_single <- sem(sem_single, data = df_base, missing = "ML", fixed.x = FALSE, warn = FALSE)
aic_single <- as.numeric(fitMeasures(fit_single, "aic"))
r2_single  <- as.numeric(inspect(fit_single, "r2")[[target_sar]])
pe_single  <- parameterEstimates(fit_single, standardized = TRUE)

b_growth_single <- pe_single %>% filter(lhs == target_sar, rhs == target_growth) %>% pull(std.all)
b_cpue_single   <- pe_single %>% filter(lhs == target_sar, rhs == target_cpue) %>% pull(std.all)

cat(sprintf("Single-Arena Baseline: AIC = %.2f | R2(SAR) = %.2f | Beta_Growth = %.2f | Beta_CPUE = %.2f\n\n",
            aic_single, r2_single, b_growth_single, b_cpue_single))

# ------------------------------------------------------------------------------
# STEP 2: FIT DUAL- AND 2-TERM TRIPLE-ARENA MODELS
# ------------------------------------------------------------------------------
sig_ocean_candidates <- peff_sweep_df %>%
  filter(years_set == "24yr", mod_beta < 0, mod_p <= 0.05) %>%
  arrange(mod_aic)

results_list <- list()

for (i in seq_len(nrow(sig_ocean_candidates))) {
  cand <- sig_ocean_candidates[i, ]
  p_col <- cand$predator
  clean_pname <- gsub("^(X08_|X09_|X10_|X11_|X15_)", "", p_col)
  
  z_peff_col_name <- paste0("z_Peff_", clean_pname)
  if (!z_peff_col_name %in% names(peff_ts_df)) {
    z_peff_col_name <- paste0("Peff_", clean_pname)
    if (!z_peff_col_name %in% names(peff_ts_df)) next
  }
  
  df_fit <- df_base
  df_fit$P_ocean_eff        <- peff_ts_df[[z_peff_col_name]]
  df_fit$InRiver_SeaLions   <- df_base[[inriver_sealions_col]]
  df_fit$F_known_removals   <- df_base[[inriver_mort_col]]
  
  sem_base <- paste0(
    target_growth, " ~ ", target_prey_ncc, "\n",
    target_cpue,   " ~ ", hake_col, "\n"
  )
  
  # Dual-Arena (NCC Smolt + Subarctic Marine)
  sem_dual <- paste0(
    sem_base,
    target_sar, " ~ ", target_growth, " + ", target_cpue, " + ", target_prey_ak, " + P_ocean_eff"
  )
  fit_dual <- sem(sem_dual, data = df_fit,fixed.x = TRUE, warn = FALSE)
  aic_dual <- as.numeric(fitMeasures(fit_dual, "aic"))
  r2_dual  <- as.numeric(inspect(fit_dual, "r2")[[target_sar]])
  pe_dual  <- parameterEstimates(fit_dual, standardized = TRUE)
  p_ocean_dual <- pe_dual %>% filter(lhs == target_sar, rhs == "P_ocean_eff")
  
  # Triple-Arena with BOTH in-river terms:
  # 1. EMB Sea Lion abundance index
  # 2. Cumulative direct mortality rate (Bird + Tailrace Pinniped + Adult Harvest)
  sem_triple <- paste0(
    sem_base,
    target_sar, " ~ ", target_growth, " + ", target_cpue, " + ", target_prey_ak, " + ",
    "P_ocean_eff + InRiver_SeaLions + F_known_removals"
  )
  
  fit_triple <- tryCatch(
    sem(sem_triple, data = df_fit,fixed.x = TRUE, warn = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(fit_triple) && lavInspect(fit_triple, "converged")) {
    aic_triple <- as.numeric(fitMeasures(fit_triple, "aic"))
    r2_triple  <- as.numeric(inspect(fit_triple, "r2")[[target_sar]])
    pe_triple  <- parameterEstimates(fit_triple, standardized = TRUE)
    
    p_ocean_trip    <- pe_triple %>% filter(lhs == target_sar, rhs == "P_ocean_eff")
    p_sealions_trip <- pe_triple %>% filter(lhs == target_sar, rhs == "InRiver_SeaLions")
    p_removals_trip <- pe_triple %>% filter(lhs == target_sar, rhs == "F_known_removals")
    
    clean_label <- gsub("DFA_|_predAK", "", clean_pname)
    clean_label <- gsub("_smoltyr|_adultyr", "", clean_label)
    
    results_list[[length(results_list) + 1]] <- tibble(
      ocean_predator       = clean_label,
      mod_mechanism        = cand$mod_variant,
      alt_prey             = cand$mod_prey,
      
      # Single-Arena
      single_aic           = round(aic_single, 2),
      single_r2            = round(r2_single, 2),
      
      # Dual-Arena
      dual_aic             = round(aic_dual, 2),
      delta_aic_dual_v_sng = round(aic_dual - aic_single, 2),
      dual_r2              = round(r2_dual, 2),
      dual_b_ocean         = round(p_ocean_dual$std.all[1], 2),
      dual_p_ocean         = round(p_ocean_dual$pvalue[1], 3),
      
      # Triple-Arena
      triple_aic           = round(aic_triple, 2),
      delta_aic_trp_v_dual = round(aic_triple - aic_dual, 2),
      triple_r2            = round(r2_triple, 2),
      triple_b_ocean       = round(p_ocean_trip$std.all[1], 2),
      triple_p_ocean       = round(p_ocean_trip$pvalue[1], 3),
      
      # In-River Stage Estimates (2 parameters)
      triple_b_sealions    = round(p_sealions_trip$std.all[1], 2),
      triple_p_sealions    = round(p_sealions_trip$pvalue[1], 3),
      triple_b_removals    = round(p_removals_trip$std.all[1], 2),
      triple_p_removals    = round(p_removals_trip$pvalue[1], 3)
    )
  }
}

arena_df <- bind_rows(results_list) %>% arrange(triple_aic)


write_csv(arena_df, file.path(data_out_dir, "obj4_data_ocean_plus_CR_predation_2term_inriver.csv"))

# ------------------------------------------------------------------------------
# STEP 3: PUBLICATION GT TABLE
# ------------------------------------------------------------------------------
cat("Compiling GT Table for 2-term in-river stage...\n")

table_obj4_gt <- arena_df %>%
  mutate(
    Dual_Ocean_Disp = paste0(sprintf("%.2f", dual_b_ocean), " (", sprintf("%.2f", dual_p_ocean), ")"),
    Trip_Ocean_Disp = paste0(sprintf("%.2f", triple_b_ocean), " (", sprintf("%.2f", triple_p_ocean), ")"),
    Trip_SL_Disp    = paste0(sprintf("%.2f", triple_b_sealions), " (", sprintf("%.2f", triple_p_sealions), ")"),
    Trip_Rem_Disp   = paste0(sprintf("%.2f", triple_b_removals), " (", sprintf("%.2f", triple_p_removals), ")"),
    
    mod_mech_clean  = mod_mechanism,
    mod_mech_clean  = gsub("Joint Temp \\+ AltPrey", "Temp+Prey", mod_mech_clean),
    mod_mech_clean  = gsub("AltPrey Buffered", "Prey", mod_mech_clean),
    mod_mech_clean  = gsub("AltPrey Only", "Prey", mod_mech_clean),
    mod_mech_clean  = gsub("Temp Overlap Only", "Temp", mod_mech_clean),
    mod_mech_clean  = gsub("Temp Only", "Temp", mod_mech_clean)
  ) %>%
  select(
    ocean_predator, mod_mech_clean,
    Dual_Ocean_Disp, dual_r2, dual_aic, delta_aic_dual_v_sng,
    Trip_Ocean_Disp, Trip_SL_Disp, Trip_Rem_Disp,
    triple_r2, triple_aic, delta_aic_trp_v_dual
  ) %>%
  gt() %>%
  tab_header(
    title = md("**Table 4. Multi-Arena Life-Stage Bottlenecks: In-River Stage -- EXPLORATORY**"),
    subtitle = "Single (NCC), Dual (+ Subarctic), and Triple Arena (+ EMB Sea Lions & Measured Removals by birds (RTR), pinnipeds in tailrace (ACOE) and Zone 1-5 harvest (Joint Staffs Report))."
  ) %>%
  tab_spanner(label = "Dual-Arena (+ Subarctic Alaska)", columns = c(Dual_Ocean_Disp, dual_r2, dual_aic, delta_aic_dual_v_sng)) %>%
  tab_spanner(label = "Triple-Arena (+ In-River Sea Lions & Combined Measured Removals)", columns = c(Trip_Ocean_Disp, Trip_SL_Disp, Trip_Rem_Disp, triple_r2, triple_aic, delta_aic_trp_v_dual)) %>%
  cols_label(
    ocean_predator       = md("**Predator**"),
    mod_mech_clean       = md("**Gating**"),
    Dual_Ocean_Disp      = md("**P<sub>ocean</sub><br>&beta; (p)**"),
    dual_r2              = md("**R<sup>2</sup>**"),
    dual_aic             = md("**AIC**"),
    delta_aic_dual_v_sng = md("**&Delta;AIC<br>(D&minus;S)**"),
    Trip_Ocean_Disp      = md("**P<sub>ocean</sub><br>&beta; (p)**"),
    Trip_SL_Disp         = md("**Sea Lions<br>&beta; (p)**"),
    Trip_Rem_Disp        = md("**F<sub>removals</sub><br>&beta; (p)**"),
    triple_r2            = md("**R<sup>2</sup>**"),
    triple_aic           = md("**AIC**"),
    delta_aic_trp_v_dual = md("**&Delta;AIC<br>(T&minus;D)**")
  ) %>%
  cols_align(align = "left", columns = c(ocean_predator, mod_mech_clean)) %>%
  cols_align(align = "center", columns = -c(ocean_predator, mod_mech_clean)) %>%
  tab_style(
    style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#1B5E20", weight = "bold")),
    locations = cells_body(
      columns = c(delta_aic_trp_v_dual),
      rows = delta_aic_trp_v_dual < -2.0
    )
  ) %>%
  tab_options(
    table.font.size            = px(8.0),
    heading.title.font.size    = px(9.5),
    heading.subtitle.font.size = px(7.5),
    column_labels.font.size    = px(7.5),
    data_row.padding           = px(1.2)
  )

# Save deliverables
gtsave(table_obj4_gt, file.path(tbl_out_dir, "obj4_table1_ocean_plus_CR_predation_2term_inriver.pdf"))
gtsave(table_obj4_gt, file.path(tbl_out_dir, "obj4_table1_ocean_plus_CR_predation_2term_inriver.png"), vwidth = 740, vheight = 420)

cat("\nAnalysis complete! Tables and data saved to Outputs/.\n")

print(table_obj4_gt)


