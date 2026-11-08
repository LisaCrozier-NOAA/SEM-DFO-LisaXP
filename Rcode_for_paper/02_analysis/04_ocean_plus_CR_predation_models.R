# ==============================================================================
# Script: 02_analysis/04_ocean_plus_CR_predation_models.R
# Purpose: Objective 4 - Multi-Basin Life-Stage Bottlenecks (24-Year Complete Cases)
#          Single: SAR ~ Growth + CPUE
#          Dual:   SAR ~ Growth + CPUE + PreyAK + P_ocean_eff
#          Triple: SAR ~ Growth + CPUE + PreyAK + P_ocean_eff + InRiver_SeaLions
# Outputs:
#   - Outputs/data/obj4_data_ocean_plus_CR_predation_24yr.csv
#   - Outputs/tables/obj4_table1_ocean_plus_CR_predation_24yr.[pdf|png]
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

raw_csv_path <- file.path(proj_dir, "metadata", "sem_altprey_data_1998_2021_10212026.csv")
if (!file.exists(raw_csv_path)) {
  stop("Missing master raw altprey dataset at: ", raw_csv_path)
}

raw_df_master <- read.csv(raw_csv_path, stringsAsFactors = FALSE, row.names = 1)
raw_df_master <- raw_df_master[, names(raw_df_master) != "" & !is.na(names(raw_df_master))]

# Complete cases for 24-year series
valid_cols <- names(raw_df_master)[colSums(!is.na(raw_df_master)) == 24]
df_raw <- raw_df_master %>% 
  select(all_of(valid_cols)) %>% 
  filter(complete.cases(.))

scale_custom_bounds <- function(x, a = 0.2, b = 1.0) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}

df_base <- df_raw %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), ~ scale_custom_bounds(.x, 0.0, 1.0)))

# Structural variable anchors
target_growth    <- "X06_DFA_IGF_mu_smoltyr"
target_prey_ncc  <- if ("X05_DFA_abundSardine_smoltyr" %in% names(df_base)) "X05_DFA_abundSardine_smoltyr" else "X04_marketsquid_GAM_smoltyr"
target_cpue      <- "X07_cpue_IntSprJunHW_smoltyr"
hake_col         <- "X09_DFA_HakeAge5Plus_smoltyr"
target_prey_ak   <- "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"
target_sar       <- "X16_DFA_sar_UCsprCh_smoltyr"
inriver_pred_col <- "X10_DFA_AllSeaLionsEMB_adultyr"

# Load Objective 3 deliverables
peff_ts_path    <- file.path(data_out_dir, "obj3_data_peff_timeseries_24yr.csv")
peff_sweep_path <- file.path(data_out_dir, "obj3_data_peff_sweep_all_models_24yr.csv")

if (!file.exists(peff_ts_path) || !file.exists(peff_sweep_path)) {
  stop("Missing Objective 3 deliverables. Run script 03_effective_predation_peff_pipeline.R first.")
}

peff_ts_df    <- read_csv(peff_ts_path, show_col_types = FALSE)
peff_sweep_df <- read_csv(peff_sweep_path, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# STEP 1: FIT BASELINE SINGLE-ARENA MODEL (NCC ONLY)
# ------------------------------------------------------------------------------
cat("Fitting Single-Arena baseline model (Early Marine NCC only)...\n")

sem_single <- paste0(
  target_growth, " ~ ", target_prey_ncc, "\n",
  target_cpue,   " ~ ", hake_col, "\n",
  target_sar,    " ~ ", target_growth, " + ", target_cpue
)
fit_single <- sem(sem_single, data = df_base, missing = "ML", warn = FALSE)
aic_single <- as.numeric(fitMeasures(fit_single, "aic"))
r2_single  <- as.numeric(inspect(fit_single, "r2")[[target_sar]])
pe_single  <- parameterEstimates(fit_single, standardized = TRUE)

b_growth_single <- pe_single %>% filter(lhs == target_sar, rhs == target_growth) %>% pull(std.all)
b_cpue_single   <- pe_single %>% filter(lhs == target_sar, rhs == target_cpue) %>% pull(std.all)

cat(sprintf("Single-Arena Baseline: AIC = %.2f | R2(SAR) = %.2f | Beta_Growth = %.2f | Beta_CPUE = %.2f\n\n",
            aic_single, r2_single, b_growth_single, b_cpue_single))

# ------------------------------------------------------------------------------
# STEP 2: FILTER STATISTICALLY SIGNIFICANT (p <= 0.05) SUBARCTIC PREDATORS
# ------------------------------------------------------------------------------
sig_ocean_candidates <- peff_sweep_df %>%
  filter(years_set == "24yr", mod_beta < 0, mod_p <= 0.05) %>%
  arrange(mod_aic)

cat(sprintf("Evaluating %d statistically significant (p <= 0.05) subarctic predators across arenas.\n", 
            nrow(sig_ocean_candidates)))

# ------------------------------------------------------------------------------
# STEP 3: FIT DUAL- AND TRIPLE-ARENA MODELS
# ------------------------------------------------------------------------------
results_list <- list()

for (i in seq_len(nrow(sig_ocean_candidates))) {
  cand <- sig_ocean_candidates[i, ]
  p_col <- cand$predator
  clean_pname <- gsub("^(X08_|X09_|X10_|X11_|X15_)", "", p_col)
  peff_col_name <- paste0("Peff_", clean_pname)
  
  if (!peff_col_name %in% names(peff_ts_df)) next
  
  df_fit <- df_base
  df_fit$P_ocean_eff  <- peff_ts_df[[peff_col_name]]
  df_fit$InRiver_Pred <- df_base[[inriver_pred_col]]
  
  sem_base <- paste0(
    target_growth, " ~ ", target_prey_ncc, "\n",
    target_cpue,   " ~ ", hake_col, "\n"
  )
  
  # Dual-Arena (NCC Smolt + Subarctic Marine)
  sem_dual <- paste0(
    sem_base,
    target_sar, " ~ ", target_growth, " + ", target_cpue, " + ", target_prey_ak, " + P_ocean_eff"
  )
  fit_dual <- sem(sem_dual, data = df_fit, missing = "ML", warn = FALSE)
  aic_dual <- as.numeric(fitMeasures(fit_dual, "aic"))
  r2_dual  <- as.numeric(inspect(fit_dual, "r2")[[target_sar]])
  pe_dual  <- parameterEstimates(fit_dual, standardized = TRUE)
  p_ocean_dual <- pe_dual %>% filter(lhs == target_sar, rhs == "P_ocean_eff")
  
  # Triple-Arena (NCC Smolt + Subarctic Marine + River Adult Sea Lions)
  sem_triple <- paste0(
    sem_base,
    target_sar, " ~ ", target_growth, " + ", target_cpue, " + ", target_prey_ak, " + P_ocean_eff + InRiver_Pred"
  )
  fit_triple <- sem(sem_triple, data = df_fit, missing = "ML", warn = FALSE)
  
  if (!is.null(fit_triple) && lavInspect(fit_triple, "converged")) {
    aic_triple <- as.numeric(fitMeasures(fit_triple, "aic"))
    r2_triple  <- as.numeric(inspect(fit_triple, "r2")[[target_sar]])
    pe_triple  <- parameterEstimates(fit_triple, standardized = TRUE)
    
    p_ocean_trip   <- pe_triple %>% filter(lhs == target_sar, rhs == "P_ocean_eff")
    p_inriver_trip <- pe_triple %>% filter(lhs == target_sar, rhs == "InRiver_Pred")
    
    clean_label <- gsub("DFA_|_predAK", "", clean_pname)
    clean_label <- gsub("_smoltyr|_adultyr", "", clean_label)
    
    results_list[[length(results_list) + 1]] <- tibble(
      ocean_predator       = clean_label,
      mod_mechanism        = cand$mod_variant,
      alt_prey             = cand$mod_prey,
      
      # Single-Arena Baseline Reference
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
      triple_b_river       = round(p_inriver_trip$std.all[1], 2),
      triple_p_river       = round(p_inriver_trip$pvalue[1], 3)
    )
  }
}

arena_df <- bind_rows(results_list) %>% arrange(triple_aic)

# Save Clean CSV
write_csv(arena_df, file.path(data_out_dir, "obj4_data_ocean_plus_CR_predation_24yr.csv"))

# ------------------------------------------------------------------------------
# STEP 4: CONSTRUCT PUBLICATION GT TABLE (TIGHT PORTRAIT FIT)
# ------------------------------------------------------------------------------
cat("Compiling publication table for Objective 4...\n")

table_obj4_gt <- arena_df %>%
  mutate(
    Dual_Ocean_Disp = paste0(sprintf("%.2f", dual_b_ocean), " (", sprintf("%.2f", dual_p_ocean), ")"),
    Trip_Ocean_Disp = paste0(sprintf("%.2f", triple_b_ocean), " (", sprintf("%.2f", triple_p_ocean), ")"),
    Trip_River_Disp = paste0(sprintf("%.2f", triple_b_river), " (", sprintf("%.2f", triple_p_river), ")"),
    
    # Clean mechanism strings: remove "Joint", "Only", and "Overlap"
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
    Trip_Ocean_Disp, Trip_River_Disp, triple_r2, triple_aic, delta_aic_trp_v_dual
  ) %>%
  gt() %>%
  tab_header(
    title = md("**Table 4. Life-Stage Arena Bottlenecks: Single-, Dual-, and Triple-Arena SEM Comparisons**"),
    subtitle = "Evaluating smolt-to-adult survival (SAR) across life stages: Coastal NCC alone, adding Subarctic Alaska marine transit, and adding Columbia River adult pinnipeds."
  ) %>%
  tab_spanner(label = "Dual-Arena (+ Subarctic Alaska)", columns = c(Dual_Ocean_Disp, dual_r2, dual_aic, delta_aic_dual_v_sng)) %>%
  tab_spanner(label = "Triple-Arena (+ Adult Columbia River Sea Lions)", columns = c(Trip_Ocean_Disp, Trip_River_Disp, triple_r2, triple_aic, delta_aic_trp_v_dual)) %>%
  cols_label(
    ocean_predator       = md("**Predator**"),
    mod_mech_clean       = md("**Gating**"),
    Dual_Ocean_Disp      = md("**P<sub>ocean</sub><br>&beta; (p)**"),
    dual_r2              = md("**R<sup>2</sup>**"),
    dual_aic             = md("**AIC**"),
    delta_aic_dual_v_sng = md("**&Delta;AIC<br>(D&minus;S)**"),
    Trip_Ocean_Disp      = md("**P<sub>ocean</sub><br>&beta; (p)**"),
    Trip_River_Disp      = md("**P<sub>river</sub><br>&beta; (p)**"),
    triple_r2            = md("**R<sup>2</sup>**"),
    triple_aic           = md("**AIC**"),
    delta_aic_trp_v_dual = md("**&Delta;AIC<br>(T&minus;D)**")
  ) %>%
  cols_align(align = "left", columns = c(ocean_predator, mod_mech_clean)) %>%
  cols_align(align = "center", columns = c(Dual_Ocean_Disp, dual_r2, dual_aic, delta_aic_dual_v_sng,
                                           Trip_Ocean_Disp, Trip_River_Disp, triple_r2, triple_aic, delta_aic_trp_v_dual)) %>%
  cols_width(
    ocean_predator       ~ px(95),
    mod_mech_clean       ~ px(65),
    Dual_Ocean_Disp      ~ px(62),
    dual_r2              ~ px(34),
    dual_aic             ~ px(42),
    delta_aic_dual_v_sng ~ px(45),
    Trip_Ocean_Disp      ~ px(62),
    Trip_River_Disp      ~ px(62),
    triple_r2            ~ px(34),
    triple_aic           ~ px(42),
    delta_aic_trp_v_dual ~ px(45)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#1B5E20", weight = "bold")),
    locations = cells_body(
      columns = c(delta_aic_dual_v_sng),
      rows = delta_aic_dual_v_sng < -2.0
    )
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#E8F5E9"), cell_text(color = "#1B5E20", weight = "bold")),
    locations = cells_body(
      columns = c(delta_aic_trp_v_dual, Trip_River_Disp),
      rows = delta_aic_trp_v_dual < -2.0
    )
  ) %>%
  tab_source_note(
    source_note = md(sprintf(
      "**Single-Arena Baseline (Coastal NCC only):** AIC = %.2f, R<sup>2</sup>(SAR) = %.2f (Growth &beta; = %.2f, CPUE &beta; = %.2f). &Delta;AIC<sub>(D&minus;S)</sub> compares Dual-Arena (adding Alaska marine) to Single-Arena; &Delta;AIC<sub>(T&minus;D)</sub> compares Triple-Arena (adding adult in-river pinnipeds) to Dual-Arena. Shaded cells denote &Delta;AIC &le; &minus;2.0.",
      aic_single, r2_single, b_growth_single, b_cpue_single
    ))
  ) %>%
  tab_options(
    table.font.size            = px(8.0),
    heading.title.font.size    = px(9.5),
    heading.subtitle.font.size = px(7.5),
    column_labels.font.size    = px(7.5),
    column_labels.font.weight  = "bold",
    row_group.font.size        = px(8.0),
    data_row.padding           = px(1.2),
    column_labels.padding      = px(2),
    source_notes.font.size     = px(6.8),
    source_notes.padding       = px(3)
  )

# Export PDF and PNG with compact width (~620px canvas)
gtsave(table_obj4_gt, file.path(tbl_out_dir, "obj4_table1_ocean_plus_CR_predation_24yr.pdf"))
gtsave(table_obj4_gt, file.path(tbl_out_dir, "obj4_table1_ocean_plus_CR_predation_24yr.png"), vwidth = 620, vheight = 420)

cat("\nObjective 4 complete! Files exported to Outputs/data and Outputs/tables.\n")