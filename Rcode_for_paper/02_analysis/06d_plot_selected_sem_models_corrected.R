



#                     Model_ID    AIC   CFI  SRMR Growth_R2 CPUE_R2 SAR_R2
# 1               M1_NCC_Prey 187.81 0.934 0.032     0.890       -      -
# 2           M2_AK_Prey_Pink  59.39 1.000 0.000         -       -  0.480
# 3         M3_AK_Prey_Latent 330.48 0.475 0.257         -       -  0.959
# 4         M4_AK_Pred_Latent 426.18 0.573 0.137         -       -  0.513
# 5 M5_Confirmatory_Synthesis 638.22 0.612 0.159     0.896   0.749  0.518


# ==============================================================================
# Script: 02_analysis/06d_plot_selected_sem_models_corrected.R
# Purpose: Fit and plot all 5 selected SEM models with:
#          - Corrected paths: Hake -> CPUE (NCC), Latent_Peff -> SAR (AK)
#          - Dynamic R^2 annotations printed directly on endogenous salmon nodes
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(semPlot)
})

# ------------------------------------------------------------------------------
# 0. SETUP & DATA INGESTION
# ------------------------------------------------------------------------------
proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")
fig_dir  <- file.path(proj_dir, "Outputs/figures/selected_sem_plots")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

scaled_data_path <- file.path(data_dir, "sem_modeling_data_screened_zscored_1998_2021.csv")
if (!file.exists(scaled_data_path)) {
  stop("Scaled data file not found at: ", scaled_data_path)
}

df_fit <- read_csv(scaled_data_path, show_col_types = FALSE)

# Harmonize Hake naming
if (!"Hake_DFA" %in% names(df_fit)) {
  hake_col <- grep("hake", names(df_fit), value = TRUE, ignore.case = TRUE)
  if (length(hake_col) > 0) df_fit$Hake_DFA <- df_fit[[hake_col[1]]]
}

# ------------------------------------------------------------------------------
# 1. SPECIFY CORRECTED MODELS
# ------------------------------------------------------------------------------
models <- list(
  # Model 1: NCC Prey -> Growth
  "M1_NCC_Prey" = "
    Latent_sardine =~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr + 
                      X02_DFA_SeaNettle_smoltyr +
                      X04_marketsquid_GAM_smoltyr + 
                      X05_DFA_abundSardine_smoltyr
    Growth ~ Latent_sardine
  ",
  
  # Model 2: AK Prey (Primary Pink Competition + CPUE -> SAR)
  "M2_AK_Prey_Pink" = "
    SAR ~ CPUE + X14_pinkSalmonNorthAmerica_smoltyr
  ",
  
  # Model 3: AK Prey (Alternative Bipolar Latent + Growth + CPUE -> SAR)
  "M3_AK_Prey_Latent" = "
    my_latent =~ X12_biomassDecaShelfSum_smoltyr + 
                 X13_EGOA_DFA_seabirds_smoltyr + 
                 X12_DFA_EGOA_DFA_lowerTrophic_smoltyr + 
                 X12_biomassMysidShelfSum_smoltyr
    SAR ~ CPUE + Growth + my_latent
  ",
  
  # Model 4: AK Predators (Subarctic Latent Peff -> SAR, along with CPUE)
  "M4_AK_Pred_Latent" = "
    Latent_Peff =~ Peff_sablefishAK_biomass3plus_smoltyr + 
                   Peff_DFA_sleeperSharkBSAI_predAK_adultyr + 
                   Peff_DFA_sleeperSharkBSAI_predAK_smoltyr +
                   Peff_Harbour_s_WS_adultyr +
                   Peff_DFA_halibutBiomassAge8plus_predAK_adultyr 
                   
    SAR ~ CPUE + Latent_Peff
  ",
  
  # Model 5: Confirmatory Synthesis (NCC Prey -> Growth; Hake -> CPUE; Growth + CPUE + Latent_Peff -> SAR)
  "M5_Confirmatory_Synthesis" = "
    # NCC Early Marine Trophic & Predation Blocks
    Latent_sardine =~ X01_DFA_sumPreyOfPrey_planktonJun_smoltyr + 
                      X02_DFA_SeaNettle_smoltyr +
                      X04_marketsquid_GAM_smoltyr + 
                      X05_DFA_abundSardine_smoltyr
    Growth ~ Latent_sardine
    
    # Coastal Predation regulates smolt outmigration density (CPUE)
    CPUE ~ Hake_DFA

    # Subarctic Offshore Predation Complex
    Latent_Peff =~ Peff_sablefishAK_biomass3plus_smoltyr + 
                   Peff_DFA_sleeperSharkBSAI_predAK_adultyr + 
                   Peff_DFA_sleeperSharkBSAI_predAK_smoltyr +
                   Peff_Harbour_s_WS_adultyr +
                   Peff_DFA_halibutBiomassAge8plus_predAK_adultyr +
                   Peff_DFA_AllSeaLionsEMB_adultyr

    # Overall Marine Survival
    SAR ~ Growth + CPUE + Latent_Peff
  "
)

# ------------------------------------------------------------------------------
# 2. BASE LABELS DICTIONARY
# ------------------------------------------------------------------------------
base_labels <- c(
  "SAR"                                              = "Chinook\nSAR",
  "Growth"                                           = "Early Marine\nGrowth",
  "CPUE"                                             = "Smolt\nCPUE",
  "Hake_DFA"                                         = "NCC Hake\nPredation",
  
  "Latent_sardine"                                   = "Latent\nSardine Regime",
  "my_latent"                                        = "Latent\nAK Prey",
  "Latent_Peff"                                      = "Latent Subarctic\nEffective Predation",
  
  "X01_DFA_sumPreyOfPrey_planktonJun_smoltyr"        = "Plankton\n(Jun)",
  "X02_DFA_SeaNettle_smoltyr"                        = "Sea\nNettle",
  "X04_marketsquid_GAM_smoltyr"                      = "Market\nSquid",
  "X05_DFA_abundSardine_smoltyr"                     = "Sardine\nAbundance",
  
  "X14_pinkSalmonNorthAmerica_smoltyr"               = "Pink Salmon\n(Competitor)",
  "X12_biomassDecaShelfSum_smoltyr"                  = "Shelf Decapod\nBiomass",
  "X13_EGOA_DFA_seabirds_smoltyr"                    = "EGOA\nSeabirds",
  "X12_DFA_EGOA_DFA_lowerTrophic_smoltyr"            = "EGOA Lower\nTrophic DFA",
  "X12_biomassMysidShelfSum_smoltyr"                 = "Shelf Mysid\nBiomass",
  
  "Peff_sablefishAK_biomass3plus_smoltyr"            = "Sablefish\nAge-3+ Peff",
  "Peff_DFA_sleeperSharkBSAI_predAK_adultyr"         = "Sleeper Shark\nAdult Peff",
  "Peff_DFA_sleeperSharkBSAI_predAK_smoltyr"         = "Sleeper Shark\nSmolt Peff",
  "Peff_Harbour_s_WS_adultyr"                        = "Harbour Seal\nAdult Peff",
  "Peff_DFA_halibutBiomassAge8plus_predAK_adultyr"   = "Halibut Age-8+\nPeff",
  "Peff_DFA_AllSeaLionsEMB_adultyr"                  = "Sea Lions\nAdult Peff"
)

# ------------------------------------------------------------------------------
# 3. HELPER: DYNAMIC NODE LABELS WITH EMBEDDED R^2
# ------------------------------------------------------------------------------
get_model_labels <- function(fit_obj) {
  all_vars <- lavNames(fit_obj)
  lbls <- base_labels[all_vars]
  lbls[is.na(lbls)] <- all_vars[is.na(lbls)]
  
  # Extract R2 values
  r2_vals <- inspect(fit_obj, "r2")
  
  # Append R2 to the three focal salmon variables if they are endogenous in this model
  salmon_vars <- c("SAR", "Growth", "CPUE")
  for (v in salmon_vars) {
    if (v %in% names(r2_vals) && !is.na(r2_vals[v]) && r2_vals[v] > 0) {
      lbls[v] <- sprintf("%s\n[R² = %.2f]", lbls[v], r2_vals[v])
    }
  }
  return(lbls)
}

# ------------------------------------------------------------------------------
# 4. PLOTTING FUNCTION FOR STANDALONE MODELS
# ------------------------------------------------------------------------------
plot_sem_model <- function(fit_obj, model_title, file_name, width = 8.5, height = 6) {
  lbls <- get_model_labels(fit_obj)
  
  png(file_name, width = width, height = height, units = "in", res = 300)
  semPaths(
    fit_obj,
    whatLabels     = "std",
    layout         = "tree2",
    rotation       = 2,
    edge.label.cex = 1.05,
    edge.color     = "black",
    posCol         = "#1b7837",    # Positive: green
    negCol         = "#b2182b",    # Negative: red
    nodeLabels     = lbls,
    sizeMan        = 9.5,
    sizeMan2       = 6,
    sizeLat        = 11.5,
    sizeLat2       = 7.5,
    color          = list(man = "#deebf7", lat = "#fff7bc"),
    style          = "lisrel",
    residuals      = FALSE,
    intercepts     = FALSE,
    fade           = FALSE,
    mar            = c(4, 5, 4, 5)
  )
  title(main = model_title, cex.main = 1.15, font.main = 2)
  dev.off()
}

# ------------------------------------------------------------------------------
# 5. FIT AND RENDER ALL MODELS
# ------------------------------------------------------------------------------
fit_objects <- list()

model_meta <- tribble(
  ~m_id,                       ~title,                                                            ~w,   ~h,
  "M1_NCC_Prey",               "Model 1: NCC Prey Block (Sardine Regime -> Growth)",              8.0,  5.5,
  "M2_AK_Prey_Pink",           "Model 2: AK Prey Primary (Pink Salmon Competition -> SAR)",       7.5,  5.0,
  "M3_AK_Prey_Latent",         "Model 3: AK Prey Alternative (Shelf & EGOA Latent -> SAR)",       8.5,  6.0,
  "M4_AK_Pred_Latent",         "Model 4: AK Predator Block (Subarctic Peff Latent -> SAR)",       9.0,  6.0,
  "M5_Confirmatory_Synthesis", "Model 5: Multi-Stage Confirmatory Synthesis Model",               11.5, 8.0
)

for (i in seq_len(nrow(model_meta))) {
  m_id  <- model_meta$m_id[i]
  m_ttl <- model_meta$title[i]
  w_val <- model_meta$w[i]
  h_val <- model_meta$h[i]
  
  cat(sprintf("Fitting: %s ...\n", m_id))
  
  fit <- tryCatch(
    sem(
      models[[m_id]],
      data         = df_fit,
      missing      = "ML",
      bounds       = "pos.var",
      optim.method = "BFGS",
      control      = list(iter.max = 2000)
    ),
    error = function(e) {
      cat("Error fitting", m_id, ":", e$message, "\n")
      return(NULL)
    }
  )
  
  if (!is.null(fit) && lavInspect(fit, "converged")) {
    fit_objects[[m_id]] <- fit
    out_file <- file.path(fig_dir, paste0("fig_", m_id, ".png"))
    plot_sem_model(fit, m_ttl, out_file, width = w_val, height = h_val)
    cat(sprintf("Saved plot: %s\n", out_file))
  }
}

# ------------------------------------------------------------------------------
# 6. GENERATE 4-PANEL REGIONAL FIGURE WITH R^2 ANNOTATIONS
# ------------------------------------------------------------------------------
cat("\nGenerating 4-panel multi-model compilation figure...\n")
comp_path <- file.path(fig_dir, "fig_compilation_regional_blocks_M1_to_M4.png")
png(comp_path, width = 16, height = 12, units = "in", res = 300)
layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))

plot_panel <- function(fit_obj, panel_title) {
  lbls <- get_model_labels(fit_obj)
  semPaths(
    fit_obj,
    whatLabels     = "std",
    layout         = "tree2",
    rotation       = 2,
    edge.label.cex = 1.1,
    edge.color     = "black",
    posCol         = "#1b7837",
    negCol         = "#b2182b",
    nodeLabels     = lbls,
    sizeMan        = 10,
    sizeMan2       = 6,
    sizeLat        = 12,
    sizeLat2       = 7.5,
    color          = list(man = "#deebf7", lat = "#fff7bc"),
    style          = "lisrel",
    residuals      = FALSE,
    intercepts     = FALSE,
    fade           = FALSE,
    mar            = c(4, 5, 4, 5)
  )
  title(main = panel_title, cex.main = 1.3, font.main = 2)
}

plot_panel(fit_objects[["M1_NCC_Prey"]],       "(A) NCC Prey Block: Coastal Upwelling / Sardine Regime")
plot_panel(fit_objects[["M2_AK_Prey_Pink"]],   "(B) AK Prey Primary: Subarctic Pink Salmon Competition")
plot_panel(fit_objects[["M3_AK_Prey_Latent"]], "(C) AK Prey Alternative: Shelf Macrozooplankton & Forage")
plot_panel(fit_objects[["M4_AK_Pred_Latent"]], "(D) AK Predator Block: Demersal Fish & Marine Mammal Peff")

dev.off()
cat(sprintf("Saved multi-panel plot: %s\n", comp_path))

# ------------------------------------------------------------------------------
# 7. PRINT CONCISE FIT & R^2 AUDIT TABLE
# ------------------------------------------------------------------------------
cat("\n====================================================================\n")
cat("SUMMARY TABLE: FIT METRICS AND R^2 FOR ALL 5 MODELS\n")
cat("====================================================================\n")

fit_table <- tibble(
  Model_ID = names(fit_objects),
  AIC      = map_dbl(fit_objects, ~ round(as.numeric(fitMeasures(.x, "aic")), 2)),
  CFI      = map_dbl(fit_objects, ~ round(as.numeric(fitMeasures(.x, "cfi")), 3)),
  SRMR     = map_dbl(fit_objects, ~ round(as.numeric(fitMeasures(.x, "srmr")), 3)),
  Growth_R2= map_chr(fit_objects, ~ {
    r2 <- inspect(.x, "r2")
    if ("Growth" %in% names(r2) && !is.na(r2["Growth"])) sprintf("%.3f", r2["Growth"]) else "-"
  }),
  CPUE_R2  = map_chr(fit_objects, ~ {
    r2 <- inspect(.x, "r2")
    if ("CPUE" %in% names(r2) && !is.na(r2["CPUE"])) sprintf("%.3f", r2["CPUE"]) else "-"
  }),
  SAR_R2   = map_chr(fit_objects, ~ {
    r2 <- inspect(.x, "r2")
    if ("SAR" %in% names(r2) && !is.na(r2["SAR"])) sprintf("%.3f", r2["SAR"]) else "-"
  })
)

print(as.data.frame(fit_table))
write_csv(fit_table, file.path(data_dir, "obj4_selected_models_fit_and_r2_summary.csv"))

#fix model 5 plot-----------
plot_sem_model_fixed <- function(fit_obj, model_title, file_name, width = 11, height = 7.5) {
  # 1. Extract exact node names in semPlot's native order via @Vars$name
  sem_mod <- semPlotModel(fit_obj)
  node_names <- sem_mod@Vars$name
  
  # 2. Extract model R2
  r2_vals <- inspect(fit_obj, "r2")
  
  # 3. Match base labels strictly to semPlot's node sequence
  lbls <- base_labels[node_names]
  lbls[is.na(lbls)] <- node_names[is.na(lbls)]
  
  # 4. Append R2 to the endogenous salmon nodes
  salmon_vars <- c("SAR", "Growth", "CPUE")
  for (i in seq_along(node_names)) {
    v <- node_names[i]
    if (v %in% salmon_vars && v %in% names(r2_vals) && !is.na(r2_vals[v]) && r2_vals[v] > 0) {
      lbls[i] <- sprintf("%s\n[R² = %.2f]", lbls[i], r2_vals[v])
    }
  }
  
  png(file_name, width = width, height = height, units = "in", res = 300)
  semPaths(
    fit_obj,
    whatLabels     = "std",
    layout         = "tree2",
    rotation       = 2,
    edge.label.cex = 1.05,
    edge.color     = "black",
    posCol         = "#1b7837",
    negCol         = "#b2182b",
    nodeLabels     = lbls,          # Guaranteed 1:1 match with sem_mod@Vars$name
    sizeMan        = 9.5,
    sizeMan2       = 6,
    sizeLat        = 11.5,
    sizeLat2       = 7.5,
    color          = list(man = "#deebf7", lat = "#fff7bc"),
    style          = "lisrel",
    residuals      = FALSE,
    intercepts     = FALSE,
    fade           = FALSE,
    mar            = c(4, 5, 4, 5)
  )
  title(main = model_title, cex.main = 1.15, font.main = 2)
  dev.off()
}

# Re-run Model 5
plot_sem_model_fixed(
  fit_objects[["M5_Confirmatory_Synthesis"]], 
  "Model 5: Multi-Stage Confirmatory Synthesis Model", 
  file.path(fig_dir, "fig_M5_Confirmatory_Synthesis.png")
)


