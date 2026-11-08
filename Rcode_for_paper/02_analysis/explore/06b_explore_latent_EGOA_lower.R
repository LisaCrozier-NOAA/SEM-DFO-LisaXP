# ==============================================================================
# Script: 02_analysis/06b_explore_latent_EGOA_lower.R
# Purpose: Exploratory single-node latent SEM for EGOA Lower Trophic constituents.
#          Holds all other nodes fixed at their top configuration (DAG 1B + Sablefish Peff).
#          Backward prunes non-significant factor loadings (=~) one at a time.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# 0. SETUP & LOAD DATA
# ------------------------------------------------------------------------------
proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")

ak_const_path <- file.path(data_dir, "ak_dfa_raw_constituents_1998_2021.csv")
meta_path     <- file.path(data_dir, "ak_dfa_constituents_metadata_summary.csv")
alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")

if (!file.exists(ak_const_path)) stop("Missing: ", ak_const_path)
if (!file.exists(meta_path))     stop("Missing: ", meta_path)
if (!file.exists(alt_data_path)) stop("Missing: ", alt_data_path)
if (!file.exists(peff_ts_path))  stop("Missing: ", peff_ts_path)

ak_const_df <- read_csv(ak_const_path, show_col_types = FALSE)
meta_df     <- read_csv(meta_path, show_col_types = FALSE)
peff_df     <- read_csv(peff_ts_path, show_col_types = FALSE)

# ==============================================================================
# Helper Function: Auto-Extract Top DFA Constituents
# ==============================================================================
get_top_dfa_indicators <- function(dfa_loadings_csv,
                                   region = c("EGOA", "WGOA"),
                                   loading_threshold = 0.20,
                                   available_cols) {
  
  region <- match.arg(region)
  
  egoa_dict <- c(
    # LTL
    "Chla biomass"               = "CHLA.BIOM",
    "Chla peak bloom"            = "CHLA.PEAK",
    "Copepod community size"     = "SECM.COP.COMM.SIZE",
    "Euphausiid abundance"       = "SECM.EUPH.DENS",
    "Gastropod abundance"        = "SECM.GAST.DENS",
    "Amphipod abundance"         = "SECM.HYP.AMPH",
    "Large copepod abundance"    = "SECM.LG.CAL.COP.DENS",
    "Mesozooplankton abundance"  = "SECM.MESOZOOP.DENS",
    "Small copepod abundance"    = "SECM.SM.CAL.COP.DENS",
    # MTL
    "Herring biomass (Craig)"    = "CRG.HERR.MATBIOM",
    "Humpback whale birth date"  = "GB.HUMP.CRUD.BR.JA",
    "Juv chinook salmon cpue"    = "SECM.CHIN.CPUE",
    "Juv chum salmon cpue"       = "SECM.CHUM.CPUE",
    "Juv coho salmon cpue"       = "SECM.COHO.CPUE",
    "Juv pink salmon cpue"       = "SECM.PINK.CPUE",
    "Juv sockeye salmon cpue"    = "SECM.SOCK.CPUE",
    "Herring biomass (Sitka)"    = "STKA.HERR.MATBIOM",
    # Seabirds
    "Common murre prod."         = "SB.REPSUC.COMU.LAZ",
    "Fork-tailed storm petrel prod." = "SB.REPSUC.FTSP.LAZ",
    "Leach's storm petrel prod." = "SB.REPSUC.LSP.LAZ",
    "Thick-billed murre prod."   = "SB.REPSUC.TBMU.LAZ",
    "Common murre hatch"         = "SBRD.MEANHATCH.COMU.LAZ",
    "Fork-tailed storm petrel hatch" = "SBRD.MEANHATCH.FTSP.LAZ",
    "Glaucous-winged gull hatch" = "SBRD.MEANHATCH.GWGU.LAZ",
    "Leach's storm petrel hatch " = "SBRD.MEANHATCH.LHSP.LAZ",
    "Rhinocerus auklet hatch"    = "SBRD.MEANHATCH.RHAU.LAZ",
    "Thick billed murre hatch"   = "SBRD.MEANHATCH.TBMU.LAZ"
  )
  
  wgoa_dict <- c(
    # LTL
    "Chla biomass"               = "CHLA.BIOM",
    "Chla peak bloom"            = "CHLA.PEAK",
    "Larval sandlance"           = "FOCI.ICH.AMMPER",
    "Larval arrowtooth"          = "FOCI.ICH.ATHSTO",
    "Larval ronquils"            = "FOCI.ICH.BATSPP",
    "Larval pollock"             = "FOCI.ICH.GADCHA",
    "Larval Alaska pollock"      = "FOCI.ICH.GADMAC",
    "Larval flathead sole"       = "FOCI.ICH.HIPELA",
    "Larval P. halibut"          = "FOCI.ICH.HIPSTE",
    "Larval rock sole"           = "FOCI.ICH.LEPBIL",
    "Larval nor. rock sole"      = "FOCI.ICH.LEPPOL",
    "Larval starry flounder"     = "FOCI.ICH.PLASTE",
    "Larval rockfishes"          = "FOCI.ICH.SEBSPP",
    "Larval nor. lampfish"       = "FOCI.ICH.STELEU",
    "Copepod community size"     = "SWDLN.COP.CMMTY.SZ.SPR",
    "Lg calanoid copepod fall"   = "SWDLN.LG.CAL.COP.FLL",
    "Lg calanoid copepod spr"    = "SWDLN.LG.CAL.COP.SPR",
    "Sm calanoid copepod fall"   = "SWDLN.SM.CAL.COP.FLL",
    "Sm calanoid copepod spr"    = "SWDLN.SM.CAL.COP.SPR",
    # MTL
    "Chiniak eulachon"           = "CHI.EULACHON.CPUE",
    "Chiniak jellyfish"          = "CHI.JELLY.CPUE",
    "Chiniak pink shrimp"        = "CHI.SHRIMP.CPUE",
    "Prop. capelin in diet"      = "MID.IL.CAPELIN",
    "Pavlof shrimp"              = "PAV.SHRIMP.CPUE",
    "Pavlof capelin"             = "PAV.CAPELIN.CPUE",
    "Pavlof Eulachon"            = "PAV.EULACHON.CPUE",
    "Pavlof Jellyfish"           = "PAV.JELLY.CPUE",
    "Sea lion pup counts"        = "SSL.WEST.PUP.PRED",
    # Seabirds
    "Common murre hatch"         = "SB.MEANHATCH.COMU.CHOW",
    "Thick-billed murre hatch"   = "SB.MEANHATCH.TBMU.CHOW",
    "Parakeet auklet hatch"      = "SB.MEANHATCH.PAAU.CHOW",
    "Horned puffin hatch"        = "SB.MEANHATCH.HOPU.CHOW",
    "Tufted puffin hatch"        = "SB.MEANHATCH.TUPU.CHOW",
    "Glaucous-winged gull hatch" = "SB.MEANHATCH.GWGU.CHOW",
    "Kittiwake hatch"            = "SB.MEANHATCH.BLKW.EAMT",
    "Common murre prod."         = "SB.REPSUC.COMU.CHOW",
    "Thick-billed murre prod."   = "SB.REPSUC.TBMU.CHOW",
    "Horned puffin prod."        = "SB.REPSUC.HPUFF.CHOW",
    "Tufted puffin prod."        = "SB.REPSUC.TPUF.CHOW",
    "Kittiwake prod."            = "SB.REPSUC.BLK.CHOW",
    "Parakeet auklet prod."      = "SB.REPSUC.PAUK.CHOW",
    "Black-legged kittiwake prod." = "SB.REPSUC.BLK.AMAT",
    "Fork-tailed storm petrel prod." = "SB.REPSUC.FTSP.AMAT"
  )
  
  dict <- if (region == "EGOA") egoa_dict else wgoa_dict
  
  df_loadings <- read_csv(dfa_loadings_csv, show_col_types = FALSE)
  load_col <- grep("^loading", names(df_loadings), value = TRUE)[1]
  
  filtered <- df_loadings %>%
    mutate(
      abs_loading = abs(.data[[load_col]]),
      raw_code    = dict[trimws(names)],
      col_name    = paste0(region, "_", raw_code)
    ) %>%
    filter(!is.na(raw_code), abs_loading >= loading_threshold) %>%
    arrange(desc(abs_loading))
  
  valid_indicators <- filtered %>%
    filter(col_name %in% available_cols)
  
  cat(sprintf("\n[Auto-Extractor] %s: Found %d indicators with |loading| >= %.2f\n",
              basename(dfa_loadings_csv), nrow(valid_indicators), loading_threshold))
  
  for (i in seq_len(nrow(valid_indicators))) {
    cat(sprintf("  -> %-32s (Loading: %6.3f) -> Column: %s\n",
                valid_indicators$names[i],
                valid_indicators[[load_col]][i],
                valid_indicators$col_name[i]))
  }
  
  return(valid_indicators$col_name)
}

# ------------------------------------------------------------------------------
# 1. SELECT CANDIDATES & ASSEMBLE Z-SCORED DATASET
# ------------------------------------------------------------------------------
# 1a. Filter to candidate EGOA Lower Trophic variables with <= 2 missing years
egoa_ltl_vars <- meta_df %>%
  filter(region == "EGOA", subguild == "Lower Trophic", missing_years <= 3) %>%
  pull(indicator_col)

cat(sprintf("Found %d candidate EGOA Lower Trophic indicators in meta table:\n", length(egoa_ltl_vars)))

df_ak_raw <- ak_const_df %>%
  select(year, all_of(egoa_ltl_vars))

# 1b. Load DAG 1B fixed anchors across full 24-year series
df_base_all <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)

yr_col <- grep("^year$", names(df_base_all), value = TRUE, ignore.case = TRUE)
df_base_all$year <- if (length(yr_col) > 0) {
  as.numeric(df_base_all[[yr_col[1]]])
} else {
  as.numeric(rownames(df_base_all))
}

df_anchors <- df_base_all %>%
  select(
    year,
    Growth  = X06_DFA_IGF_mu_smoltyr,
    CPUE    = X07_cpue_IntSprJunHW_smoltyr,
    SAR     = X16_DFA_sar_UCsprCh_smoltyr,
    PreyNCC = X05_DFA_abundSardine_smoltyr,
    PredNCC = X09_DFA_HakeAge5Plus_smoltyr
  ) %>%
  filter(year >= 1998 & year <= 2021)

df_anchors$P_ocean_eff <- peff_df$Peff_sablefishAK_biomass3plus_smoltyr

cat(sprintf("Anchors loaded: %d years (%d - %d)\n", 
            nrow(df_anchors), min(df_anchors$year), max(df_anchors$year)))

# 1c. Merge with AK constituents and apply z-scoring
df_fit <- df_anchors %>%
  left_join(df_ak_raw, by = "year") %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

# 1d. Auto-extract top DFA indicators (> 0.20 loading)
egoa_ltl_csv <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Ferris-DFAIndicators-goa/MARSS results/EGOA/egoa_ltl_loadings.csv"

current_indicators <- get_top_dfa_indicators(
  dfa_loadings_csv  = egoa_ltl_csv,
  region            = "EGOA",
  loading_threshold = 0.20,
  available_cols    = names(df_fit)
)

# ------------------------------------------------------------------------------
# 2. STEPWISE BACKWARD PRUNING ENGINE (Focal: Latent_PreyAK)
# ------------------------------------------------------------------------------
build_syntax <- function(indicators) {
  meas_block <- paste0("Latent_PreyAK =~ ", paste(indicators, collapse = " + "))
  
  struct_block <- paste0(
    "Growth ~ PreyNCC\n",
    "CPUE   ~ PredNCC\n",
    "SAR    ~ Growth + CPUE + Latent_PreyAK + P_ocean_eff"
  )
  
  paste(meas_block, struct_block, sep = "\n")
}

step <- 1
step_history <- list()

cat("\n====================================================================\n")
cat("STARTING STEPWISE PRUNING FOR: Latent_PreyAK (EGOA Lower Trophic)\n")
cat("====================================================================\n")

repeat {
  cat(sprintf("\n--- STEP %d: Testing %d Indicators ---\n", step, length(current_indicators)))
  print(current_indicators)
  
  syntax <- build_syntax(current_indicators)
  
  fit_err <- NULL
  fit <- tryCatch(
    sem(
      syntax, 
      data    = df_fit, 
      std.lv  = TRUE, 
      missing = "ML", 
      bounds  = "pos.var",          # Prevents Heywood cases / negative variance crashes
      optim.method = "BFGS",        # More forgiving on small-sample surfaces
      control = list(iter.max = 2000)
    ),
    error = function(e) {
      fit_err <<- e$message
      return(NULL)
    }
  )
  
  if (is.null(fit)) {
    cat("Model failed to fit. Error message:\n", fit_err, "\n")
    break
  }
  
  if (!lavInspect(fit, "converged")) {
    cat("Warning: Model failed to converge at step", step, "\n")
    break
  }
  
  # Fit indices
  aic_val  <- round(as.numeric(fitMeasures(fit, "aic")), 2)
  cfi_val  <- round(as.numeric(fitMeasures(fit, "cfi")), 3)
  srmr_val <- round(as.numeric(fitMeasures(fit, "srmr")), 3)
  p_chi2   <- round(as.numeric(fitMeasures(fit, "pvalue")), 4)
  r2_sar   <- round(as.numeric(inspect(fit, "r2")["SAR"]), 3)
  
  pe <- parameterEstimates(fit, standardized = TRUE)
  
  # Structural path: Latent_PreyAK -> SAR
  sar_path  <- pe %>% filter(lhs == "SAR", rhs == "Latent_PreyAK")
  path_beta <- round(sar_path$std.all[1], 3)
  path_p    <- round(sar_path$pvalue[1], 4)
  
  cat(sprintf("Fit: AIC = %.2f | p(Chi2) = %.4f | CFI = %.3f | SRMR = %.3f | SAR R2 = %.3f\n",
              aic_val, p_chi2, cfi_val, srmr_val, r2_sar))
  cat(sprintf("Structural Path: SAR ~ Latent_PreyAK : beta = %.3f (p = %.4f)\n", path_beta, path_p))
  
  # Factor loadings
  loadings <- pe %>%
    filter(op == "=~", lhs == "Latent_PreyAK") %>%
    select(rhs, std.all, pvalue) %>%
    arrange(desc(pvalue))
  
  cat("\nFactor Loadings (std.all and p-value):\n")
  print(as.data.frame(loadings))
  
  # Record step
  step_history[[step]] <- tibble(
    step         = step, 
    n_indicators = length(current_indicators),
    AIC          = aic_val, 
    CFI          = cfi_val, 
    SRMR         = srmr_val, 
    p_chi2       = p_chi2,
    SAR_R2       = r2_sar, 
    Latent_Beta  = path_beta, 
    Latent_p     = path_p,
    worst_ind    = loadings$rhs[1], 
    worst_p      = round(loadings$pvalue[1], 4)
  )
  
  # Stop if fewer than 2 indicators left
  if (length(current_indicators) <= 2) {
    cat("\nReached minimum 2 indicators for latent variable. Pruning stops.\n")
    break
  }
  
  # Check for non-significant loadings (p >= 0.05)
  worst_ind <- loadings$rhs[1]
  worst_p   <- loadings$pvalue[1]
  
  if (is.na(worst_p) || worst_p < 0.05) {
    cat(sprintf("\n[SUCCESS] All remaining factor loadings are statistically significant (p < 0.05)!\n"))
    break
  }
  
  cat(sprintf("\n>>> Dropping least significant loading: '%s' (p = %.4f)\n", worst_ind, worst_p))
  current_indicators <- setdiff(current_indicators, worst_ind)
  step <- step + 1
}

# ------------------------------------------------------------------------------
# 3. SUMMARY OF PRUNING TRAJECTORY
# ------------------------------------------------------------------------------
prune_summary <- bind_rows(step_history)
cat("\n====================================================================\n")
cat("PRUNING TRAJECTORY SUMMARY\n")
cat("====================================================================\n")
print(as.data.frame(prune_summary))

out_summary_path <- file.path(data_dir, "obj4_explore_EGOA_lower_pruning_history.csv")
write_csv(prune_summary, out_summary_path)
cat(sprintf("\nSaved pruning history to: %s\n", out_summary_path))