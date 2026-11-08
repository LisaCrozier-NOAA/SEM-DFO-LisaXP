

#Figure out which prey to include from importance matrix,  
#then see how correlated they are with each other (id groups of corr variables)
#then combine sets of the most correlated variables into latent nodes 
#prune until everything is significant
#compare pruned models by aic and model validity
#results: pink no growth no latent is best by far

val_ranks_wide <- read_csv("Rcode_for_paper/Outputs/data/obj2_data_variable_importance_wide.csv") 
print(val_ranks_wide %>% filter(grepl("X14|X13|X12",LisaName)) %>% arrange(desc(DAG1B_long)),n=Inf)

akprey_cor<-read.csv(file.path(data_out_dir,"AKPrey_cor_matrix.csv"),row.names=1);akprey_cor
print(round(akprey_cor, 3))
#                  Pink Decapods Mysids EGOA_lower WGOA_mid WGOA_lower EGOA_seabirds Pollock_age1_2
# Pink            1.000   -0.139 -0.098     -0.450   -0.309      0.432        -0.477          0.407
# Decapods       -0.139    1.000 -0.200     -0.327    0.535     -0.570         0.198         -0.218
# Mysids         -0.098   -0.200  1.000      0.572    0.286      0.010        -0.371         -0.327
# EGOA_lower     -0.450   -0.327  0.572      1.000   -0.074      0.120        -0.191         -0.285
# WGOA_mid       -0.309    0.535  0.286     -0.074    1.000     -0.777         0.331         -0.163
# WGOA_lower      0.432   -0.570  0.010      0.120   -0.777      1.000        -0.602         -0.015
# EGOA_seabirds  -0.477    0.198 -0.371     -0.191    0.331     -0.602         1.000          0.111
# Pollock_age1_2  0.407   -0.218 -0.327     -0.285   -0.163     -0.015         0.111          1.000

# ==============================================
#   MODEL SUMMARY REPORT               
# ==============================================
#   
#   --- Model Fit Metrics ---
#   pvalue    cfi    aic  rmsea 
# NA  1.000 54.457  0.000 
# 
# --- R-Squared (Endogenous Variables) ---
#   No target R2 variables found in this model.
# 
# --- Regression Path Ranking (By P-Value) ---
#   lhs op         rhs    est    se      z pvalue
# SAR  ~        Pink -0.288 0.128 -2.246  0.025
# SAR  ~ P_ocean_eff -0.388 0.137 -2.834  0.005
# SAR  ~        CPUE  0.506 0.137  3.681  0.000


# ==============================================================================
# Script: 02_analysis/11b_fit_dag1b_two_latent_bipolar_model.R
# Purpose: Fit DAG 1B with Growth restored, testing the two bipolar latent
#          clusters identified from the correlation matrix:
#          - Latent_WGOA: Decapods + WGOA_mid + WGOA_lower
#          - Latent_EGOA: Mysids + EGOA_lower + EGOA_seabirds + Pollock_age1_2
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(semPlot)
})

# ------------------------------------------------------------------------------
# 0. SETUP & DATA INGESTION----------
# ------------------------------------------------------------------------------
proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")
fig_dir  <- file.path(proj_dir, "Outputs/figures")

alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")

df_raw  <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)
peff_df <- read_csv(peff_ts_path, show_col_types = FALSE)

yr_col <- grep("^year$", names(df_raw), value = TRUE, ignore.case = TRUE)
df_raw$year <- if (length(yr_col) > 0) as.numeric(df_raw[[yr_col[1]]]) else as.numeric(rownames(df_raw))
df_24 <- df_raw %>% filter(year >= 1998 & year <= 2021)

find_col <- function(pattern, df) {
  matches <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
  if (length(matches) == 0) stop("Could not find column matching: ", pattern)
  matches[1]
}

# Assemble full DAG 1B dataset
df_dag1b <- tibble(
  year           = df_24$year,
  SAR            = df_24[[find_col("sar_UCsprCh", df_24)]],
  Growth         = df_24[[find_col("IGF_mu", df_24)]],
  CPUE           = df_24[[find_col("cpue_IntSprJunHW", df_24)]],
  P_ocean_eff    = peff_df$Peff_sablefishAK_biomass3plus_smoltyr,
  P_ocean_eff_sable    = peff_df$Peff_sablefishAK_biomass3plus_smoltyr,
  P_ocean_eff_ssl    = peff_df$Peff_ssl_seak_pup_pred_smoltyr,
  
  # Group 1: Western Gulf Shelf Cluster
  Decapods       = df_24[[find_col("biomassDecaShelfSum", df_24)]],
  WGOA_mid       = df_24[[find_col("WGOA.*mid", df_24)]],
  WGOA_lower     = df_24[[find_col("WGOA.*lower", df_24)]],
  
  # Group 2: Eastern Gulf Pelagic Cluster
  Mysids         = df_24[[find_col("biomassMysidShelfSum", df_24)]],
  EGOA_lower     = df_24[[find_col("EGOA.*lower", df_24)]],
  EGOA_seabirds  = df_24[[find_col("EGOA.*seabird|EGOA.*sb", df_24)]],
  Pollock_age1_2 = df_24[[find_col("pollockGOA_age1and2", df_24)]],
  
  # Group 3: Pink
  Pink         = df_24[[find_col("pinkSalmonNorthAmerica", df_24)]]
  
  )

# Standardize variables (z-scores)
df_fit <- df_dag1b %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))


# ------------------------------------------------------------------------------
# 1a. SPECIFY DAG 1B WITH ONE BIPOLAR LATENT VARIABLES -- wgoa----
# ------------------------------------------------------------------------------
syntax_dag1b_wgoa_latent <- "
  # --- Measurement Model ---
  # Group 1: WGOA Shelf Bipolar Axis (Decapods, Mid-trophic vs. Lower-trophic)
  Latent_WGOA =~ Decapods + WGOA_mid + WGOA_lower

 
  # --- Structural Model (DAG 1B Core Backbone) ---
  SAR ~ Growth + CPUE + Latent_WGOA  + P_ocean_eff 
"

fit_dag1b_wgoa <- sem(
  syntax_dag1b_wgoa_latent,
  data         = df_fit,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_dag1b_wgoa)
#prune out growth----------
syntax_dag1b_wgoa_latent <- "
  # --- Measurement Model ---
  # Group 1: WGOA Shelf Bipolar Axis (Decapods, Mid-trophic vs. Lower-trophic)
  Latent_WGOA =~ Decapods + WGOA_mid + WGOA_lower

 
  # --- Structural Model (DAG 1B Core Backbone) ---
  SAR ~  CPUE + Latent_WGOA  + P_ocean_eff 
"

fit_dag1b_wgoa <- sem(
  syntax_dag1b_wgoa_latent,
  data         = df_fit,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_dag1b_wgoa)
#all signif, but wgoa_mid has neg var

#====================
# ------------------------------------------------------------------------------
# 1b. SPECIFY DAG 1B WITH ONE BIPOLAR LATENT VARIABLES -- egoa------
# ------------------------------------------------------------------------------
syntax_dag1b_egoa_latent <- "
  # --- Measurement Model ---
  # Group 2: EGOA Pelagic / Gadid Axis (Mysids, Lower-trophic, Seabirds, Pollock)
  Latent_EGOA =~ Mysids + EGOA_lower + EGOA_seabirds + Pollock_age1_2

 
  # --- Structural Model (DAG 1B Core Backbone) ---
  SAR ~ Growth + CPUE + Latent_EGOA  + P_ocean_eff 
"

fit_dag1b_egoa <- sem(
  syntax_dag1b_egoa_latent,
  data         = df_fit,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_dag1b_egoa)
#prune out growth----------
syntax_dag1b_egoa_latent <- "
  # --- Measurement Model ---
  # Group 2: EGOA Pelagic / Gadid Axis (Mysids, Lower-trophic, Seabirds, Pollock)
  Latent_EGOA =~ Mysids + EGOA_lower + EGOA_seabirds + Pollock_age1_2

 
  # --- Structural Model (DAG 1B Core Backbone) ---
  SAR ~  CPUE + Latent_EGOA  + P_ocean_eff 
"
#prune out seabirds------
#prune out growth----------
syntax_dag1b_egoa_latent <- "
  # --- Measurement Model ---
  # Group 2: EGOA Pelagic / Gadid Axis (Mysids, Lower-trophic, Seabirds, Pollock)
  Latent_EGOA =~ Mysids + EGOA_lower +  Pollock_age1_2

 
  # --- Structural Model (DAG 1B Core Backbone) ---
  SAR ~  CPUE + Latent_EGOA  + P_ocean_eff 
"
#prune out pollock-------
syntax_dag1b_egoa_latent <- "
  # --- Measurement Model ---
  # Group 2: EGOA Pelagic / Gadid Axis (Mysids, Lower-trophic, Seabirds, Pollock)
  Latent_EGOA =~ Mysids + EGOA_lower 

 
  # --- Structural Model (DAG 1B Core Backbone) ---
  SAR ~  CPUE + Latent_EGOA  + P_ocean_eff 
"

fit_dag1b_egoa <- sem(
  syntax_dag1b_egoa_latent,
  data         = df_fit,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_dag1b_egoa)

#all signif, but wgoa_mid has neg var
sem_output_summary_fxn(fit_dag1b_egoa)
sem_output_summary_fxn(fit_dag1b_wgoa)
#Latent_EGOA   much better by aic

#add pink-----------
syntax_dag1b_egoa_latent_pink <- "
  # --- Measurement Model ---
  # Group 2: EGOA Pelagic / Gadid Axis (Mysids, Lower-trophic, Seabirds, Pollock)
  Latent_EGOA =~ Mysids + EGOA_lower 

 
  # --- Structural Model (DAG 1B Core Backbone) ---
  SAR ~  CPUE + Latent_EGOA + Pink + P_ocean_eff 
"
#prine out EGOA_lower-----
syntax_dag1b_egoa_latent_pink <- "

  # --- Structural Model (DAG 1B Core Backbone) ---
  SAR ~  CPUE + Mysids + Pink + P_ocean_eff 
"


fit_dag1b_egoa_pink <- sem(
  syntax_dag1b_egoa_latent_pink,
  data         = df_fit,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_dag1b_egoa_pink)
sem_output_summary_fxn(fit_dag1b_egoa_pink)

#all signif, no latent, but much lower aic for fit_dag1b_egoa_pink, 
#SAR ~  CPUE + Mysids + Pink + P_ocean_eff 
sem_output_summary_fxn(fit_dag1b_wgoa)
sem_output_summary_fxn(fit_dag1b_egoa)
sem_output_summary_fxn(fit_dag1b_egoa_pink)

#run3 PINK==============
print(round(akprey_cor[,1:2] %>%arrange(desc(abs(Pink))), 3))
#                   Pink Decapods
# Pink            1.000   -0.139
# EGOA_seabirds  -0.477    0.198
# EGOA_lower     -0.450   -0.327
# WGOA_lower      0.432   -0.570
# Pollock_age1_2  0.407   -0.218
# WGOA_mid       -0.309    0.535

# Decapods       -0.139    1.000
# Mysids         -0.098   -0.200

# ------------------------------------------------------------------------------
# 1c. SPECIFY DAG 1B WITH everything strongly corr w/ pink
# ------------------------------------------------------------------------------
syntax_dag1b_pink_latent <- "
  # --- Measurement Model ---
  # Group 2: EGOA Pelagic / Gadid Axis (Mysids, Lower-trophic, Seabirds, Pollock)
  Latent_EGOA =~ Pink   +WGOA_lower  #Pollock_age1_2#+ EGOA_lower+ EGOA_seabirds+  WGOA_mid

 
  # --- Structural Model (DAG 1B Core Backbone) ---
  SAR ~  CPUE + Latent_EGOA  + P_ocean_eff 
"
syntax_dag1b_pink_latent <- "
  # --- Structural Model (DAG 1B Core Backbone) ---
  SAR ~   CPUE + Pink  + P_ocean_eff 
"

fit_dag1b_pink <- sem(
  syntax_dag1b_pink_latent,
  data         = df_fit,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

#prune  EGOA_lower, replace w/ pollock ------
#prune  pollock, replace w/ WGOA_mid ------
#prune  growth ------
#all sign, but neg var
#prune  seabirds ------
#prune  mid ------
#prune  latent, just leave pink  ------
#tried putting growth back in, but it is still highly insign

summary(fit_dag1b_pink)
