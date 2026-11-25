
# ==============================================================================
# Script: 02_analysis/06e_explore_free_param_core_model.R
# ==============================================================================

#Missing variables: (plus NHL, BSAI, murreJSOES)
        #                                     column n_valid n_missing first_year last_year 
        # X01_DFA_sumPreyOfPrey_planktonJun_smoltyr      23         1       1999      2021
        #                X06_Lmu_IntSprMayW_smoltyr      23         1       1999      2021
        #                  X06_StomFull_May_smoltyr      23         1       1999      2021
        #   X15_DFA_sleeperSharkBSAI_predAK_adultyr      23         1       1998      2020
        #         X15_salmonSharkGoA_predAK_adultyr      23         1       1998      2020
        #            X10_ssl_est_wholerange_adultyr      22         2       1998      2019
 

# ==============================================================================
# Step 2: Reload Master Dataset, Filter (<= 2 Missing Yrs), & Run Backbone 1B
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

proj_dir    <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
data_dir    <- file.path(proj_dir, "Outputs/data")
master_path <- file.path(data_dir, "sem_master_peff_zscored.csv")

# ------------------------------------------------------------------------------
# 1. RELOAD & FILTER TO TIME SERIES WITH <= 2 MISSING YEARS (1998-2024)
# ------------------------------------------------------------------------------
sem_master <- read_csv(master_path, show_col_types = FALSE) %>%
  filter(year >= 1998 & year <= 2024)

# Keep year + columns with at most 2 missing values in the 1998-2024 window
valid_cols <- sem_master %>%
  select(-year) %>%
  keep(~ sum(is.na(.x)) <= 2) %>%
  names()

sem_data <- sem_master %>%
  select(year, all_of(valid_cols))

cat("Years in sem_data:", paste(range(sem_data$year), collapse = " - "),
    "(n =", nrow(sem_data), "years)\n")
cat("Surviving columns (<= 2 missing years):", length(valid_cols), "\n\n")
print(valid_cols)

# ------------------------------------------------------------------------------
# 2. DEFINE & FIT BACKBONE MODEL "1B" (PURELY STRUCTURAL)
#    Central block: sar ~ growth + abundance
#    Outer DFA/composite nodes:
#      - preyNCC   -> growth    (X05_DFA_abundSardine_smoltyr)
#      - predNCC   -> abundance (X09_DFA_HakeAge5Plus_smoltyr)
#      - preyAK    -> sar       (X12_biomassDecaShelfSum_smoltyr)
#      - predAK    -> sar       (sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z)
# ------------------------------------------------------------------------------
#X05_DFA_abundSardine_smoltyr +  X12_biomassDecaShelfSum_smoltyr makes growth->sar NS, and cpue->SAR NS
#I'm using the top DAG1B_short model, because we are allowing <2 missing years
#it uses seanettle (02.CompJellyNCC_DFA1) and EGOA lower trophic, which are both DFAs (12.ZooPreyAK_DFA1)
model_1B <- "
  # Central salmon regressions
  X06_DFA_IGF_mu_smoltyr       ~ X02_DFA_SeaNettle_smoltyr
  X07_cpue_IntSprJunHW_smoltyr ~ X09_DFA_HakeAge5Plus_smoltyr
  X16_DFA_sar_UCsprCh_smoltyr  ~ X06_DFA_IGF_mu_smoltyr + 
                                 X07_cpue_IntSprJunHW_smoltyr + 
                                 X12_DFA_EGOA_DFA_lowerTrophic_smoltyr + 
                                 sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z
                                 
  # Exogenous covariances among all 4 prey/predator drivers (6 parameters)
  X02_DFA_SeaNettle_smoltyr ~~ X09_DFA_HakeAge5Plus_smoltyr + 
                               X12_DFA_EGOA_DFA_lowerTrophic_smoltyr + 
                               sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z
  X09_DFA_HakeAge5Plus_smoltyr ~~ X12_DFA_EGOA_DFA_lowerTrophic_smoltyr + 
                                  sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z
  X12_DFA_EGOA_DFA_lowerTrophic_smoltyr ~~ sablefishAK_biomass3plus_smoltyr__peff_TempPrey2_z
"


fit_1B <- sem(
  model_1B,
  data         = sem_data ,
  missing      = "ML",
  bounds       = "pos.var",
  fixed.x      = TRUE,
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_1B, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

cat("\nFree parameters estimated (npar):", fitMeasures(fit_1B, "npar"),
    "| Sample size (N):", lavInspect(fit_1B, "ntotal"), "\n")

