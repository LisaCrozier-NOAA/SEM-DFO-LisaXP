# ==============================================================================
# Script: 11_explore_competition_sweep_results.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Load and rank competition sweep results across NCC/AK and 24yr/19yr
# ==============================================================================

library(tidyverse)

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
data_out_dir <- file.path(proj_dir, "Routput_for_paper", "data")

# ------------------------------------------------------------------------------
# 1. LOAD DATASETS
# ------------------------------------------------------------------------------
load_and_rank <- function(region, years) {
  f_path <- file.path(data_out_dir, paste0("Goal3_Competition_Sweep_", region, "_", years, ".csv"))
  if (!file.exists(f_path)) {
    warning("File not found: ", f_path)
    return(tibble())
  }
  read_csv(f_path, show_col_types = FALSE) %>%
    arrange(abs_final_aic)
}

ncc24 <- load_and_rank("NCC", "24yr")
ak24  <- load_and_rank("AK",  "24yr")
ncc19 <- load_and_rank("NCC", "19yr")
ak19  <- load_and_rank("AK",  "19yr")

# ------------------------------------------------------------------------------
# 2. QUICK INSPECTION SUMMARY
# ------------------------------------------------------------------------------
# Key columns for diagnostic evaluation
display_cols <- c(
  "predator", "base_aic", "abs_final_aic", "abs_delta_aic",
  "base_r2", "abs_final_r2", "base_beta", "abs_beta_peff",
  "abs_variant_type", "abs_prey", "Q1_aic_improved", "Q3_sign_flipped"
)

cat("\n=== TOP NCC (24-Year) ===\n")
ncc24 %>% select(any_of(display_cols)) %>% print(n = 15)

cat("\n=== TOP NCC (24-Year) ===\n")
ncc24 %>% select(any_of(display_cols)) %>% print(n = 15)

cat("\n=== TOP AK (24-Year) ===\n")
ak24 %>% select(any_of(display_cols)) %>% print(n = 15)

cat("\n=== TOP NCC (19-Year) ===\n")
ncc19 %>% select(any_of(display_cols)) %>% print(n = 15)

cat("\n=== TOP AK (19-Year) ===\n")
ak19 %>% select(any_of(display_cols)) %>% print(n = 15)


#Lisa's summary -- 
ncc24 #hake alone still wins!!!!! not only were aic for the best other models higher, but the coef were still positive for all but cormorants, and they had zero explanatory value
ncc19 #hake alone wins again! technically loons+2prey had daic of 1, but since that doesn't count the 2 extra parameters, hake wins
ak24  #no hands down winner by aic, all w/in daic of 7, and 9 preds had daic<3
      #

#there was never a case when I would choose neg model over abs model
#there are 2 predators I can eliminate because even w/ altprey, peff is still positive -- X10_harbour_seal_DFO_WCVI_smoltyr & X15_pollockGOA_SSB_lag1_smoltyr
#7 candidates remain with similar performance, with the strongest impact factor and swapped signs coming from sablefish3plus
#sleeper sharks, harbor seals (raw or sst), sablefish3+, seak ssl (altprey), and maybe a suite of predatory fish (+/-)
ak24 %>% filter(abs_pvalue<0.05) %>% select(predator,base_aic,abs_final_aic,neg_final_aic,base_r2,abs_final_r2,neg_final_r2,base_beta,abs_beta_peff) 
ak24 %>% filter(abs_beta_peff<0,abs_pvalue<0.05) %>% select(predator,base_aic,abs_final_aic,base_beta,abs_beta_peff,base_r2,abs_final_r2,abs_variant_type,abs_prey) 
ak24 %>% filter(abs_beta_peff<0,abs_pvalue<0.05,Q4_impact_increased==TRUE) %>% select(predator,base_aic,abs_final_aic,base_beta,abs_beta_peff,abs_beta_std,base_r2,abs_final_r2,abs_variant_type,abs_prey) 

ak24 %>% filter(predator=="X11_ssl_seak_pup_pred_smoltyr"|predator=="X15_DFA_halibutBiomassAge8plus_predAK_adultyr") %>% select(predator,abs_variant_type,abs_prey,abs_final_r2,abs_final_aic)

# predator                                      base_aic abs_final_aic base_beta abs_beta_peff base_r2 abs_final_r2 abs_variant_type       abs_prey                                                   
# 1 X15_DFA_sleeperSharkBSAI_predAK_smoltyr          -26.4         -29.8    -0.246        -0.593   0.481        0.556 Joint Temp + 1 AltPrey capelin_WGoA_smoltyr                                       
# 2 X10_Harbour_s_WS_adultyr                         -29.0         -29.0    -0.530        -0.530   0.516        0.516 Raw Baseline           None                                                       
# 3 X15_sablefishAK_biomass4plus_smoltyr             -25.6         -28.8     0.420        -2.03    0.466        0.540 Joint Temp + 1 AltPrey sitkaHerring_EGoA_smoltyr                                  
# 4 X15_sablefishAK_biomass3plus_smoltyr             -27.6         -28.8     0.504        -2.41    0.519        0.534 Joint Temp + 2 AltPrey X13_pollockGOA_age1and2_smoltyr + sitkaHerring_EGoA_smoltyr
# 5 X10_harbour_seal_DFO_Coastwide_smoltyr           -27.3         -28.7    -0.314        -0.486   0.506        0.534 Temp Only              None                                                       
# 6 X11_ssl_seak_pup_pred_smoltyr                    -23.1         -27.1     0.122        -0.583   0.393        0.491 Joint Temp + 2 AltPrey sitkaHerring_EGoA_smoltyr + capelin_WGoA_smoltyr           
# 7 X15_DFA_halibutBiomassAge8plus_predAK_adultyr    -25.1         -26.9    -0.287        -0.509   0.457        0.500 Joint Temp + 1 AltPrey capelin_WGoA_adultyr 

ak19 %>% filter(abs_beta_peff<0,abs_pvalue<0.05,Q4_impact_increased==TRUE) %>% select(predator,base_aic,abs_final_aic,base_beta,abs_beta_peff,base_r2,abs_final_r2,abs_variant_type,abs_prey) 
ak19 %>% filter(abs_beta_peff<0,abs_pvalue<0.05,Q4_impact_increased==TRUE) %>% select(predator,base_aic,abs_final_aic,base_beta,abs_beta_peff,base_r2,abs_final_r2,abs_variant_type,abs_prey) 
# A tibble: 9 × 9
# predator                                      base_aic abs_final_aic base_beta abs_beta_peff base_r2 abs_final_r2 abs_variant_type       abs_prey                                                   
# 1 X15_DFA_sleeperSharkBSAI_predAK_adultyr          -32.4         -33.8    -0.439        -0.662   0.816        0.832 AltPrey Only (2 Prey)  sitkaHerring_EGoA_adultyr + capelin_WGoA_adultyr           
# 2 X15_spinyDogfishBSAI_predAK_smoltyr              -23.4         -33.3     0.162        -1.34    0.659        0.819 Joint Temp + 2 AltPrey capelin_WGoA_smoltyr + X12_egoa_krill_smoltyr              
# 3 X15_DFA_sleeperSharkBSAI_predAK_smoltyr          -31.6         -33.0    -0.314        -0.494   0.808        0.822 Joint Temp + 2 AltPrey X13_pollockGOA_age1and2_smoltyr + capelin_WGoA_smoltyr     
# 4 X10_Harbour_s_WS_adultyr                         -30.7         -31.5    -0.462        -0.580   0.792        0.803 AltPrey Only (1 Prey)  X05_eulachon_during_chinook_adultyr                        
# 5 X15_DFA_halibutBiomassAge8plus_predAK_adultyr    -25.0         -31.5    -0.222        -0.418   0.709        0.822 AltPrey Only (1 Prey)  capelin_WGoA_adultyr                                       
# 6 X15_sablefishAK_biomass4plus_smoltyr             -24.9         -30.5    -0.193        -0.703   0.707        0.786 Joint Temp + 2 AltPrey sitkaHerring_EGoA_smoltyr + capelin_WGoA_smoltyr           
# 7 X10_harbour_seal_DFO_Coastwide_smoltyr           -27.4         -29.6    -0.289        -0.420   0.754        0.782 Temp Only              None                                                       
# 8 X15_sablefishAK_biomass3plus_smoltyr             -23.2         -29.1    -0.118        -1.07    0.664        0.765 Joint Temp + 2 AltPrey X13_pollockGOA_age1and2_smoltyr + sitkaHerring_EGoA_smoltyr
# 9 X15_salmonSharkBSAI_predAK_smoltyr               -26.1         -28.6    -0.362        -1.12    0.721        0.756 Joint Temp + 2 AltPrey sitkaHerring_EGoA_smoltyr + capelin_WGoA_smoltyr           

names(ak24)

ak24 %>% select(predator,abs_variant_type,abs_prey,abs_final_r2,abs_final_aic) %>% filter(grepl("X15_DFA_sleeperSharkBSAI_predAK_smoltyr|sablefishAK_biomass4plus|X11_ssl_seak_pup", predator))

                                                                                 