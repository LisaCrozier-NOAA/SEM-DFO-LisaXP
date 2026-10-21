library(tidyverse)
library(lavaan)

# 1. Define Variables
target_cpue <- "X07_DFA_cpue_IntSprJunHW"
pred_col    <- "X09_DFA_HakeAge5Plus_smoltyr"
hci_col     <- "X12_DFA_biomassEuphShelfSum_smoltyr"
sst_col     <- "X01_habCompInd_smoltyr"
prey1_col   <- "X04_marketsquid_GAM_smoltyr"
prey2_col   <- "X05_anchovy_GAM_smoltyr"

# 2. Extract Vectors from minmax_complete_data
P_norm <- minmax_complete_data[[pred_col]]
HCI    <- minmax_complete_data[[hci_col]]
SST    <- minmax_complete_data[[sst_col]]
A1     <- minmax_complete_data[[prey1_col]]
A2     <- minmax_complete_data[[prey2_col]]

# 3. Fit Raw Unattenuated Baseline (CPUE ~ Raw Hake)
fit_base <- sem(paste0(target_cpue, " ~ ", pred_col), data = minmax_complete_data)
summary(fit_base, fit.measures = TRUE, rsquare = TRUE)

# 4. Construct Peff with Specific Dampening Weights (e.g., w_krill=0.5, w_squid=0.5, w_anchovy=0.5)
# (Adjust weights below as needed)
w_krill <- 1
w_prey1 <- 1
w_prey2 <- 1
w_sst <- 1
# w_krill <- 0.5
# w_prey1 <- 0.5
# w_prey2 <- 0.5

df_single <- minmax_complete_data %>%
  mutate(P_eff_0 = P_norm,
         P_eff_krill = P_norm * (1 - w_krill * HCI),
         P_eff_hci = P_norm * (1 - w_krill * SST),
         P_eff_squid = P_norm * (1 - w_krill * A1))
# * (1 - w_prey1 * A1) #* (1 - w_sst * SST))

# 5. Fit Single Attenuated Model (CPUE ~ Peff)
fit_peff_0 <- sem(paste0(target_cpue, " ~ P_eff_0"), data = df_single)
fit_peff_krill <- sem(paste0(target_cpue, " ~ P_eff_krill"), data = df_single)
fit_peff_sst <- sem(paste0(target_cpue, " ~ P_eff_hci"), data = df_single)
fit_peff_squid <- sem(paste0(target_cpue, " ~ P_eff_squid"), data = df_single)
AIC(fit_peff_0,fit_peff_krill,fit_peff_sst,fit_peff_squid)
sem_output_summary_fxn(fit_peff_0)
sem_output_summary_fxn(fit_peff_krill)
sem_output_summary_fxn(fit_peff_sst)
sem_output_summary_fxn(fit_peff_squid)
sem_output_summary_fxn(fit_peff_0)

# df         AIC
# fit_peff_0      2 -11.1396444
# fit_peff_krill  2  -4.7712639
# fit_peff_sst    2  -5.9541087
# fit_peff_squid  2  -0.8931579
# 
