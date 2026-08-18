#Fill in missing phenology for 2018-2023 using a seasonal cycle plus the presence of eulachon


library(tidyverse)
library(mgcv)

# -----------------------------------------------------------------------------
# 1. Prepare Data for Dynamic Phenology GAM
# -----------------------------------------------------------------------------
csl_eul_data <- week_all %>%
  filter(month %in% c(4, 5, 6), year >= 2011, year <= 2024) %>%
  mutate(
    year_factor = factor(year),
    # Ensure Eulachon is non-NA for fitting (use the interpolated/spline eulachon if available)
    eulachon_input = eulachon_ssb_4week_est
  )

# -----------------------------------------------------------------------------
# 2. Fit GAM Driven by Weekly Eulachon Biomass
# -----------------------------------------------------------------------------
# csl ~ s(week) + s(eulachon) allows sea lion abundance to scale dynamically 
# up or down depending on how much Eulachon is in the river that specific week.

csl_eul_gam <- gam(
  csl_nonpup_total_emb ~ s(week, k = 5) + 
    s(eulachon_input, k = 5) + 
    s(year_factor, bs = "re"), # Random effect for annual baseline variations
  data = csl_eul_data %>% filter(!is.na(csl_nonpup_total_emb) & !is.na(eulachon_input)),
  family = quasipoisson(link = "log"),
  method = "REML"
)

# Inspect model summary—check if s(eulachon_input) is statistically significant!
summary(csl_eul_gam)

# -----------------------------------------------------------------------------
# 3. Predict Missing Sea Lion Weeks (Including 2019) Driven by Eulachon
# -----------------------------------------------------------------------------
# For 2019, because 'year_factor' was never observed, we exclude the random 
# effect (exclude = "s(year_factor)") so the prediction relies directly on 
# week + Eulachon biomass!

week_eul_predicted <- csl_eul_data

week_eul_predicted$csl_gam_pred <- predict(
  csl_eul_gam,
  newdata = week_eul_predicted,
  type = "response",
  exclude = "s(year_factor)" # Uses population-average intercept + Eulachon drive for 2019
)

# -----------------------------------------------------------------------------
# 4. Construct Final Reconstructed Dataset & Overlap Index
# -----------------------------------------------------------------------------
week_final_eul_driven <- week_eul_predicted %>%
  mutate(
    # Use observed sea lions if present; otherwise use Eulachon-driven GAM prediction
    csl_final = if_else(is.na(csl_nonpup_total_emb), csl_gam_pred, as.numeric(csl_nonpup_total_emb)),
    csl_final = pmax(0, csl_final),
    
    # Recalculate weekly interaction product
    weekly_overlap = csl_final * eulachon_input
  )

# -----------------------------------------------------------------------------
# 5. Calculate Annual Integrated Overlap Index
# -----------------------------------------------------------------------------
safe_sum <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)

pheno_eul_overlap_summary <- week_final_eul_driven %>%
  group_by(year) %>%
  summarise(
    obs_csl_weeks     = sum(!is.na(csl_nonpup_total_emb)),
    obs_eul_weeks     = sum(!is.na(eulachon_ssb_4week_est)),
    
    # Reconstructed sea lion exposure (AUC)
    total_csl_exposure = safe_sum(csl_final),
    
    # Integrated Eulachon x Sea Lion Overlap Index
    i_csl_eulachon_overlap = safe_sum(weekly_overlap)
  )

print(pheno_eul_overlap_summary, n = 20)