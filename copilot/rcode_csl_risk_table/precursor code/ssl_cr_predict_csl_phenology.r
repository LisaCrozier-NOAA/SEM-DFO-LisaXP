#calculate overlap with salmon after filling in CSL data based on eulachon timing


library(tidyverse)
library(mgcv)

# -----------------------------------------------------------------------------
# 1. Load Data, Filter Apr-Jun, and Calculate 2-Week Back-Tracked Chinook
# -----------------------------------------------------------------------------
week_all <- read.csv(
  paste0(jake_path, "Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.all.results.csv"),
  row.names = NULL
)

# A. Prepare Estuary Chinook by Back-Tracking Bonneville Passage by 2 Weeks
chinook_estuary_df <- week_all %>%
  select(year, week, Spring_Achin_bonn_pass) %>%
  mutate(
    # Shift week back by 2 weeks to represent arrival in lower river / estuary
    week_estuary = week - 2
  ) %>%
  filter(week_estuary >= 1) %>%
  select(year, week = week_estuary, chinook_estuary_pass = Spring_Achin_bonn_pass)

# B. Assemble Primary Dataset
week_data_prep <- week_all %>%
  filter(month %in% c(4, 5, 6), year >= 2011, year <= 2024) %>%
  select(year, week, month, date, csl_nonpup_total_emb, eulachon_ssb_4week_est) %>%
  # Join the 2-week back-tracked Chinook counts
  left_join(chinook_estuary_df, by = c("year", "week")) %>%
  mutate(
    year_factor   = factor(year),
    # Fill NAs in Chinook with 0 if outside run window
    chinook_estuary_pass = if_else(is.na(chinook_estuary_pass), 0, chinook_estuary_pass)
  )

# -----------------------------------------------------------------------------
# 2. Fit Phenological GAM: Sea Lions ~ s(Week) + s(Eulachon Biomass)
# -----------------------------------------------------------------------------
# Predicts CSL weekly presence based on seasonal timing and Eulachon availability
csl_eul_gam <- gam(
  csl_nonpup_total_emb ~ s(week, k = 5) + 
    s(eulachon_ssb_4week_est, k = 5) + 
    s(year_factor, bs = "re"),
  data = week_data_prep %>% filter(!is.na(csl_nonpup_total_emb) & !is.na(eulachon_ssb_4week_est)),
  family = quasipoisson(link = "log"),
  method = "REML"
)

# -----------------------------------------------------------------------------
# 3. Predict CSL Phenology & Calculate Chinook x CSL Overlap
# -----------------------------------------------------------------------------
week_reconstructed <- week_data_prep

# Predict CSL abundance across all weeks (including missing weeks and 2019)
week_reconstructed$csl_gam_pred <- predict(
  csl_eul_gam,
  newdata = week_reconstructed,
  type = "response",
  exclude = "s(year_factor)" # Population-average prediction for missing year levels
)

week_final_overlap <- week_reconstructed %>%
  mutate(
    # Use observed CSL counts if present, otherwise use Eulachon-driven prediction
    csl_final = if_else(is.na(csl_nonpup_total_emb), csl_gam_pred, as.numeric(csl_nonpup_total_emb)),
    csl_final = pmax(0, csl_final),
    
    # WEEKLY OVERLAP PRODUCT: Reconstructed Sea Lions x 2-Week Back-Tracked Chinook
    weekly_csl_chinook_overlap = csl_final * chinook_estuary_pass
  )

# -----------------------------------------------------------------------------
# 4. Summary: Annual Integrated Salmon-Sea Lion Overlap Index
# -----------------------------------------------------------------------------
safe_sum <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)

salmon_csl_annual_index <- week_final_overlap %>%
  group_by(year) %>%
  summarise(
    weeks_eval             = n(),
    csl_obs_weeks          = sum(!is.na(csl_nonpup_total_emb)),
    
    # Total Seasonal Densities
    total_estuary_chinook  = safe_sum(chinook_estuary_pass),
    total_csl_exposure     = safe_sum(csl_final),
    
    # ANNUAL INTEGRATED SALMON-SEA LION OVERLAP INDEX
    i_csl_chinook_overlap  = safe_sum(weekly_csl_chinook_overlap)
  )

cat("--- ANNUAL SALMON x CSL PREDATION OVERLAP INDEX ---\n")
print(salmon_csl_annual_index, n = 20)
