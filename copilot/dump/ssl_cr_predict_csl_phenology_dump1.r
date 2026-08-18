#Fill in missing phenology for 2018-2023

library(tidyverse)
library(zoo)

# -----------------------------------------------------------------------------
# 1. Read in the Complete Weekly Dataset (week_all)
# -----------------------------------------------------------------------------
week_all <- read.csv(
  paste0(jake_path, "Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.all.results.csv"),
  row.names = NULL
)

# Quick check of structure
head(week_all)

# -----------------------------------------------------------------------------
# 2. Filter Apr-Jun, Interpolate CSL & Eulachon Missing Weeks
# -----------------------------------------------------------------------------
week_eulachon_imputed <- week_all %>%
  # Filter to Apr-Jun (months 4, 5, 6) across target years
  filter(month %in% c(4, 5, 6), year >= 2011, year <= 2024) %>%
  arrange(year, week) %>%
  group_by(year) %>%
  mutate(
    # Check data availability count per year
    n_csl  = sum(!is.na(csl_nonpup_total_emb)),
    n_eul  = sum(!is.na(eulachon_ssb_4week_est)),
    
    # A. Within-year linear interpolation for Sea Lions (if >= 4 weeks present)
    csl_interp = if_else(
      n_csl >= 4,
      zoo::na.approx(csl_nonpup_total_emb, x = week, na.rm = FALSE, rule = 2),
      csl_nonpup_total_emb # Leaves entirely missing years (2019) as NA
    ),
    
    # B. Within-year linear interpolation for Eulachon (if >= 4 weeks present)
    eul_interp = if_else(
      n_eul >= 4,
      zoo::na.approx(eulachon_ssb_4week_est, x = week, na.rm = FALSE, rule = 2),
      eulachon_ssb_4week_est
    ),
    
    # Enforce non-negativity constraints
    csl_interp = pmax(0, csl_interp),
    eul_interp = pmax(0, eul_interp),
    
    # C. Calculate weekly interaction term
    weekly_csl_eul_interaction = csl_interp * eul_interp
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 3. Calculate Annual Overlap Summary Metrics
# -----------------------------------------------------------------------------
# Safe wrapper functions for max and sum over entirely NA groups
safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_sum <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  sum(x, na.rm = TRUE)
}

# Updated Annual Overlap Summary without Warnings
eulachon_csl_annual_overlap <- week_eulachon_imputed %>%
  group_by(year) %>%
  summarise(
    weeks_evaluated = n(),
    csl_weeks_obs   = sum(!is.na(csl_nonpup_total_emb)),
    eul_weeks_obs   = sum(!is.na(eulachon_ssb_4week_est)),
    
    # Annual Integrated Predation Overlap Index (Sum of Weekly Products)
    i_csl_eulachon_overlap = safe_sum(weekly_csl_eul_interaction),
    
    # Peak Seasonal Densities
    csl_peak_apr_jun = safe_max(csl_interp),
    eul_peak_apr_jun = safe_max(eul_interp)
  )


# Display annual index table
cat("--- ANNUAL EULACHON x CSL OVERLAP INDEX SUMMARY ---\n")
print(eulachon_csl_annual_overlap, n = 20)