
# RESULTS
write.csv(annual_risk_table, file.path(output_dir, "chinook_quantile_csl_exposure_risk_1998_2024.csv"), row.names = FALSE)

print(annual_risk_table,n=Inf)
# A tibble: 27 × 11
# year csl_annual_peak csl_during_chinook csl_weekly_avg_chin total_chinook_estuary eulachon_during_chinook shad_during_chinook csl_exposure_threat eulachon_mitigation shad_mitigation overall_chinook_risk
# 1  1998           114                 582.                44.8                 51566                 161163.             1852775 Low Exposure        Moderate Buffer     Weak Late Buff… LOW                 
# 2  1999            31.5               238.                19.8                 51504                 166093.             1521323 Low Exposure        Moderate Buffer     Weak Late Buff… LOW                 
# 3  2000            44                 280.                28.0                197185                 279254.              454228 Low Exposure        Strong Buffer       Weak Late Buff… LOW                 
# 4  2001           111                 866.                78.7                437687                2043469.             1539014 Moderate Exposure   Strong Buffer       Weak Late Buff… MODERATE-LOW        
# 5  2002           190                1183.                98.6                355237                3531401.             3069988 Moderate Exposure   Strong Buffer       Strong Late Bu… MODERATE-LOW        
# 6  2003           299                1897.               136.                 265924                8295187.             4462176 Moderate Exposure   Strong Buffer       Strong Late Bu… MODERATE-LOW        
# 7  2004            89                 652                 54.3                233502                1178617.             5206649 Low Exposure        Strong Buffer       Strong Late Bu… LOW                 
# 8  2005           139                 841                 76.5                126012                  44592.             4094470 Low Exposure        Moderate Buffer     Strong Late Bu… LOW                 
# 9  2006            88                 535                 59.4                167739                  12157.             3790099 Low Exposure        Low Buffer          Strong Late Bu… LOW                 
# 10  2007            91                 804                 67                   98730                 105591.             2480456 Low Exposure        Moderate Buffer     Weak Late Buff… LOW                 
# 11  2008           134                 822                 68.5                184679                 148793.             1923455 Low Exposure        Moderate Buffer     Weak Late Buff… LOW                 
# 12  2009           120                 761                 69.2                176472                 100704.             1346660 Low Exposure        Moderate Buffer     Weak Late Buff… LOW                 
# 13  2010           181                1450                121.                 318815                      0               998134 Moderate Exposure   Low Buffer          Weak Late Buff… MODERATE            
# 14  2011           182                 899                 99.9                242662                  70238.              787066 Moderate Exposure   Moderate Buffer     Weak Late Buff… MODERATE            
# 15  2012           163                1017                102.                 215325                   2357.             2349979 Moderate Exposure   Low Buffer          Weak Late Buff… MODERATE            
# 16  2013           628                4179                418.                 146616                 146334.             3690625 Moderate Exposure   Moderate Buffer     Strong Late Bu… MODERATE            
# 17  2014          1129                4603                418.                 262705                1620949.             2577160 High Exposure       Strong Buffer       Weak Late Buff… HIGH                
# 18  2015          1971               10640                887.                 323714                  33026.             1672433 High Exposure       Low Buffer          Weak Late Buff… VERY HIGH           
# 19  2016          2484                7880.               716.                 216999                    109.             1594859 High Exposure       Low Buffer          Weak Late Buff… VERY HIGH           
# 20  2017          1019                6112                679.                 144990                      0              2702859 High Exposure       Low Buffer          Strong Late Bu… VERY HIGH           
# 21  2018           893.               5227.               523.                 133024                      0              5736863 High Exposure       Low Buffer          Strong Late Bu… VERY HIGH           
# 22  2019           900.               5298.               530.                  90114                   2269.             7124495 High Exposure       Low Buffer          Strong Late Bu… VERY HIGH           
# 23  2020           907.               5197.               520.                 115527                      0              5441196 High Exposure       Low Buffer          Strong Late Bu… VERY HIGH           
# 24  2021           677                4464.               406.                 119514                  81831.             5326254 High Exposure       Moderate Buffer     Strong Late Bu… HIGH                
# 25  2022           580                2725.               272.                 221394                 943052.             5565838 Moderate Exposure   Strong Buffer       Strong Late Bu… MODERATE-LOW        
# 26  2023          1362                5013.               501.                 167641                6221179.             4144253 High Exposure       Strong Buffer       Strong Late Bu… HIGH                
# 27  2024           538                3114.               283.                 136095                 175042.             2811783 Moderate Exposure   Strong Buffer       Strong Late Bu… MODERATE-LOW   



# ==============================================================================
# Script: Absolute CSL Exposure & Predation Risk Assessment
# Replaces Pearson correlation with absolute CSL counts during Chinook run
# ==============================================================================

library(tidyverse)

output_dir <- "copilot/outputs_csl_cr"

# Helper function to avoid -Inf warnings on all-NA years
safe_max <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(NA_real_)
  max(x_clean)
}

# -----------------------------------------------------------------------------
# 1. Load Master Dataset & Apply 2-Week Estuary Shift
# -----------------------------------------------------------------------------

master_all <- read.csv(file.path(output_dir, "master_weekly_all_species_1998_2024.csv"))

master_estuary <- master_all %>%
  select(year, week, csl_final, eulachon_final, Chin_weekly_count, Shad_weekly_count) %>%
  group_by(year) %>%
  mutate(
    # 2-Week Shift for Estuary Arrival Timing
    chin_estuary_count = lead(Chin_weekly_count, 2, default = 0),
    shad_estuary_count = lead(Shad_weekly_count, 2, default = 0)
  ) %>%
  # Define within-year proportions to identify main run timing
  mutate(
    chin_max  = safe_max(chin_estuary_count),
    chin_prop = if_else(is.na(chin_max) | chin_max <= 0, 0, chin_estuary_count / chin_max),
    
    # Flag weeks where Chinook presence is >= 10% of seasonal peak (Active Run Window)
    is_chinook_active = chin_prop >= 0.10
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 2. Compute Absolute CSL Exposure During Active Chinook Window
# -----------------------------------------------------------------------------

annual_risk_table <- master_estuary %>%
  group_by(year) %>%
  summarise(
    # 1. Primary Risk Drivers
    csl_annual_peak      = safe_max(csl_final),
    
    # CSL counts integrated specifically across weeks when Chinook are present (>=10% peak)
    csl_during_chinook   = sum(csl_final[is_chinook_active], na.rm = TRUE),
    csl_weekly_avg_chin  = mean(csl_final[is_chinook_active], na.rm = TRUE),
    
    # 2. Total Chinook Abundance in Estuary
    total_chinook_estuary = sum(chin_estuary_count, na.rm = TRUE),
    
    # 3. Alternate Prey Buffers Available DURING Active Chinook Window
    eulachon_during_chinook = sum(eulachon_final[is_chinook_active], na.rm = TRUE),
    shad_during_chinook     = sum(shad_estuary_count[is_chinook_active], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Clean up NAs/zero handling
    csl_during_chinook = replace_na(csl_during_chinook, 0),
    eulachon_during_chinook = replace_na(eulachon_during_chinook, 0),
    shad_during_chinook = replace_na(shad_during_chinook, 0),
    
    # -------------------------------------------------------------------------
    # Categorize Risk Tiers Based on Absolute Exposure & Buffers
    # -------------------------------------------------------------------------
    
    # Tier 1: CSL Exposure Threat (Absolute CSL count during run)
    csl_exposure_threat = case_when(
      csl_during_chinook >= quantile(csl_during_chinook, 0.66) ~ "High Exposure",
      csl_during_chinook >= quantile(csl_during_chinook, 0.33) ~ "Moderate Exposure",
      TRUE ~ "Low Exposure"
    ),
    
    # Tier 2: Early Spring Eulachon Buffer Strength
    eulachon_mitigation = case_when(
      eulachon_during_chinook >= quantile(eulachon_during_chinook, 0.66) ~ "Strong Buffer",
      eulachon_during_chinook >= quantile(eulachon_during_chinook, 0.33) ~ "Moderate Buffer",
      TRUE ~ "Low Buffer"
    ),
    
    # Tier 3: Late Spring Shad Buffer Strength
    shad_mitigation = case_when(
      shad_during_chinook >= quantile(shad_during_chinook, 0.50) ~ "Strong Late Buffer",
      TRUE ~ "Weak Late Buffer"
    ),
    
    # Composite Qualitative Risk Ranking
    overall_chinook_risk = case_when(
      csl_exposure_threat == "High Exposure" & eulachon_mitigation == "Low Buffer" ~ "VERY HIGH",
      csl_exposure_threat == "High Exposure" ~ "HIGH",
      csl_exposure_threat == "Moderate Exposure" & eulachon_mitigation == "Strong Buffer" ~ "MODERATE-LOW",
      csl_exposure_threat == "Moderate Exposure" ~ "MODERATE",
      TRUE ~ "LOW"
    )
  )

# -----------------------------------------------------------------------------
# 3. Print & Save Results
# -----------------------------------------------------------------------------

cat("--- Absolute CSL Exposure & Risk Assessment Table (1998–2024) ---\n")
print(
  annual_risk_table %>% 
    select(year, csl_annual_peak, csl_during_chinook, eulachon_mitigation, shad_mitigation, overall_chinook_risk) %>% 
    print(n = 30)
)

write.csv(annual_risk_table, file.path(output_dir, "chinook_quantile_csl_exposure_risk_1998_2024.csv"), row.names = FALSE)
