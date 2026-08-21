
#RESULTS
# Save Risk Table
write.csv(annual_absolute_risk, file.path(output_dir, "chinook_numeric_csl_exposure_risk_1998_2024.csv"), row.names = FALSE)
out_dir <- "copilot/outputs_6"
write.csv(annual_absolute_risk, file.path(out_dir, "chinook_numeric_csl_exposure_risk_1998_2024.csv"), row.names = FALSE)

# year csl_during_chinook eulachon_during_chinook shad_during_chinook csl_threat_tier               eulachon_buffer_tier      overall_chinook_risk
# 1  1998               582.                 161163.             1852775 Low CSL Exposure (<1k)        Substantial Buffer (>50k) LOW RISK            
# 2  1999               238.                 166093.             1521323 Low CSL Exposure (<1k)        Substantial Buffer (>50k) LOW RISK            
# 3  2000               280.                 279254.              454228 Low CSL Exposure (<1k)        Substantial Buffer (>50k) LOW RISK            
# 4  2001               866.                2043469.             1539014 Low CSL Exposure (<1k)        Substantial Buffer (>50k) LOW RISK            
# 5  2002              1183.                3531401.             3069988 Moderate CSL Exposure (1k-3k) Substantial Buffer (>50k) LOW-MODERATE RISK   
# 6  2003              1897.                8295187.             4462176 Moderate CSL Exposure (1k-3k) Substantial Buffer (>50k) LOW-MODERATE RISK   
# 7  2004               652                 1178617.             5206649 Low CSL Exposure (<1k)        Substantial Buffer (>50k) LOW RISK            
# 8  2005               841                   44592.             4094470 Low CSL Exposure (<1k)        Moderate Buffer (10k-50k) LOW RISK            
# 9  2006               535                   12157.             3790099 Low CSL Exposure (<1k)        Moderate Buffer (10k-50k) LOW RISK            
# 10  2007               804                  105591.             2480456 Low CSL Exposure (<1k)        Substantial Buffer (>50k) LOW RISK            
# 11  2008               822                  148793.             1923455 Low CSL Exposure (<1k)        Substantial Buffer (>50k) LOW RISK            
# 12  2009               761                  100704.             1346660 Low CSL Exposure (<1k)        Substantial Buffer (>50k) LOW RISK            
# 13  2010              1450                       0               998134 Moderate CSL Exposure (1k-3k) Minimal Buffer (<10k)     MODERATE RISK       
# 14  2011               899                   70238.              787066 Low CSL Exposure (<1k)        Substantial Buffer (>50k) LOW RISK            
# 15  2012              1017                    2357.             2349979 Moderate CSL Exposure (1k-3k) Minimal Buffer (<10k)     MODERATE RISK       
# 16  2013              4179                  146334.             3690625 High CSL Exposure (>3k)       Substantial Buffer (>50k) HIGH RISK           
# 17  2014              4603                 1620949.             2577160 High CSL Exposure (>3k)       Substantial Buffer (>50k) HIGH RISK           
# 18  2015             10640                   33026.             1672433 High CSL Exposure (>3k)       Moderate Buffer (10k-50k) HIGH RISK           
# 19  2016              7880.                    109.             1594859 High CSL Exposure (>3k)       Minimal Buffer (<10k)     VERY HIGH RISK      
# 20  2017              6112                       0              2702859 High CSL Exposure (>3k)       Minimal Buffer (<10k)     VERY HIGH RISK      
# 21  2018              5227.                      0              5736863 High CSL Exposure (>3k)       Minimal Buffer (<10k)     VERY HIGH RISK      
# 22  2019              5298.                   2269.             7124495 High CSL Exposure (>3k)       Minimal Buffer (<10k)     VERY HIGH RISK      
# 23  2020              5197.                      0              5441196 High CSL Exposure (>3k)       Minimal Buffer (<10k)     VERY HIGH RISK      
# 24  2021              4464.                  81831.             5326254 High CSL Exposure (>3k)       Substantial Buffer (>50k) HIGH RISK           
# 25  2022              2725.                 943052.             5565838 Moderate CSL Exposure (1k-3k) Substantial Buffer (>50k) LOW-MODERATE RISK   
# 26  2023              5013.                6221179.             4144253 High CSL Exposure (>3k)       Substantial Buffer (>50k) HIGH RISK           
# 27  2024              3114.                 175042.             2811783 High CSL Exposure (>3k)       Substantial Buffer (>50k) HIGH RISK    




# ==============================================================================
# Script: Absolute CSL Exposure & Risk Assessment (Estuary 2-Week Lagged)
# Uses clean master_estuary dataset with absolute counts during Chinook window
# ==============================================================================

library(tidyverse)

output_dir <- "copilot/outputs_csl_cr"

# Safe max helper function
safe_max <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(NA_real_)
  max(x_clean, na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# 1. Read Clean Estuary-Shifted Dataset
# -----------------------------------------------------------------------------

# Reads the file saved at the end of your master processing script
master_estuary <- read.csv(file.path(output_dir, "master_estuary_2wklagweekly_all_species_1998_2024.csv"))

# Ensure active Chinook flag is set (weeks >= 10% of peak estuary arrival)
master_estuary <- master_estuary %>%
  group_by(year) %>%
  mutate(
    chin_est_max  = safe_max(chin_estuary_count),
    chin_est_prop = if_else(is.na(chin_est_max) | chin_est_max <= 0, 0, chin_estuary_count / chin_est_max),
    is_chinook_active = chin_est_prop >= 0.10
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 2. Compute Absolute Exposure & Alternate Prey Volumes During Active Chinook Window
# -----------------------------------------------------------------------------

annual_absolute_risk <- master_estuary %>%
  group_by(year) %>%
  summarise(
    # Primary Threat: Absolute integrated CSL counts during active Chinook weeks
    csl_during_chinook      = sum(csl_final[is_chinook_active], na.rm = TRUE),
    csl_weekly_avg_chin     = mean(csl_final[is_chinook_active], na.rm = TRUE),
    csl_annual_peak         = safe_max(csl_final),
    
    # Total Chinook volume in estuary during active window
    chinook_estuary_total   = sum(chin_estuary_count[is_chinook_active], na.rm = TRUE),
    
    # Absolute Alternate Prey Volumes available DURING active Chinook weeks
    eulachon_during_chinook = sum(eulachon_final[is_chinook_active], na.rm = TRUE),
    shad_during_chinook     = sum(shad_estuary_count[is_chinook_active], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Clean NAs
    csl_during_chinook      = replace_na(csl_during_chinook, 0),
    eulachon_during_chinook = replace_na(eulachon_during_chinook, 0),
    shad_during_chinook     = replace_na(shad_during_chinook, 0),
    
    # -------------------------------------------------------------------------
    # Absolute Biological Magnitude Tiers
    # -------------------------------------------------------------------------
    
    # 1. Sea Lion Exposure Tier (Absolute counts during Chinook run)
    csl_threat_tier = case_when(
      csl_during_chinook >= 3000 ~ "High CSL Exposure (>3k)",
      csl_during_chinook >= 1000 ~ "Moderate CSL Exposure (1k-3k)",
      TRUE                       ~ "Low CSL Exposure (<1k)"
    ),
    
    # 2. Early Spring Eulachon Buffer Scale
    eulachon_buffer_tier = case_when(
      eulachon_during_chinook >= 50000 ~ "Substantial Buffer (>50k)",
      eulachon_during_chinook >= 10000 ~ "Moderate Buffer (10k-50k)",
      TRUE                             ~ "Minimal Buffer (<10k)"
    ),
    
    # 3. Late Spring Shad Buffer Scale
    shad_buffer_tier = case_when(
      shad_during_chinook >= 100000 ~ "Major Late Buffer (>100k)",
      shad_during_chinook >= 20000  ~ "Moderate Late Buffer (20k-100k)",
      TRUE                          ~ "Minimal Late Buffer (<20k)"
    ),
    
    # Composite Qualitative Risk Ranking
    overall_chinook_risk = case_when(
      csl_threat_tier == "High CSL Exposure (>3k)" & eulachon_buffer_tier == "Minimal Buffer (<10k)" ~ "VERY HIGH RISK",
      csl_threat_tier == "High CSL Exposure (>3k)" ~ "HIGH RISK",
      csl_threat_tier == "Moderate CSL Exposure (1k-3k)" & eulachon_buffer_tier == "Substantial Buffer (>50k)" ~ "LOW-MODERATE RISK",
      csl_threat_tier == "Moderate CSL Exposure (1k-3k)" ~ "MODERATE RISK",
      TRUE ~ "LOW RISK"
    )
  )

# -----------------------------------------------------------------------------
# 3. Display & Save Final Risk Table
# -----------------------------------------------------------------------------

cat("\n====================================================================\n")
cat("   ABSOLUTE CSL EXPOSURE & ALTERNATE PREY RISK TABLE (1998–2024)   \n")
cat("====================================================================\n\n")

print(
  annual_absolute_risk %>% 
    select(year, csl_during_chinook, eulachon_during_chinook, shad_during_chinook, csl_threat_tier, eulachon_buffer_tier, overall_chinook_risk),
  n = 30
)

# Save Risk Table
write.csv(annual_absolute_risk, file.path(output_dir, "chinook_numeric_csl_exposure_risk_1998_2024.csv"), row.names = FALSE)
