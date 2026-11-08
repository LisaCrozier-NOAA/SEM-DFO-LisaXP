

# ==============================================================================
# Script: 02_analysis/06a_extract_goa_dfa_constituents.R
# Purpose: Extract raw time-series constituents of EGOA and WGOA DFAs 
#          (LTL, MTL, and Seabirds) for latent subguild testing in lavaan (FIML).
# Output:
#   - Outputs/data/ak_dfa_raw_constituents_1998_2021.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------------------------
# STEP 0: DIRECTORY SETUP & PATHS
# ------------------------------------------------------------------------------
goa_dir      <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Ferris-DFAIndicators-goa/data"
proj_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper"
data_out_dir <- file.path(proj_dir, "Outputs/data")
dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

egoa_data_path <- file.path(goa_dir, "EGOA_EcoState_Data_Jan2023.csv")
wgoa_data_path <- file.path(goa_dir, "WGOA_EcoState_Data_Jan2023.csv")

if (!file.exists(egoa_data_path)) stop("Missing EGOA raw data at: ", egoa_data_path)
if (!file.exists(wgoa_data_path)) stop("Missing WGOA raw data at: ", wgoa_data_path)

# ------------------------------------------------------------------------------
# STEP 1: EXTRACT EGOA CONSTITUENTS (LTL, MTL, SEABIRDS)
# ------------------------------------------------------------------------------
cat("Extracting EGOA constituents (LTL, MTL, Seabirds)...\n")
e_data_raw <- read.csv(egoa_data_path, na.strings = c("NA", "na", ""))
names(e_data_raw)[1] <- "year"

# 1. Lower Trophic Level
egoa_ltl_codes <- c(
  "CHLA.BIOM", "CHLA.PEAK", "SECM.COP.COMM.SIZE", "SECM.EUPH.DENS", 
  "SECM.GAST.DENS", "SECM.HYP.AMPH", "SECM.LG.CAL.COP.DENS", 
  "SECM.MESOZOOP.DENS", "SECM.SM.CAL.COP.DENS"
)

# 2. Mid Trophic Level
egoa_mtl_codes <- c(
  "CRG.HERR.MATBIOM", "GB.HUMP.CRUD.BR.JA", "SECM.CHIN.CPUE", 
  "SECM.CHUM.CPUE", "SECM.COHO.CPUE", "SECM.PINK.CPUE", 
  "SECM.SOCK.CPUE", "STKA.HERR.MATBIOM"
)

# 3. Seabirds (Productivity and Mean Hatch Dates)
egoa_sb_codes <- c(
  "SB.REPSUC.COMU.LAZ", "SB.REPSUC.FTSP.LAZ", "SB.REPSUC.LSP.LAZ", 
  "SB.REPSUC.TBMU.LAZ", "SBRD.MEANHATCH.COMU.LAZ", "SBRD.MEANHATCH.FTSP.LAZ", 
  "SBRD.MEANHATCH.GWGU.LAZ", "SBRD.MEANHATCH.LHSP.LAZ", "SBRD.MEANHATCH.RHAU.LAZ", 
  "SBRD.MEANHATCH.TBMU.LAZ"
)

egoa_target_codes <- unique(c(egoa_ltl_codes, egoa_mtl_codes, egoa_sb_codes))

df_egoa_wide <- e_data_raw %>%
  filter(year >= 1998, year <= 2021) %>%
  select(year, any_of(egoa_target_codes)) %>%
  rename_with(~ paste0("EGOA_", .x), -year)

# ------------------------------------------------------------------------------
# STEP 2: EXTRACT WGOA CONSTITUENTS (LTL, MTL, SEABIRDS)
# ------------------------------------------------------------------------------
cat("Extracting WGOA constituents (LTL, MTL, Seabirds)...\n")
w_data_raw <- read.csv(wgoa_data_path, na.strings = c("NA", "na", ""))
names(w_data_raw)[1] <- "year"

# 1. Lower Trophic Level
wgoa_ltl_codes <- c(
  'CHLA.BIOM', 'CHLA.PEAK', 'FOCI.ICH.AMMPER', 'FOCI.ICH.ATHSTO', 
  'FOCI.ICH.BATSPP', 'FOCI.ICH.GADCHA', 'FOCI.ICH.GADMAC', 
  'FOCI.ICH.HIPELA', 'FOCI.ICH.HIPSTE', 'FOCI.ICH.LEPBIL', 
  'FOCI.ICH.LEPPOL', 'FOCI.ICH.PLASTE', 'FOCI.ICH.SEBSPP', 
  'FOCI.ICH.STELEU', 'SWDLN.COP.CMMTY.SZ.SPR', 'SWDLN.LG.CAL.COP.FLL', 
  'SWDLN.LG.CAL.COP.SPR', 'SWDLN.SM.CAL.COP.FLL', 'SWDLN.SM.CAL.COP.SPR'
)

# 2. Mid Trophic Level
wgoa_mtl_codes <- c(
  'CHI.EULACHON.CPUE', 'CHI.JELLY.CPUE', 'CHI.SHRIMP.CPUE', 
  'MID.IL.CAPELIN', 'PAV.SHRIMP.CPUE', 'PAV.CAPELIN.CPUE', 
  'PAV.EULACHON.CPUE', 'PAV.JELLY.CPUE'
)

# 3. Seabirds (Productivity and Mean Hatch Dates)
wgoa_sb_codes <- c(
  'SB.MEANHATCH.COMU.CHOW', 'SB.MEANHATCH.TBMU.CHOW', 'SB.MEANHATCH.PAAU.CHOW', 
  'SB.MEANHATCH.HOPU.CHOW', 'SB.MEANHATCH.TUPU.CHOW', 'SB.MEANHATCH.GWGU.CHOW', 
  'SB.MEANHATCH.BLKW.EAMT', 'SB.REPSUC.COMU.CHOW', 'SB.REPSUC.TBMU.CHOW', 
  'SB.REPSUC.HPUFF.CHOW', 'SB.REPSUC.TPUF.CHOW', 'SB.REPSUC.BLK.CHOW', 
  'SB.REPSUC.PAUK.CHOW', 'SB.REPSUC.BLK.AMAT', 'SB.REPSUC.FTSP.AMAT'
)

wgoa_target_codes <- unique(c(wgoa_ltl_codes, wgoa_mtl_codes, wgoa_sb_codes))

df_wgoa_wide <- w_data_raw %>%
  filter(year >= 1998, year <= 2021) %>%
  select(year, any_of(wgoa_target_codes)) %>%
  rename_with(~ paste0("WGOA_", .x), -year)

# ------------------------------------------------------------------------------
# STEP 3: MERGE & REPORT MISSINGNESS PROFILES FOR FIML
# ------------------------------------------------------------------------------
df_goa_constituents <- full_join(df_egoa_wide, df_wgoa_wide, by = "year") %>%
  arrange(year)

# Calculate missingness counts across the 24-year complete frame (1998–2021)
na_counts <- colSums(is.na(df_goa_constituents[, names(df_goa_constituents) != "year"]))

cat("\n====================================================================\n")
cat("CONSTITUENT EXTRACTION SUMMARY (1998-2021, N = 24 Years)\n")
cat("====================================================================\n")
cat(sprintf("Total Indicators Extracted: %d\n", length(na_counts)))
cat(sprintf("  - 100%% Complete (0 missing years): %d\n", sum(na_counts == 0)))
cat(sprintf("  - 1 to 2 missing years:            %d\n", sum(na_counts >= 1 & na_counts <= 2)))
cat(sprintf("  - 3 to 5 missing years:            %d\n", sum(na_counts >= 3 & na_counts <= 5)))
cat(sprintf("  - > 5 missing years:               %d\n", sum(na_counts > 5)))

# Export all extracted series along with a missingness diagnostic summary
out_file <- file.path(data_out_dir, "ak_dfa_raw_constituents_1998_2021.csv")
write_csv(df_goa_constituents, out_file)

meta_summary <- tibble(
  indicator_col = names(na_counts),
  missing_years = as.integer(na_counts),
  region        = ifelse(grepl("^EGOA_", names(na_counts)), "EGOA", "WGOA"),
  subguild      = case_when(
    grepl("SB|SBRD", names(na_counts))  ~ "Seabird",
    grepl("CHLA|COP|EUPH|GAST|AMPH|FOCI|SWDLN", names(na_counts)) ~ "Lower Trophic",
    TRUE ~ "Mid Trophic"
  )
) %>% arrange(missing_years, subguild)

write_csv(meta_summary, file.path(data_out_dir, "ak_dfa_constituents_metadata_summary.csv"))

cat(sprintf("\nMaster constituents table saved to:\n  -> %s\n", out_file))
cat(sprintf("Metadata summary saved to:\n  -> %s\n", file.path(data_out_dir, "ak_dfa_constituents_metadata_summary.csv")))

