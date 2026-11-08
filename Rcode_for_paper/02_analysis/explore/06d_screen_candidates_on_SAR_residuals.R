


#RESULTS--------
#using sablefish as peff (df_baseline$P_ocean_eff <- peff_df$Peff_sablefishAK_biomass3plus_smoltyr)
# TOP 15 CANDIDATES PREDICTING UNEXPLAINED SURVIVAL RESIDUALS
#                    indicator r_with_SAR r_with_resid abs_r_resid region      subguild missing_years
# 1      EGOA_CRG.HERR.MATBIOM      0.238        0.491       0.491   EGOA   Mid Trophic             0
# 2     EGOA_STKA.HERR.MATBIOM      0.244        0.480       0.480   EGOA   Mid Trophic             0
# 3       WGOA_PAV.SHRIMP.CPUE      0.158        0.433       0.433   WGOA   Mid Trophic             2
# 4        EGOA_SECM.CHIN.CPUE     -0.163       -0.350       0.350   EGOA   Mid Trophic             0
# 5        EGOA_SECM.EUPH.DENS      0.369        0.338       0.338   EGOA Lower Trophic             0
# 6  WGOA_SWDLN.SM.CAL.COP.FLL     -0.159       -0.185       0.185   WGOA Lower Trophic             0
# 7        EGOA_SECM.GAST.DENS     -0.301       -0.175       0.175   EGOA Lower Trophic             0
# 8     EGOA_SB.REPSUC.LSP.LAZ     -0.214       -0.150       0.150   EGOA       Seabird             2
# 9             EGOA_CHLA.BIOM     -0.120       -0.141       0.141   EGOA Lower Trophic             0
# 10    WGOA_PAV.EULACHON.CPUE      0.154       -0.140       0.140   WGOA   Mid Trophic             2
# 11   EGOA_SB.REPSUC.TBMU.LAZ     -0.171       -0.137       0.137   EGOA       Seabird             2
# 12     WGOA_PAV.CAPELIN.CPUE     -0.091        0.131       0.131   WGOA   Mid Trophic             2
# 13       EGOA_SECM.CHUM.CPUE     -0.004        0.129       0.129   EGOA   Mid Trophic             0
# 14       EGOA_SECM.COHO.CPUE      0.041       -0.128       0.128   EGOA   Mid Trophic             0
# 15            WGOA_CHLA.BIOM      0.083        0.105       0.105   WGOA Lower Trophic             0


#using sslseak as peff
# indicator r_with_SAR r_with_resid abs_r_resid region      subguild missing_years
# 1  WGOA_SWDLN.LG.CAL.COP.SPR     -0.439        0.412       0.412   WGOA Lower Trophic             0
# 2        EGOA_SECM.CHIN.CPUE     -0.163       -0.334       0.334   EGOA   Mid Trophic             0
# 3  WGOA_SWDLN.SM.CAL.COP.SPR     -0.107        0.292       0.292   WGOA Lower Trophic             0
# 4    EGOA_SB.REPSUC.TBMU.LAZ     -0.171       -0.284       0.284   EGOA       Seabird             2
# 5      WGOA_PAV.CAPELIN.CPUE     -0.091        0.229       0.229   WGOA   Mid Trophic             2
# 6  WGOA_SWDLN.LG.CAL.COP.FLL      0.273        0.215       0.215   WGOA Lower Trophic             0
# 7             WGOA_CHLA.BIOM      0.083        0.212       0.212   WGOA Lower Trophic             0
# 8        EGOA_SECM.PINK.CPUE      0.019        0.208       0.208   EGOA   Mid Trophic             0
# 9        EGOA_SECM.SOCK.CPUE     -0.021        0.195       0.195   EGOA   Mid Trophic             0
# 10            EGOA_CHLA.BIOM     -0.120       -0.184       0.184   EGOA Lower Trophic             0
# 11       EGOA_SECM.COHO.CPUE      0.041        0.184       0.184   EGOA   Mid Trophic             0
# 12   EGOA_SECM.COP.COMM.SIZE      0.102        0.183       0.183   EGOA Lower Trophic             0
# 13    EGOA_SB.REPSUC.LSP.LAZ     -0.214       -0.177       0.177   EGOA       Seabird             2
# 14       WGOA_PAV.JELLY.CPUE      0.224        0.149       0.149   WGOA   Mid Trophic             2
# 15   EGOA_SB.REPSUC.FTSP.LAZ     -0.196       -0.142       0.142   EGOA       Seabird             2
 


# ==============================================================================
# Script: 02_analysis/06d_screen_candidates_on_SAR_residuals.R
# Purpose: Screen candidate EGOA/WGOA constituents against the unexplained 
#          residual survival variance from the baseline model (DAG 1B without AKprey).
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------------------------
# 0. PATH SETUP & DATA LOADING
# ------------------------------------------------------------------------------
proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")

ak_const_path <- file.path(data_dir, "ak_dfa_raw_constituents_1998_2021.csv")
meta_path     <- file.path(data_dir, "ak_dfa_constituents_metadata_summary.csv")
alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")

ak_const_df <- read_csv(ak_const_path, show_col_types = FALSE)
meta_df     <- read_csv(meta_path, show_col_types = FALSE)
df_base_all <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)
peff_df     <- read_csv(peff_ts_path, show_col_types = FALSE)

# Safe year extraction
yr_col <- grep("^year$", names(df_base_all), value = TRUE, ignore.case = TRUE)
df_base_all$year <- if (length(yr_col) > 0) as.numeric(df_base_all[[yr_col[1]]]) else as.numeric(rownames(df_base_all))

# ------------------------------------------------------------------------------
# 1. ASSEMBLE BASELINE REGRESSION & EXTRACT RESIDUALS
# ------------------------------------------------------------------------------
df_baseline <- df_base_all %>%
  select(
    year,
    Growth  = X06_DFA_IGF_mu_smoltyr,
    CPUE    = X07_cpue_IntSprJunHW_smoltyr,
    SAR     = X16_DFA_sar_UCsprCh_smoltyr
  ) %>%
  filter(year >= 1998 & year <= 2021)

#df_baseline$P_ocean_eff <- peff_df$Peff_sablefishAK_biomass3plus_smoltyr
df_baseline$P_ocean_eff <- peff_df$Peff_ssl_seak_pup_pred_smoltyr

# Fit baseline regression (pure z-scores)
df_scaled <- df_baseline %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

base_lm <- lm(SAR ~ Growth + CPUE + P_ocean_eff, data = df_scaled)
base_r2 <- summary(base_lm)$r.squared

cat("\n====================================================================\n")
cat(sprintf("BASELINE MODEL (No AKprey): R2 = %.3f | Residual Variance = %.3f\n", 
            base_r2, 1 - base_r2))
cat("====================================================================\n")
print(summary(base_lm)$coefficients)

# Extract isolated unexplained survival variance
df_scaled$SAR_resid <- residuals(base_lm)

# ------------------------------------------------------------------------------
# 2. SCREEN CONSTITUENTS AGAINST SAR RESIDUALS
# ------------------------------------------------------------------------------
# Filter indicators with <= 2 missing years
robust_vars <- meta_df %>%
  filter(missing_years <= 2) %>%
  pull(indicator_col)

df_analysis <- df_scaled %>%
  select(year, SAR, SAR_resid) %>%
  left_join(ak_const_df %>% select(year, all_of(robust_vars)), by = "year")

# Calculate correlations with raw SAR vs. SAR_resid
screening_results <- lapply(robust_vars, function(v) {
  r_raw   <- cor(df_analysis[[v]], df_analysis$SAR, use = "pairwise.complete.obs")
  r_resid <- cor(df_analysis[[v]], df_analysis$SAR_resid, use = "pairwise.complete.obs")
  
  tibble(
    indicator    = v,
    r_with_SAR   = round(r_raw, 3),
    r_with_resid = round(r_resid, 3),
    abs_r_resid  = abs(round(r_resid, 3))
  )
}) %>%
  bind_rows() %>%
  left_join(meta_df %>% select(indicator = indicator_col, region, subguild, missing_years), by = "indicator") %>%
  arrange(desc(abs_r_resid))

cat("\n====================================================================\n")
cat("TOP 15 CANDIDATES PREDICTING UNEXPLAINED SURVIVAL RESIDUALS\n")
cat("====================================================================\n")
print(as.data.frame(head(screening_results, 15)))

#out_screen_path <- file.path(data_dir, "obj4_candidates_ranked_by_SAR_residual_sablefish.csv")
out_screen_path <- file.path(data_dir, "obj4_candidates_ranked_by_SAR_residual_sslseak.csv")
write_csv(screening_results, out_screen_path)

# ------------------------------------------------------------------------------
# 3. INTER-CORRELATION MATRIX OF TOP RESIDUAL EXPLAINERS
# ------------------------------------------------------------------------------
top_vars <- head(screening_results$indicator, 8)
cat("\nInter-correlation matrix among top 8 residual explainers:\n")
cat("(Look for pairs/trios with high shared correlation to form a stable latent)\n\n")
cor_mat <- cor(df_analysis[, top_vars], use = "pairwise.complete.obs")
print(round(cor_mat, 2))