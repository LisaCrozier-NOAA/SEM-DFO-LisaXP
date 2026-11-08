

#RESULTS----
# print(as.data.frame(head(sar_ranking, 15)))
#                     indicator r_with_SAR abs_r region      subguild missing_years
# 1          WGOA_MID.IL.CAPELIN      0.554 0.554   WGOA   Mid Trophic             0
# 2    WGOA_SWDLN.LG.CAL.COP.SPR     -0.439 0.439   WGOA Lower Trophic             0
# 3          EGOA_SECM.EUPH.DENS      0.369 0.369   EGOA Lower Trophic             0
# 4  WGOA_SWDLN.COP.CMMTY.SZ.SPR     -0.325 0.325   WGOA Lower Trophic             0
# 5          EGOA_SECM.GAST.DENS     -0.301 0.301   EGOA Lower Trophic             0
# 6      EGOA_SECM.MESOZOOP.DENS      0.274 0.274   EGOA   Mid Trophic             0
# 7    WGOA_SWDLN.LG.CAL.COP.FLL      0.273 0.273   WGOA Lower Trophic             0
# 8       EGOA_STKA.HERR.MATBIOM      0.244 0.244   EGOA   Mid Trophic             0
# 9        EGOA_CRG.HERR.MATBIOM      0.238 0.238   EGOA   Mid Trophic             0
# 10         WGOA_PAV.JELLY.CPUE      0.224 0.224   WGOA   Mid Trophic             2
# 11      EGOA_SB.REPSUC.LSP.LAZ     -0.214 0.214   EGOA       Seabird             2
# 12          EGOA_SECM.HYP.AMPH      0.210 0.210   EGOA Lower Trophic             0
# 13     EGOA_SB.REPSUC.FTSP.LAZ     -0.196 0.196   EGOA       Seabird             2
# 14   EGOA_SECM.LG.CAL.COP.DENS      0.191 0.191   EGOA Lower Trophic             0
# 15     EGOA_SB.REPSUC.TBMU.LAZ     -0.171 0.171   EGOA       Seabird             2

# ==============================================================================
# Script: 02_analysis/06c_screen_constituent_correlations.R
# Purpose: Build a correlation matrix across all EGOA/WGOA constituents (<= 2 NAs)
#          and SAR to identify viable latent clusters before fitting SEM.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(corrplot)
})

proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")

ak_const_path <- file.path(data_dir, "ak_dfa_raw_constituents_1998_2021.csv")
meta_path     <- file.path(data_dir, "ak_dfa_constituents_metadata_summary.csv")
alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")

ak_const_df <- read_csv(ak_const_path, show_col_types = FALSE)
meta_df     <- read_csv(meta_path, show_col_types = FALSE)
df_base_all <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)

# Extract SAR
yr_col <- grep("^year$", names(df_base_all), value = TRUE, ignore.case = TRUE)
df_base_all$year <- if (length(yr_col) > 0) as.numeric(df_base_all[[yr_col[1]]]) else as.numeric(rownames(df_base_all))

df_sar <- df_base_all %>%
  select(year, SAR = X16_DFA_sar_UCsprCh_smoltyr) %>%
  filter(year >= 1998 & year <= 2021)

# Filter candidate indicators with <= 2 missing years
robust_vars <- meta_df %>%
  filter(missing_years <= 2) %>%
  pull(indicator_col)

# Merge
df_all <- df_sar %>%
  left_join(ak_const_df %>% select(year, all_of(robust_vars)), by = "year") %>%
  select(-year)

# 1. Correlations with SAR
sar_cors <- sapply(robust_vars, function(v) {
  cor(df_all[[v]], df_all$SAR, use = "pairwise.complete.obs")
})

sar_ranking <- tibble(
  indicator = names(sar_cors),
  r_with_SAR = round(as.numeric(sar_cors), 3),
  abs_r = abs(r_with_SAR)
) %>%
  left_join(meta_df %>% select(indicator = indicator_col, region, subguild, missing_years), by = "indicator") %>%
  arrange(desc(abs_r))

cat("\n====================================================================\n")
cat("TOP 15 CONSTITUENTS CORRELATED WITH SAR (EGOA & WGOA)\n")
cat("====================================================================\n")
print(as.data.frame(head(sar_ranking, 15)))

# 2. Inter-correlation matrix of top candidates
top_candidates <- head(sar_ranking$indicator, 12)
cor_mat <- cor(df_all[, top_candidates], use = "pairwise.complete.obs")

cat("\nInter-correlation matrix among top 12 candidates:\n")
print(round(cor_mat, 2))

# Export tables
write_csv(sar_ranking, file.path(data_dir, "obj4_goa_constituent_correlations_with_SAR.csv"))
