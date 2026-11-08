
#RESULTS:
#Correltaion matrix among akprey indices
#                     SAR Decapods Mysids EGOA_lower WGOA_mid WGOA_lower EGOA_seabirds Pollock_age1_2
# SAR             1.000   -0.612  0.219      0.424   -0.145      0.028        -0.048          0.172
# Decapods       -0.612    1.000 -0.200     -0.327    0.535     -0.570         0.198         -0.218
# Mysids          0.219   -0.200  1.000      0.572    0.286      0.010        -0.371         -0.327
# EGOA_lower      0.424   -0.327  0.572      1.000   -0.074      0.120        -0.191         -0.285
# WGOA_mid       -0.145    0.535  0.286     -0.074    1.000     -0.777         0.331         -0.163
# WGOA_lower      0.028   -0.570  0.010      0.120   -0.777      1.000        -0.602         -0.015
# EGOA_seabirds  -0.048    0.198 -0.371     -0.191    0.331     -0.602         1.000          0.111
# Pollock_age1_2  0.172   -0.218 -0.327     -0.285   -0.163     -0.015         0.111          1.000

#                 Decapods Mysids EGOA_lower WGOA_mid WGOA_lower EGOA_seabirds Pollock_age1_2
# Decapods          1.000 -0.200     -0.327    0.535     -0.570         0.198         -0.218
# Mysids           -0.200  1.000      0.572    0.286      0.010        -0.371         -0.327
# EGOA_lower       -0.327  0.572      1.000   -0.074      0.120        -0.191         -0.285
# WGOA_mid          0.535  0.286     -0.074    1.000     -0.777         0.331         -0.163
# WGOA_lower       -0.570  0.010      0.120   -0.777      1.000        -0.602         -0.015
# EGOA_seabirds     0.198 -0.371     -0.191    0.331     -0.602         1.000          0.111
# Pollock_age1_2   -0.218 -0.327     -0.285   -0.163     -0.015         0.111          1.000

#from the correlation matrix, it looks like there might be 2 groupings -- 
#1 with decapods, wgoa_mid and wgoa_lower, and 
#another with mysids, egoa_lower, egoa _seabirds, and pollock1and2. 

# ==============================================================================
# Script: 02_analysis/11a_corr_grouping_workflows.R
# Purpose: Benchmark two formal SEM starting strategies:
#          1. Performance-First (grouped by DAG 1B Akaike weights)
#          2. Coherence-First (grouped by empirical pairwise correlation)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")

alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")

df_raw  <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)
peff_df <- read_csv(peff_ts_path, show_col_types = FALSE)

yr_col <- grep("^year$", names(df_raw), value = TRUE, ignore.case = TRUE)
df_raw$year <- if (length(yr_col) > 0) as.numeric(df_raw[[yr_col[1]]]) else as.numeric(rownames(df_raw))
df_24 <- df_raw %>% filter(year >= 1998 & year <= 2021)

find_col <- function(pattern, df) {
  matches <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
  if (length(matches) == 0) stop("Could not find column matching: ", pattern)
  matches[1]
}

# ------------------------------------------------------------------------------
# 1. ASSEMBLE CANDIDATE POOL (7 CANDIDATES FROM RANKING TABLE)
# ------------------------------------------------------------------------------
df_candidates <- tibble(
  year           = df_24$year,
  SAR            = df_24[[find_col("sar_UCsprCh", df_24)]],
  CPUE           = df_24[[find_col("cpue_IntSprJunHW", df_24)]],
  P_ocean_eff    = peff_df$Peff_sablefishAK_biomass3plus_smoltyr,
  
  # The 7 candidates from val_ranks_wide:
  Pink       = df_24[[find_col("pinkSalmonNorthAmerica", df_24)]],
  Decapods       = df_24[[find_col("biomassDecaShelfSum", df_24)]],
  Mysids         = df_24[[find_col("biomassMysidShelfSum", df_24)]],
  EGOA_lower     = df_24[[find_col("EGOA.*lower", df_24)]],
  WGOA_mid       = df_24[[find_col("WGOA.*mid", df_24)]],
  WGOA_lower     = df_24[[find_col("WGOA.*lower", df_24)]],
  EGOA_seabirds  = df_24[[find_col("EGOA.*seabird|EGOA.*sb", df_24)]],
  Pollock_age1_2 = df_24[[find_col("pollockGOA_age1and2", df_24)]]
)

df_fit <- df_candidates %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

cand_cols <- c("Pink","Decapods", "Mysids", "EGOA_lower", "WGOA_mid", "WGOA_lower", "EGOA_seabirds", "Pollock_age1_2")

cor_mat <- cor(df_fit[, c(cand_cols)], use = "pairwise.complete.obs")

cat("\n========================================================================================\n")
cat("STRATEGY 2 DIAGNOSTIC: PAIRWISE CORRELATION MATRIX ACROSS ALL 7 CANDIDATES\n")
cat("========================================================================================\n")
print(round(cor_mat, 3))
write.csv(as.data.frame(cor_mat),file.path(data_out_dir,"AKPrey_cor_matrix.csv"),row.names = TRUE)

#BELOW NOT USED --DUMP=====================
# # ------------------------------------------------------------------------------
# # 3. HELPER: FIT AND EVALUATE LATENT SEM
# # ------------------------------------------------------------------------------
# fit_latent_block <- function(indicators, strategy_name) {
#   if (length(indicators) == 2) {
#     meas_block <- sprintf("Latent_AK =~ a*%s + a*%s", indicators[1], indicators[2])
#   } else {
#     meas_block <- paste0("Latent_AK =~ ", paste(indicators, collapse = " + "))
#   }
#   
#   struct_block <- "SAR ~ CPUE + Latent_AK + P_ocean_eff"
#   syntax       <- paste(meas_block, struct_block, sep = "\n")
#   
#   fit <- tryCatch(
#     sem(
#       syntax, 
#       data         = df_fit, 
#       missing      = "ML", 
#       bounds       = "pos.var",
#       optim.method = "BFGS",
#       control      = list(iter.max = 1000)
#     ),
#     error = function(e) NULL
#   )
#   
#   if (is.null(fit) || !lavInspect(fit, "converged")) {
#     return(tibble(
#       Strategy = strategy_name, n_ind = length(indicators), Indicators = paste(indicators, collapse = ", "),
#       AIC = NA_real_, CFI = NA_real_, SRMR = NA_real_, SAR_R2 = NA_real_,
#       Latent_Beta = NA_character_, Heywood = "DID NOT CONVERGE"
#     ))
#   }
#   
#   pe <- parameterEstimates(fit, standardized = TRUE)
#   loadings <- pe %>% filter(op == "=~")
#   ov_vars  <- pe %>% filter(op == "~~", lhs == rhs, lhs %in% indicators)
#   sar_path <- pe %>% filter(lhs == "SAR", rhs == "Latent_AK")
#   
#   has_heywood <- any(ov_vars$est <= 0.001) || any(abs(loadings$std.all) > 1.0)
#   
#   tibble(
#     Strategy    = strategy_name,
#     n_ind       = length(indicators),
#     Indicators  = paste(indicators, collapse = ", "),
#     AIC         = round(as.numeric(fitMeasures(fit, "aic")), 2),
#     CFI         = round(as.numeric(fitMeasures(fit, "cfi")), 3),
#     SRMR        = round(as.numeric(fitMeasures(fit, "srmr")), 3),
#     SAR_R2      = round(as.numeric(inspect(fit, "r2")["SAR"]), 3),
#     Latent_Beta = sprintf("%.3f (p=%.3f)", sar_path$std.all[1], sar_path$pvalue[1]),
#     Heywood     = ifelse(has_heywood, "YES", "NO")
#   )
# }
# 
# # ------------------------------------------------------------------------------
# # 4. BENCHMARKING BOTH STRATEGIES
# # ------------------------------------------------------------------------------
# runs <- list()
# 
# # Strategy 1: Top performers by DAG1B_long weights
# # Rank 1: Decapods (0.719), Rank 2: Mysids (0.095), Rank 3: EGOA_lower (0.038), Rank 4: WGOA_mid (0.037)
# runs[[1]] <- fit_latent_block(c("Decapods", "Mysids"), "Strat1_Top2_Performers")
# runs[[2]] <- fit_latent_block(c("Decapods", "Mysids", "EGOA_lower"), "Strat1_Top3_Performers")
# runs[[3]] <- fit_latent_block(c("Decapods", "Mysids", "EGOA_lower", "WGOA_mid"), "Strat1_Top4_Performers")
# 
# # Strategy 2: Coherence-First (Clusters identified from cor_mat)
# # Group A: Crustacean Shelf Cluster (Decapods + Mysids)
# runs[[4]] <- fit_latent_block(c("Decapods", "Mysids"), "Strat2_Crustacean_Cluster")
# 
# # Group B: Regional Lower Trophic Pelagic Pair (EGOA_lower + WGOA_lower)
# runs[[5]] <- fit_latent_block(c("EGOA_lower", "WGOA_lower"), "Strat2_Regional_LowerTrophic")
# 
# # Group C: Offshore Forage & Seabird Guild (WGOA_mid + EGOA_seabirds)
# runs[[6]] <- fit_latent_block(c("WGOA_mid", "EGOA_seabirds"), "Strat2_MidTrophic_Seabirds")
# 
# # Group D: Multi-Trophic Pelagic Triangle (WGOA_lower + WGOA_mid + EGOA_seabirds)
# runs[[7]] <- fit_latent_block(c("WGOA_lower", "WGOA_mid", "EGOA_seabirds"), "Strat2_GOA_DFA_Triangle")
# 
# comparison_df <- bind_rows(runs)
# 
# cat("\n========================================================================================\n")
# cat("BENCHMARK SUMMARY: PERFORMANCE-FIRST vs. COHERENCE-FIRST STARTING POINTS\n")
# cat("========================================================================================\n")
# print(as.data.frame(comparison_df))
# 
# write_csv(comparison_df, file.path(data_dir, "obj4_systematic_starting_points_benchmark.csv"))