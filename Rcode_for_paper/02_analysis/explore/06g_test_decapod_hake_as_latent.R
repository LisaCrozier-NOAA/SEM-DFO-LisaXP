
#WCVI should have been in NCC?
# SAR   CPUE   Hake Decapods   Pink  P_eff
# SAR       1.000  0.618 -0.630   -0.612 -0.254 -0.579
# CPUE      0.618  1.000 -0.865   -0.628  0.093 -0.359
# Hake     -0.630 -0.865  1.000    0.570  0.066  0.447
# Decapods -0.612 -0.628  0.570    1.000 -0.139  0.131
# Pink     -0.254  0.093  0.066   -0.139  1.000  0.034
# P_eff    -0.579 -0.359  0.447    0.131  0.034  1.000


# Model    AIC p_Chi2   CFI  SRMR SAR_R2       Latent_Load                                                                                                            SAR_Paths
# 1 Mod1_Latent_WarmPush_with_CPUE 183.41 0.0000 0.486 0.251  0.785 load=0.718 (p=NA) CPUE: b=0.081 (p=0.683) | Latent_WarmPush: b=-0.675 (p=0.008) | Pink: b=-0.318 (p=0.007) | P_eff: b=-0.437 (p=0.001)
# 2   Mod2_Latent_WarmPush_No_CPUE 181.58 0.1193 0.905 0.131  0.818 load=0.720 (p=NA)                           Latent_WarmPush: b=-0.726 (p=0.000) | Pink: b=-0.304 (p=0.008) | P_eff: b=-0.435 (p=0.000)
# 3    Mod3_Observed_Decapod_Triad  45.81     NA 1.000 0.000  0.728     Observed Only                                  Decapods: b=-0.592 (p=0.000) | Pink: b=-0.320 (p=0.003) | P_eff: b=-0.491 (p=0.000)
# 4    Mod4_Both_Separate_Observed  47.35     NA 1.000 0.000  0.734     Observed Only       Hake: b=-0.100 (p=0.494) | Decapods: b=-0.538 (p=0.000) | Pink: b=-0.307 (p=0.005) | P_eff: b=-0.454 (p=0.000)
> 
# ==============================================================================
# Script: 02_analysis/06g_test_decapod_hake_as_latent.R
# Purpose: Test correlation between decapod shelf biomass, Hake biomass,
#          and the other marine survival predictors across the 24-year series.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(corrplot)
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

df_warm_push <- tibble(
  year      = df_24$year,
  SAR       = df_24[[find_col("sar_UCsprCh", df_24)]],
  CPUE      = df_24[[find_col("cpue_IntSprJunHW", df_24)]],
  Hake      = df_24[[find_col("HakeAge5Plus", df_24)]],
  Decapods  = df_24[[find_col("biomassDecaShelfSum", df_24)]],
  Pink      = df_24[[find_col("pinkSalmonNorthAmerica", df_24)]],
  P_eff     = peff_df$Peff_sablefishAK_biomass3plus_smoltyr
)

cor_mat <- cor(df_warm_push %>% select(-year), use = "pairwise.complete.obs")

cat("\n====================================================================\n")
cat("CORRELATION MATRIX: WARM WATER ADVECTION & SURVIVAL\n")
cat("====================================================================\n")
print(round(cor_mat, 3))


# ==============================================================================
# Script: 02_analysis/06h_evaluate_warm_push_latent_vs_observed.R
# Purpose: Test if Hake and Decapods form a coherent latent construct
#          (Latent_WarmPush) or if Decapods alone dominates the structural model.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# Assemble dataset from your previous df_warm_push or df_fit
# Scale all variables (z-scores)
df_eval <- df_warm_push %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

# ------------------------------------------------------------------------------
# 1. SPECIFY MODEL ARCHITECTURES
# ------------------------------------------------------------------------------

# Model 1: Latent Warm-Water Push (Hake + Decapods with equated loadings)
# Testing if warm-water advection is an underlying unobserved factor
syntax_mod1 <- "
  Latent_WarmPush =~ a*Hake + a*Decapods
  SAR ~ CPUE + Latent_WarmPush + Pink + P_eff
"

# Model 2: Latent Warm-Water Push WITHOUT CPUE (since CPUE collides with Hake)
syntax_mod2 <- "
  Latent_WarmPush =~ a*Hake + a*Decapods
  SAR ~ Latent_WarmPush + Pink + P_eff
"

# Model 3: Observed Triad (Decapods alone without Hake or CPUE)
syntax_mod3 <- "
  SAR ~ Decapods + Pink + P_eff
"

# Model 4: Both Hake and Decapods as separate observed predictors
syntax_mod4 <- "
  SAR ~ Hake + Decapods + Pink + P_eff
"

models <- list(
  "Mod1_Latent_WarmPush_with_CPUE" = syntax_mod1,
  "Mod2_Latent_WarmPush_No_CPUE"   = syntax_mod2,
  "Mod3_Observed_Decapod_Triad"    = syntax_mod3,
  "Mod4_Both_Separate_Observed"    = syntax_mod4
)

# ------------------------------------------------------------------------------
# 2. FIT AND EXTRACT DIAGNOSTICS
# ------------------------------------------------------------------------------
eval_summary <- list()

for (m_name in names(models)) {
  fit <- tryCatch(
    sem(
      models[[m_name]], 
      data    = df_eval, 
      missing = "ML", 
      bounds  = "pos.var",
      control = list(iter.max = 1000)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !lavInspect(fit, "converged")) {
    cat(sprintf("[-] %s failed to converge.\n", m_name))
    next
  }
  
  pe <- parameterEstimates(fit, standardized = TRUE)
  
  # Format structural regressions
  paths <- pe %>%
    filter(lhs == "SAR", op == "~") %>%
    mutate(desc = sprintf("%s: b=%.3f (p=%.3f)", rhs, std.all, pvalue)) %>%
    pull(desc) %>%
    paste(collapse = " | ")
  
  # Check latent loading if applicable
  loading_str <- if (grepl("Latent", m_name)) {
    load_rows <- pe %>% filter(op == "=~")
    sprintf("load=%.3f (p=%.3f)", load_rows$std.all[1], load_rows$pvalue[1])
  } else {
    "Observed Only"
  }
  
  eval_summary[[m_name]] <- tibble(
    Model        = m_name,
    AIC          = round(as.numeric(fitMeasures(fit, "aic")), 2),
    p_Chi2       = round(as.numeric(fitMeasures(fit, "pvalue")), 4),
    CFI          = round(as.numeric(fitMeasures(fit, "cfi")), 3),
    SRMR         = round(as.numeric(fitMeasures(fit, "srmr")), 3),
    SAR_R2       = round(as.numeric(inspect(fit, "r2")["SAR"]), 3),
    Latent_Load  = loading_str,
    SAR_Paths    = paths
  )
}

results_df <- bind_rows(eval_summary)
cat("\n====================================================================\n")
cat("EVALUATION: WARM-WATER PUSH ARCHITECTURES\n")
cat("====================================================================\n")
print(as.data.frame(results_df))
