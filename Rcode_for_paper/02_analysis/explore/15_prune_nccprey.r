

#                   NCCsardine SeaNettle Anchovy_gam JSOES_plankton marketsquid_GAM
# NCCsardine           1.000     0.686      -0.241         -0.931          -0.904
# SeaNettle            0.686     1.000       0.275         -0.636          -0.652
# Anchovy_gam         -0.241     0.275       1.000          0.368           0.398
# JSOES_plankton      -0.931    -0.636       0.368          1.000           0.972
# marketsquid_GAM     -0.904    -0.652       0.398          0.972           1.000

#Figure out which prey to include from importance matrix,  
#then see how correlated they are with each other (id groups of corr variables)
#then combine sets of the most correlated variables into latent nodes 
#prune until everything is significant
#compare pruned models by aic and model validity
#results: pink no growth no latent is best by far
suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(semPlot)
})

proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")
fig_dir  <- file.path(proj_dir, "Outputs/figures")
data_out_dir <- file.path(proj_dir, "Outputs/data")


val_ranks_wide <- read_csv("Rcode_for_paper/Outputs/data/obj2_data_variable_importance_wide.csv") 
print(val_ranks_wide %>% filter(grepl("X01|X02|X03|X04|X05",LisaName)) %>% arrange(desc(DAG1B_long)),n=Inf)


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
  Growth           = df_24[[find_col("IGF", df_24)]],
  P_ocean_eff    = peff_df$Peff_sablefishAK_biomass3plus_smoltyr,
  P_ocean_eff_sable    = peff_df$Peff_sablefishAK_biomass3plus_smoltyr,
  P_ocean_eff_ssl    = peff_df$Peff_ssl_seak_pup_pred_smoltyr,
  
  # Candidates from val_ranks_wide:
  NCCsardine       = df_24[[find_col("X05_DFA_abundSardine", df_24)]],
  SeaNettle       = df_24[[find_col("X02_DFA_SeaNettle", df_24)]],
  Anchovy_gam         = df_24[[find_col("X05_anchovy_GAM", df_24)]],
  JSOES_plankton  = df_24[[find_col("X01_DFA_sumPreyOfPrey_planktonJun", df_24)]],
  marketsquid_GAM     = df_24[[find_col("X04_marketsquid_GAM", df_24)]]

)

df_preds_scaled <- df_candidates %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

cand_cols <- c("NCCsardine","SeaNettle", "Anchovy_gam", "JSOES_plankton", "marketsquid_GAM")

cor_mat <- cor(df_fit[, c(cand_cols)], use = "pairwise.complete.obs")

cat("\n========================================================================================\n")
cat("STRATEGY 2 DIAGNOSTIC: PAIRWISE CORRELATION MATRIX ACROSS ALL 7 CANDIDATES\n")
cat("========================================================================================\n")
print(round(cor_mat, 3))
write.csv(as.data.frame(cor_mat),file.path(data_out_dir,"NCCPrey_cor_matrix.csv"),row.names = TRUE)



nccprey_cor<-read.csv(file.path(data_out_dir,"NCCPrey_cor_matrix.csv"),row.names=1);nccprey_cor
print(round(nccprey_cor, 3))


# ------------------------------------------------------------------------------
# 1b. SPECIFY DAG 1B WITH ONE BIPOLAR LATENT VARIABLES -- mammals------
# ------------------------------------------------------------------------------
# Model 1: Group 1 Alone (Mammal Latent)
sardine_latent <- "
  Latent_sardine =~ NCCsardine + 
                   JSOES_plankton + 
                   marketsquid_GAM + 
                   SeaNettle
  Growth ~   Latent_sardine

  SAR ~ CPUE +  Growth + P_ocean_eff
"

sardine_latent <- "
  Latent_sardine =~ NCCsardine + 
                   JSOES_plankton + 
                   marketsquid_GAM + 
                   SeaNettle
  Growth ~   Latent_sardine

"


fit_sardine <- sem(
  sardine_latent,
  data         = df_preds_scaled,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_sardine)
sem_output_summary_fxn(fit_sardine)

#prune growth--------

#all latent var p=0, but Latent_Mammal  p=0.45 -- NS!!!
#prune 1 at a time?
Latent_Mammal =~ DFA_Steller_s_l_WS_adultyr + 
  harbour_seal_DFO_WCVI_smoltyr + 
  Killer_whales_NR_BC_adultyr + 
  ssl_seak_pup_pred_smoltyr
mammal_latent <- "
  SAR ~ CPUE + Pink 
"

fit_mammal <- sem(
  mammal_latent,
  data         = df_preds_scaled,
  missing      = "ML",
  bounds       = "pos.var",
  optim.method = "BFGS",
  control      = list(iter.max = 2000)
)

summary(fit_mammal)
sem_output_summary_fxn(fit_mammal)
