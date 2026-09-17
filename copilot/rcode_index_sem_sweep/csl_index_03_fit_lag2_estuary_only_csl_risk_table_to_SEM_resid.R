
#RESULT
write.csv(guild_updated, "copilot/outputs_6/guild_dfas1_24yr_ishark_issl_icsl.csv", row.names = FALSE)



# idxE10_csl_chinook_ratio  = safe_z(log1p(csl_during_chinook) - log1p(chinook_estuary_total)),
# print(head(res_cmp, 12))
# # A tibble: 12 × 12
# csl_term                     n sar_r2   aic   bic   cfi  rmsea csl_est     csl_p delta_r2_vs_base delta_aic_vs_base delta_bic_vs_base
# <chr>                    <int>  <dbl> <dbl> <dbl> <dbl>  <dbl>   <dbl>     <dbl>            <dbl>             <dbl>             <dbl>
#   1 idxE10_csl_chinook_ratio    24  0.820  101.  111. 1     0       -0.615 0.0000608           0.119            -10.3              -9.13 
# 2 idxE8_mult_weekly           24  0.762  108.  119. 0.999 0.0235  -0.285 0.0174              0.0611            -3.08             -1.90 
# 3 idxE2_csl_only_weekly       24  0.761  108.  119. 1     0       -0.356 0.0230              0.0604            -2.68             -1.50 
# 4 csl_week_z                  24  0.761  108.  119. 1     0       -0.356 0.0230              0.0604            -2.68             -1.50 
# 5 idxE11_csl_eul_ratio        24  0.759  110.  121. 0.850 0.260   -0.353 0.0846              0.0586            -0.803             0.375
# 6 est_risk_eul_heavy          24  0.756  107.  118. 0.966 0.120    0.343 0.00998             0.0551            -3.86             -2.68 
# 7 idxE1_csl_only_count        24  0.754  109.  120. 1     0       -0.304 0.0384              0.0532            -1.95             -0.767
# 8 csl_z                       24  0.754  109.  120. 1     0       -0.304 0.0384              0.0532            -1.95             -0.767
# 9 idxE7_mult_count            24  0.749  109.  120. 0.995 0.0425  -0.242 0.0455              0.0486            -1.70             -0.521
# 10 idxE9_mult_peak             24  0.730  111.  121. 0.983 0.0803  -0.188 0.112               0.0297            -0.400             0.778
# 11 idxE3_csl_only_peak         24  0.728  111.  122. 1     0       -0.208 0.138               0.0279            -0.109             1.07 
# 12 csl_peak_z                  24  0.728  111.  122. 1     0       -0.208 0.138               0.0279            -0.109             1.07 

#Explanation:
# This is a log predator-to-prey pressure index (standardized):
#   
#   log1p(csl_during_chinook) increases with more sea lions.
# log1p(chinook_estuary_total) increases with more salmon.
# Subtracting them means the index is high when:
#   CSL are high and/or
# Chinook are low.
# So it approximates:
#   log(  1  +    CSL
#         --------------------
#         1  +  Chinook
#       )
# then z-scored.
# 
# That’s why it can outperform raw CSL count: it captures relative predator pressure per available salmon, not just predator abundance.
# 
# Why coefficient is negative in your SEM result
# You got csl_est = -0.615 (significant).
# Interpretation depends on residual structure/path model context, but mechanically it means:
#   
#   higher predator:prey pressure index predicts lower x16_sar (or lower than expected after other terms), which is directionally plausible.

# ---------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

out_dir <- "copilot/outputs_6"

# 0) Paths
# ---------------------------
path_guild<-  file.path(out_dir, "guild_dfas1_24yr.csv")
path_csl   <- file.path(out_dir, "chinook_numeric_csl_exposure_risk_1998_2024.csv")

if (!file.exists(path_guild)) stop("Missing: ", path_guild)
if (!file.exists(path_csl))   stop("Missing: ", path_csl)

# ---------------------------
# 1) Load
# ---------------------------
guild <- read.csv(path_guild) %>%
  clean_names() %>%
  select(-any_of("x10_harbor_seal_cr_2yr_lead"))

csl <- read.csv(path_csl) %>% clean_names()

need_csl <- c(
  "year","csl_during_chinook","csl_weekly_avg_chin","csl_annual_peak",
  "chinook_estuary_total","eulachon_during_chinook","shad_during_chinook"
)
miss <- setdiff(need_csl, names(csl))
if (length(miss) > 0) stop("CSL estuary file missing: ", paste(miss, collapse = ", "))

safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - m) / s)
}

# ---------------------------
# 2) Build estuary-only candidate indices
# (higher index = higher expected predation pressure)
# eulachon weighted more than shad
# ---------------------------
csl2 <- csl %>%
  mutate(
    csl_z       = safe_z(csl_during_chinook),
    csl_week_z  = safe_z(csl_weekly_avg_chin),
    csl_peak_z  = safe_z(csl_annual_peak),
    chinook_z   = safe_z(chinook_estuary_total),
    eul_z       = safe_z(eulachon_during_chinook),
    shad_z      = safe_z(shad_during_chinook),
    
    # prey buffers (higher prey => lower risk), eulachon emphasized
    est_buf_eul_heavy = 0.75 * eul_z + 0.25 * shad_z,
    est_risk_eul_heavy = -(0.75 * eul_z + 0.25 * shad_z),
    
    # candidate risk indices
    idxE1_csl_only_count      = csl_z,
    idxE2_csl_only_weekly     = csl_week_z,
    idxE3_csl_only_peak       = csl_peak_z,
    idxE4_csl_minus_buffer    = safe_z(csl_z + est_risk_eul_heavy),
    idxE5_csl_week_minus_buf  = safe_z(csl_week_z + est_risk_eul_heavy),
    idxE6_csl_peak_minus_buf  = safe_z(csl_peak_z + est_risk_eul_heavy),
    idxE7_mult_count          = safe_z(csl_z * (1 + est_risk_eul_heavy)),
    idxE8_mult_weekly         = safe_z(csl_week_z * (1 + est_risk_eul_heavy)),
    idxE9_mult_peak           = safe_z(csl_peak_z * (1 + est_risk_eul_heavy)),
    idxE10_csl_chinook_ratio  = safe_z(log1p(csl_during_chinook) - log1p(chinook_estuary_total)),
    idxE11_csl_eul_ratio      = safe_z(log1p(csl_during_chinook) - log1p(eulachon_during_chinook + 1)),
    idxE12_csl_shad_ratio     = safe_z(log1p(csl_during_chinook) - log1p(shad_during_chinook + 1))
  ) %>%
  select(year, starts_with("idxE"), csl_z, csl_week_z, csl_peak_z, est_risk_eul_heavy)

# ---------------------------
# 3) Apply 2-year lead alignment for adult return
# For SAR year t, use CSL metrics from t+2
# ---------------------------
# confirm candidate exists
if (!"est_risk_eul_heavy" %in% names(csl2)) {
  stop("est_risk_eul_heavy was not created in csl2.")
}

# build lag2 table including all idxE* + est_risk_eul_heavy
csl_lag2_for_sar <- csl2 %>%
  transmute(
    year = year - 2,
    est_risk_eul_heavy,
    idxE1_csl_only_count,
    idxE2_csl_only_weekly,
    idxE3_csl_only_peak,
    idxE4_csl_minus_buffer,
    idxE5_csl_week_minus_buf,
    idxE6_csl_peak_minus_buf,
    idxE7_mult_count,
    idxE8_mult_weekly,
    idxE9_mult_peak,
    idxE10_csl_chinook_ratio,
    idxE11_csl_eul_ratio,
    idxE12_csl_shad_ratio
  )

# then candidate terms:
cand_terms <- names(csl_lag2_for_sar) %>% setdiff("year")

# optional: prioritize display by AIC and also by R2
# after res_cmp is built:
res_cmp_aic <- res_cmp %>% arrange(aic)
res_cmp_r2  <- res_cmp %>% arrange(desc(delta_r2_vs_base), delta_aic_vs_base)

write.csv(res_cmp_aic, file.path(out_dir, "csl_estuary_lag2_model_comparison_sortedAIC.csv"), row.names = FALSE)
write.csv(res_cmp_r2,  file.path(out_dir, "csl_estuary_lag2_model_comparison_sortedR2.csv"), row.names = FALSE)

cat("\nTop by AIC:\n")
print(head(res_cmp_aic, 12))

cat("\nTop by ΔR2:\n")
print(head(res_cmp_r2, 12))
# ---------------------------
# 4) Build modeling table
# ---------------------------
req_guild <- c(
  "year","x07_dfa_cpue_int_spr_jun_hw","x09_dfa_hake_age5plus","x16_sar",
  "x15_dfa_sleeper_shark_bsai_pred_ak","x15_issl_z"
)
missg <- setdiff(req_guild, names(guild))
if (length(missg) > 0) stop("Guild missing: ", paste(missg, collapse = ", "))

d0 <- guild %>%
  left_join(csl_lag2_for_sar, by = "year") %>%
  filter(year >= 1998, year <= 2021)

# ---------------------------
# 5) SEM fit helper
# ---------------------------
fit_one <- function(df, csl_term = NULL) {
  base_keep <- c("year","x07_dfa_cpue_int_spr_jun_hw","x09_dfa_hake_age5plus","x16_sar",
                 "x15_dfa_sleeper_shark_bsai_pred_ak","x15_issl_z")
  keep <- if (is.null(csl_term)) base_keep else c(base_keep, csl_term)
  
  dd <- df %>% select(all_of(keep)) %>% drop_na()
  if (nrow(dd) < 20) return(NULL)
  
  rhs <- "x07_dfa_cpue_int_spr_jun_hw + x15_dfa_sleeper_shark_bsai_pred_ak + x15_issl_z"
  if (!is.null(csl_term)) rhs <- paste(rhs, "+", csl_term)
  
  model_txt <- paste0(
    "x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus\n",
    "x16_sar ~ ", rhs, "\n"
  )
  
  fit <- tryCatch(sem(model_txt, data = dd, std.lv = TRUE, missing = "ML", warn = FALSE),
                  error = function(e) NULL)
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  
  r2 <- inspect(fit, "r2")[["x16_sar"]]
  fm <- fitMeasures(fit, c("aic","bic","cfi","rmsea"))
  pe <- parameterEstimates(fit)
  
  est <- NA_real_; p <- NA_real_
  if (!is.null(csl_term)) {
    row <- pe %>% filter(lhs == "x16_sar", op == "~", rhs == csl_term)
    if (nrow(row) > 0) {
      est <- row$est[1]
      p <- row$pvalue[1]
    }
  }
  
  tibble(
    csl_term = ifelse(is.null(csl_term), "BASE_NO_CSL", csl_term),
    n = nrow(dd),
    sar_r2 = as.numeric(r2),
    aic = as.numeric(fm["aic"]),
    bic = as.numeric(fm["bic"]),
    cfi = as.numeric(fm["cfi"]),
    rmsea = as.numeric(fm["rmsea"]),
    csl_est = est,
    csl_p = p
  )
}

# ---------------------------
# 6) Fit base + estuary variants
# ---------------------------
cand_terms <- names(csl_lag2_for_sar) %>% setdiff("year")

res_base <- fit_one(d0, NULL)
res_vars <- map_dfr(cand_terms, ~fit_one(d0, .x))
res <- bind_rows(res_base, res_vars)

if (nrow(res) == 0) stop("No models converged.")

base <- res %>% filter(csl_term == "BASE_NO_CSL") %>% slice(1)

res_cmp <- res %>%
  mutate(
    delta_r2_vs_base = sar_r2 - base$sar_r2,
    delta_aic_vs_base = aic - base$aic,
    delta_bic_vs_base = bic - base$bic
  ) %>%
  arrange(desc(delta_r2_vs_base), delta_aic_vs_base)

write.csv(res_cmp, file.path(out_dir, "csl_estuary_lag2_model_comparison.csv"), row.names = FALSE)

best_aug <- res_cmp %>%
  filter(csl_term != "BASE_NO_CSL") %>%
  arrange(desc(delta_r2_vs_base), delta_aic_vs_base) %>%
  slice(1)

# ---------------------------
# 7) Prediction plot for base vs TWO selected CSL models
#    (keep prior script unchanged up to here)
# ---------------------------

refit_preds <- function(df, csl_term = NULL) {
  base_keep <- c("year","x07_dfa_cpue_int_spr_jun_hw","x09_dfa_hake_age5plus","x16_sar",
                 "x15_dfa_sleeper_shark_bsai_pred_ak","x15_issl_z")
  keep <- if (is.null(csl_term)) base_keep else c(base_keep, csl_term)
  
  dd <- df %>% select(all_of(keep)) %>% drop_na()
  
  rhs <- "x07_dfa_cpue_int_spr_jun_hw + x15_dfa_sleeper_shark_bsai_pred_ak + x15_issl_z"
  if (!is.null(csl_term)) rhs <- paste(rhs, "+", csl_term)
  
  model_txt <- paste0(
    "x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus\n",
    "x16_sar ~ ", rhs, "\n"
  )
  
  fit <- sem(model_txt, data = dd, std.lv = TRUE, missing = "ML", warn = FALSE)
  pred <- lavPredictY(fit) %>% as.data.frame()
  
  dd %>% transmute(year, sar_obs = x16_sar, sar_hat = pred$x16_sar)
}

# base + best index (from res_cmp) + forced eulachon-heavy model
pb <- refit_preds(d0, NULL) %>% 
  rename(sar_hat_base = sar_hat)

pp_best <- refit_preds(d0, best_aug$csl_term) %>% 
  rename(sar_hat_best = sar_hat)

pp_eul <- refit_preds(d0, "est_risk_eul_heavy") %>% 
  rename(sar_hat_est_risk_eul_heavy = sar_hat)

# combine by year
pred_compare <- pb %>%
  select(year, sar_obs, sar_hat_base) %>%
  left_join(pp_best %>% select(year, sar_hat_best), by = "year") %>%
  left_join(pp_eul %>% select(year, sar_hat_est_risk_eul_heavy), by = "year")

# rename best column to explicit term name
best_colname <- paste0("sar_hat_", best_aug$csl_term)
names(pred_compare)[names(pred_compare) == "sar_hat_best"] <- best_colname

write.csv(pred_compare,
          file.path(out_dir, "csl_estuary_lag2_base_best_and_eul_predictions.csv"),
          row.names = FALSE)

# plot all lines
plot_df <- pred_compare %>%
  pivot_longer(
    cols = -year,
    names_to = "series",
    values_to = "value"
  )

p <- ggplot(plot_df, aes(year, value, color = series, linetype = series)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray65") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  theme_bw() +
  labs(
    title = "SAR observed vs predicted: base + best CSL + est_risk_eul_heavy",
    subtitle = paste0("Best term: ", best_aug$csl_term,
                      " | ΔR2=", round(best_aug$delta_r2_vs_base, 3),
                      " | ΔAIC=", round(best_aug$delta_aic_vs_base, 2)),
    x = "Year", y = "x16_sar", color = NULL, linetype = NULL
  )
p
ggsave(file.path(out_dir, "csl_estuary_lag2_base_best_and_eul_timeseries.png"),
       p, width = 11, height = 6, dpi = 170)

# ---------------------------
# 8) Add prediction columns to guild_dfas1_24yr.csv and re-save
# ---------------------------
guild_updated <- guild %>%
  left_join(
    csl_lag2_for_sar %>%
      transmute(
        year,
        x15_icsl_cslchin_ratio = idxE10_csl_chinook_ratio,
        x15_icsl_risk_eulachon75_shad25 = est_risk_eul_heavy
      ),
    by = "year"
  ) 

# preserve your exact input file path and overwrite it intentionally
write.csv(guild_updated, "copilot/outputs_6/guild_dfas1_24yr_ishark_issl_icsl.csv", row.names = FALSE)
write.csv(guild_updated, "data_Lisa/guild_dfas1_24yr_ishark_issl_icsl.csv", row.names = FALSE)

# also write a safety copy
# write.csv(guild_updated,
#           file.path(out_dir, "guild_dfas1_24yr_with_csl_predictions.csv"),
#           row.names = FALSE)

# ---------------------------
# 9) Console summary
# ---------------------------
cat("\n=== Base model ===\n")
print(base)

cat("\n=== Top 12 estuary-only lag2 variants ===\n")
print(head(res_cmp, 12))

cat("\n=== Best augmented model ===\n")
print(best_aug)

cat("\nAdded prediction columns:\n")
print(names(pred_compare))

message("\nWrote:")
message(" - ", file.path(out_dir, "csl_estuary_lag2_base_best_and_eul_predictions.csv"))
message(" - ", file.path(out_dir, "csl_estuary_lag2_base_best_and_eul_timeseries.png"))
message(" - ", path_guild, " (updated)")
message(" - ", file.path(out_dir, "guild_dfas1_24yr_with_csl_predictions.csv"))
