suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

# =========================================================
# CSL vs residual alignment from top SEM model
# Repo: LisaCrozier-NOAA/SEM-DFO-LisaXP
# =========================================================

out_dir <- "copilot/outputs_6"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------------------------
# 0) Paths
# ---------------------------
#path_guild <- "data_Lisa/guild.dfasAK.csv"
path_guild<-  "copilot/outputs_6/guild_dfas1_24yr.csv"

path_csl   <- file.path(out_dir, "master_csl_predation_risk_table_5ocean_prey_1998_2024.csv")

if (!file.exists(path_guild)) stop("Missing: ", path_guild)
if (!file.exists(path_csl))   stop("Missing: ", path_csl)

# ---------------------------
# 1) Load and prep data
# ---------------------------
guild <- read.csv(path_guild) %>%
  clean_names() %>%
  select(-any_of("x10_harbor_seal_cr_2yr_lead"))

csl <- read.csv(path_csl) %>%
  clean_names()

# ensure year
if (!"year" %in% names(guild)) stop("guild missing year")
if (!"year" %in% names(csl)) stop("csl table missing year")

# ---------------------------
# 2) Refit top SEM model (1998-2021)
# top model: x15_dfa_sleeper_shark_bsai_pred_ak + x15_issl_z
# ---------------------------
req_sem <- c(
  "year",
  "x07_dfa_cpue_int_spr_jun_hw",
  "x09_dfa_hake_age5plus",
  "x16_sar",
  "x15_dfa_sleeper_shark_bsai_pred_ak",
  "x15_issl_z"
)
miss_sem <- setdiff(req_sem, names(guild))
if (length(miss_sem) > 0) {
  stop("guild missing SEM columns: ", paste(miss_sem, collapse = ", "))
}

df_sem <- guild %>%
  filter(year >= 1998, year <= 2021) %>%
  select(all_of(req_sem)) %>%
  drop_na()

if (nrow(df_sem) < 20) stop("Too few complete rows for SEM fit.")

sem_text <- '
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + x15_dfa_sleeper_shark_bsai_pred_ak + x15_issl_z
'

fit <- sem(sem_text, data = df_sem, std.lv = TRUE, missing = "ML", warn = FALSE)
if (!lavInspect(fit, "converged")) stop("SEM did not converge.")

# Extract fitted values for x16_sar equation and residuals
pred_all <- lavPredictY(fit) %>% as.data.frame()
if (!"x16_sar" %in% names(pred_all)) stop("Could not extract predicted x16_sar from lavaan fit.")

resid_df <- df_sem %>%
  mutate(
    sar_obs = x16_sar,
    sar_hat = pred_all$x16_sar,
    sar_resid = sar_obs - sar_hat,
    sar_resid_sign = sign(sar_resid)
  ) %>%
  select(year, sar_obs, sar_hat, sar_resid, sar_resid_sign)

# ---------------------------
# 3) Join CSL table + quantitative features
# ---------------------------
d <- resid_df %>%
  left_join(csl, by = "year")

# Required CSL quantitative fields to start
req_csl <- c("csl_during_chinook", "ocean_buffer_score")
miss_csl <- setdiff(req_csl, names(d))
if (length(miss_csl) > 0) stop("Joined data missing CSL fields: ", paste(miss_csl, collapse = ", "))

# helper z
safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - m) / s)
}

# build quantitative candidates
d <- d %>%
  mutate(
    ocean_risk_factor = 5 - ocean_buffer_score,
    
    # raw or lightly transformed components
    csl_log = log1p(csl_during_chinook),
    eulachon_log = if ("eulachon_during_chinook" %in% names(.)) log1p(eulachon_during_chinook) else NA_real_,
    shad_log     = if ("shad_during_chinook" %in% names(.)) log1p(shad_during_chinook) else NA_real_
  )

# standardized components
d <- d %>%
  mutate(
    csl_z = safe_z(csl_during_chinook),
    csl_log_z = safe_z(csl_log),
    ocean_buffer_z = safe_z(ocean_buffer_score),
    ocean_risk_z = safe_z(ocean_risk_factor),
    eulachon_z = if ("eulachon_during_chinook" %in% names(.)) safe_z(eulachon_during_chinook) else NA_real_,
    shad_z     = if ("shad_during_chinook" %in% names(.)) safe_z(shad_during_chinook) else NA_real_
  )

# Optional 5-ocean-prey quantitative composite if columns exist
ocean5_cols <- intersect(
  c("x09_dfa_hake_age5plus_doug", "x05_dfa_abund_sardine_doug", "sitka_herring_e_go_a_doug",
    "x05_herring_ncc_doug", "x05_anchovy_gam_doug"),
  names(d)
)

if (length(ocean5_cols) > 0) {
  d <- d %>%
    mutate(
      ocean5_mean = rowMeans(select(., all_of(ocean5_cols)), na.rm = TRUE),
      ocean5_mean_z = safe_z(ocean5_mean)
    )
} else {
  d$ocean5_mean <- NA_real_
  d$ocean5_mean_z <- NA_real_
}

# Candidate CSL indices (higher = higher predicted predation pressure)
d <- d %>%
  mutate(
    idx1_csl_only          = csl_z,
    idx2_csl_oceanrisk     = csl_z + ocean_risk_z,
    idx3_csl_minus_buffer  = csl_z - ocean_buffer_z,
    idx4_csl_buffer_mult   = csl_z * (1 + ocean_risk_z),
    idx5_csl_estuary_prey  = csl_z - 0.5*eulachon_z - 0.5*shad_z,
    idx6_full_quant        = csl_z + ocean_risk_z - 0.5*eulachon_z - 0.5*shad_z,
    idx7_with_ocean5       = csl_z - 0.5*ocean5_mean_z
  )

# ---------------------------
# 4) Evaluate alignment with residuals
# ---------------------------
candidate_cols <- c(
  "idx1_csl_only",
  "idx2_csl_oceanrisk",
  "idx3_csl_minus_buffer",
  "idx4_csl_buffer_mult",
  "idx5_csl_estuary_prey",
  "idx6_full_quant",
  "idx7_with_ocean5",
  "csl_z",
  "ocean_risk_z",
  "ocean_buffer_z"
)

candidate_cols <- candidate_cols[candidate_cols %in% names(d)]

score_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8) {
    return(tibble(
      n = sum(ok), spearman = NA_real_, pearson = NA_real_,
      lm_slope = NA_real_, lm_p = NA_real_,
      sign_match = NA_real_, sign_match_nonzero = NA_real_
    ))
  }
  
  xx <- x[ok]
  yy <- y[ok]
  
  # sign match: predicted risk should align with negative residual?
  # If higher risk => more negative residual, then compare sign(x) to sign(-y)
  sm_all <- mean(sign(xx) == sign(-yy), na.rm = TRUE)
  
  nz <- (xx != 0 & yy != 0)
  sm_nz <- if (sum(nz) > 0) mean(sign(xx[nz]) == sign(-yy[nz]), na.rm = TRUE) else NA_real_
  
  fit_lm <- lm(yy ~ xx)
  co <- summary(fit_lm)$coefficients
  
  tibble(
    n = sum(ok),
    spearman = suppressWarnings(cor(xx, yy, method = "spearman")),
    pearson  = suppressWarnings(cor(xx, yy, method = "pearson")),
    lm_slope = unname(co["xx", "Estimate"]),
    lm_p     = unname(co["xx", "Pr(>|t|)"]),
    sign_match = sm_all,
    sign_match_nonzero = sm_nz
  )
}

eval_tbl <- map_dfr(candidate_cols, function(cc) {
  score_one(d[[cc]], d$sar_resid) %>% mutate(candidate = cc, .before = 1)
}) %>%
  mutate(
    # for interpretation where higher risk should mean lower residual, negative corr is "better"
    expected_direction = "negative_correlation_with_residual",
    score_rank = rank(-sign_match_nonzero, ties.method = "min")
  ) %>%
  arrange(desc(sign_match_nonzero), spearman)

# ---------------------------
# 5) Save outputs
# ---------------------------
write.csv(eval_tbl, file.path(out_dir, "csl_residual_alignment_candidate_scores.csv"), row.names = FALSE)
write.csv(d, file.path(out_dir, "csl_residual_alignment_joined_timeseries.csv"), row.names = FALSE)

fit_measures <- fitMeasures(fit, c("aic", "bic", "cfi", "rmsea"))
fit_tbl <- tibble(
  metric = names(fit_measures),
  value = as.numeric(fit_measures)
)
write.csv(fit_tbl, file.path(out_dir, "top_model_refit_fitmeasures.csv"), row.names = FALSE)

# ---------------------------
# 6) Plots
# ---------------------------
# choose best candidate by sign match nonzero
best_cand <- eval_tbl$candidate[1]

p_ts <- d %>%
  select(year, sar_resid, all_of(best_cand)) %>%
  pivot_longer(cols = -year, names_to = "series", values_to = "value") %>%
  ggplot(aes(year, value, color = series)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray55") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.6) +
  theme_bw() +
  labs(
    title = paste0("Top-model SAR residual vs best CSL candidate: ", best_cand),
    subtitle = "Interpretation target: higher CSL risk aligns with more negative residuals",
    x = "Year", y = "Value", color = NULL
  )

ggsave(file.path(out_dir, "csl_best_candidate_vs_residual_timeseries.png"), p_ts, width = 10, height = 5.6, dpi = 170)

p_sc <- ggplot(d, aes(x = .data[[best_cand]], y = sar_resid)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray60") +
  geom_vline(xintercept = 0, linetype = 2, color = "gray60") +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue3") +
  theme_bw() +
  labs(
    title = paste0("Residual alignment scatter: ", best_cand),
    subtitle = "Expected slope/correlation is negative if risk index is directionally consistent",
    x = best_cand, y = "SAR residual (obs - fitted)"
  )

ggsave(file.path(out_dir, "csl_best_candidate_vs_residual_scatter.png"), p_sc, width = 7, height = 5.4, dpi = 170)

p_rank <- eval_tbl %>%
  mutate(candidate = fct_reorder(candidate, sign_match_nonzero, .desc = TRUE)) %>%
  ggplot(aes(candidate, sign_match_nonzero)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "CSL candidate indices ranked by sign-match with -residual",
    x = NULL, y = "Sign match (nonzero years)"
  )

ggsave(file.path(out_dir, "csl_candidate_signmatch_rank.png"), p_rank, width = 8.5, height = 5.8, dpi = 170)

# ---------------------------
# 7) Console summary
# ---------------------------
cat("\n=== Top-model refit fit measures ===\n")
print(fit_tbl)

cat("\n=== CSL candidate ranking (top 10) ===\n")
print(head(eval_tbl, 10))

cat("\nBest candidate: ", best_cand, "\n")
message("Outputs written to: ", out_dir)


#examine r2, add best csl index to SEM----------
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

out_dir <- "copilot/outputs_6"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------- paths --------
#path_guild <- "data_Lisa/guild.dfasAK.csv"
path_guild<-  "copilot/outputs_6/guild_dfas1_24yr.csv"
path_csl   <- file.path(out_dir, "master_csl_predation_risk_table_5ocean_prey_1998_2024.csv")

# -------- load --------
guild <- read.csv(path_guild) %>% clean_names() %>%
  select(-any_of("x10_harbor_seal_cr_2yr_lead"))
csl <- read.csv(path_csl) %>% clean_names()

# -------- rebuild idx4 exactly as prior script --------
safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - m)/s)
}

tmp <- csl %>%
  mutate(
    ocean_risk_factor = 5 - ocean_buffer_score,
    csl_z = safe_z(csl_during_chinook),
    ocean_risk_z = safe_z(ocean_risk_factor),
    idx4_csl_buffer_mult = csl_z * (1 + ocean_risk_z)
  ) %>%
  select(year, idx4_csl_buffer_mult)

dat <- guild %>%
  left_join(tmp, by = "year") %>%
  filter(year >= 1998, year <= 2021) %>%
  select(
    year,
    x07_dfa_cpue_int_spr_jun_hw,
    x09_dfa_hake_age5plus,
    x16_sar,
    x15_dfa_sleeper_shark_bsai_pred_ak,
    x15_issl_z,
    idx4_csl_buffer_mult
  ) %>%
  drop_na()

if (nrow(dat) < 20) stop("Too few complete rows after join.")

# -------- model specs --------
model_base <- '
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + x15_dfa_sleeper_shark_bsai_pred_ak + x15_issl_z
'

model_plus <- '
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + x15_dfa_sleeper_shark_bsai_pred_ak + x15_issl_z + idx4_csl_buffer_mult
'

fit_base <- sem(model_base, data = dat, std.lv = TRUE, missing = "ML", warn = FALSE)
fit_plus <- sem(model_plus, data = dat, std.lv = TRUE, missing = "ML", warn = FALSE)

if (!lavInspect(fit_base, "converged")) stop("Base model did not converge.")
if (!lavInspect(fit_plus, "converged")) stop("Augmented model did not converge.")

# -------- R2 for SAR --------
r2_base <- inspect(fit_base, "r2")[["x16_sar"]]
r2_plus <- inspect(fit_plus, "r2")[["x16_sar"]]
delta_r2 <- r2_plus - r2_base

# -------- fit metrics --------
fm_base <- fitMeasures(fit_base, c("aic","bic","cfi","rmsea"))
fm_plus <- fitMeasures(fit_plus, c("aic","bic","cfi","rmsea"))

cmp <- tibble(
  model = c("base", "plus_csl_idx4"),
  n = nrow(dat),
  sar_r2 = c(r2_base, r2_plus),
  aic = c(fm_base["aic"], fm_plus["aic"]),
  bic = c(fm_base["bic"], fm_plus["bic"]),
  cfi = c(fm_base["cfi"], fm_plus["cfi"]),
  rmsea = c(fm_base["rmsea"], fm_plus["rmsea"])
) %>%
  mutate(
    delta_r2_vs_base = sar_r2 - sar_r2[model == "base"],
    delta_aic_vs_base = aic - aic[model == "base"],
    delta_bic_vs_base = bic - bic[model == "base"]
  )

write.csv(cmp, file.path(out_dir, "sar_r2_base_vs_plus_csl.csv"), row.names = FALSE)

# -------- observed vs predicted --------
pred_base <- lavPredictY(fit_base) %>% as.data.frame()
pred_plus <- lavPredictY(fit_plus) %>% as.data.frame()

pred_df <- dat %>%
  transmute(
    year,
    sar_obs = x16_sar,
    sar_hat_base = pred_base$x16_sar,
    sar_hat_plus = pred_plus$x16_sar
  )

write.csv(pred_df, file.path(out_dir, "sar_observed_vs_pred_base_plus.csv"), row.names = FALSE)

p <- pred_df %>%
  pivot_longer(-year, names_to = "series", values_to = "value") %>%
  ggplot(aes(year, value, color = series)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  theme_bw() +
  labs(
    title = "SAR observed vs predicted: base vs +CSL index",
    subtitle = sprintf("R2 base=%.3f, R2 plus=%.3f, ΔR2=%.3f", r2_base, r2_plus, delta_r2),
    x = "Year", y = "x16_sar", color = NULL
  )

ggsave(file.path(out_dir, "sar_base_vs_plus_timeseries.png"), p, width = 10, height = 5.6, dpi = 170)

cat("\n=== SAR prediction comparison ===\n")
print(cmp)
