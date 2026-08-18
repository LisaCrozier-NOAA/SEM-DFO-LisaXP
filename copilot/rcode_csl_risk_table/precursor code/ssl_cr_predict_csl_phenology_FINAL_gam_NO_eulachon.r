

#GAM for csl phenology: Refactored R Script (Without Eulachon)

# ==============================================================================
# CSL Reconstructed Weekly Time Series (1998–2024, Weeks 10–26)
# Model: GAM with s(week) + offset(log_annual_scale)
# ==============================================================================

library(tidyverse)
library(mgcv)

output_dir <- "copilot/outputs_csl_cr"

# -----------------------------------------------------------------------------
# 1. Read and Clean Input Data
# -----------------------------------------------------------------------------

# CSL Weekly Census Data
week_all_census <- read.csv(file.path(output_dir, "csl_week_all_census.csv"), row.names = NULL) %>%
  filter(year >= 1998, week >= 10, week <= 26)

# CSL Annual Baseline Data (1998–2024)
csl_annual <- read.csv(file.path(output_dir, "csl_annual_baseline_1976_2024.csv"), row.names = NULL) %>%
  filter(year >= 1998) %>%
  select(year, csl_annual_mean)

# -----------------------------------------------------------------------------
# 2. Prep Annual Scaling Term
# -----------------------------------------------------------------------------

csl_annual_prep <- csl_annual %>%
  mutate(
    # Fallback to mean across all years if annual mean is zero or missing
    csl_annual_mean_clean = if_else(is.na(csl_annual_mean) | csl_annual_mean <= 0, 
                                    mean(csl_annual_mean, na.rm = TRUE), 
                                    csl_annual_mean),
    log_annual_scale = log(csl_annual_mean_clean)
  )

# -----------------------------------------------------------------------------
# 3. Construct Complete Master Grid (1998–2024, Weeks 10–26)
# -----------------------------------------------------------------------------

master_grid <- expand_grid(
  year = 1998:2024,
  week = 10:26
) %>%
  left_join(week_all_census, by = c("year", "week")) %>%
  left_join(csl_annual_prep, by = "year") %>%
  arrange(year, week) %>%
  mutate(
    # Continuous time index for plotting multi-year series
    time_index = year + (week - 10) / 17
  )

# -----------------------------------------------------------------------------
# 4. Fit Simplified GAM Model (Training on 2011–2024 Observed CSL)
# -----------------------------------------------------------------------------

gam_train_data <- master_grid %>%
  filter(
    year >= 2011, year <= 2024,
    !is.na(csl_nonpup_total_emb_lisa)
  )

csl_gam <- gam(
  csl_nonpup_total_emb_lisa ~ s(week, k = 5) + offset(log_annual_scale),
  data = gam_train_data,
  family = quasipoisson(link = "log"),
  method = "REML"
)

cat("--- GAM Summary (s(week) only) ---\n")
print(summary(csl_gam))

# -----------------------------------------------------------------------------
# 5. Predict and Reconstruct Full Weekly Time Series (1998–2024)
# -----------------------------------------------------------------------------

gam_preds <- predict(csl_gam, newdata = master_grid, type = "link", se.fit = TRUE)

master_reconstructed <- master_grid %>%
  mutate(
    fit_link     = gam_preds$fit,
    se_link      = gam_preds$se.fit,
    
    # Back-transform from log link to count scale
    csl_gam_pred = exp(fit_link),
    csl_lwr_95   = exp(fit_link - 1.96 * se_link),
    csl_upr_95   = exp(fit_link + 1.96 * se_link),
    
    # Hierarchy: Always use observed if present, otherwise use GAM prediction
    data_source  = if_else(!is.na(csl_nonpup_total_emb_lisa), "Observed", "Predicted"),
    csl_final    = if_else(!is.na(csl_nonpup_total_emb_lisa), 
                           as.numeric(csl_nonpup_total_emb_lisa), 
                           csl_gam_pred)
  )

# -----------------------------------------------------------------------------
# 6. Diagnostic Plots
# -----------------------------------------------------------------------------

# Plot 1: Faceted View for Training Years (2011–2024)
p_recent <- ggplot(master_reconstructed %>% filter(year >= 2011), aes(x = week)) +
  geom_ribbon(aes(ymin = csl_lwr_95, ymax = csl_upr_95), fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = csl_gam_pred, color = "GAM Model Fit"), linewidth = 0.8) +
  geom_point(aes(y = csl_nonpup_total_emb_lisa, color = "Observed Data"), size = 2, na.rm = TRUE) +
  facet_wrap(~ year, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("Observed Data" = "black", "GAM Model Fit" = "steelblue")) +
  labs(
    title = "Weekly CSL Phenology Model Fit: s(week) + offset(log_annual_scale)",
    subtitle = "Fitted on observed 2011–2024 weekly counts",
    x = "Week (10–26)",
    y = "CSL Non-Pup Count",
    color = "Legend"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_recent)

# Plot 2: Full Time Series (1998–2024)
p_full_series <- ggplot(master_reconstructed, aes(x = time_index)) +
  geom_ribbon(aes(ymin = csl_lwr_95, ymax = csl_upr_95), fill = "steelblue", alpha = 0.15) +
  geom_line(aes(y = csl_final), color = "gray40", linewidth = 0.5) +
  geom_point(aes(y = csl_final, color = data_source), size = 1.2) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "firebrick")) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 2)) +
  labs(
    title = "Reconstructed CSL Weekly Time Series (1998–2024, Weeks 10–26)",
    subtitle = "Black = Observed surveys | Red = GAM predicted values",
    x = "Year",
    y = "CSL Count",
    color = "Data Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_full_series)

# -----------------------------------------------------------------------------
# 7. Save Output Files--------
# -----------------------------------------------------------------------------

write.csv(master_reconstructed, file.path(output_dir, "csl_reconstructed_weekly_no_eulachon_1998_2024.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "csl_weekly_fit_no_eulachon_2011_2024.png"), p_recent, width = 11, height = 8)
ggsave(file.path(output_dir, "csl_weekly_reconstructed_no_eulachon_1998_2024.png"), p_full_series, width = 12, height = 6)


# -----------------------------------------------------------------------------
# 8. GAM with and without eulachon--------------
# -----------------------------------------------------------------------------
# ==============================================================================
# Model Comparison: With Eulachon vs. Without Eulachon
# ==============================================================================
master_grid <- read.csv(file.path(output_dir,"eulachon_csl_weekly_reconstruct_master_grid_1998_2024.csv"))


# Ensure both models are fitted on the exact same dataset for fair comparison
train_subset <- master_grid %>%
  filter(
    year >= 2011, year <= 2024,
    !is.na(csl_nonpup_total_emb_lisa),
    !is.na(eulachon_reconstructed)
  )

# Model 1: With Eulachon
gam_with_eul <- gam(
  csl_nonpup_total_emb_lisa ~ s(week, k = 5) + 
    s(eulachon_reconstructed, k = 5) + 
    offset(log_annual_scale),
  data = train_subset,
  family = quasipoisson(link = "log"),
  method = "REML"
)

# Model 2: Without Eulachon (Week Only)
gam_no_eul <- gam(
  csl_nonpup_total_emb_lisa ~ s(week, k = 5) + 
    offset(log_annual_scale),
  data = train_subset,
  family = quasipoisson(link = "log"),
  method = "REML"
)

# -----------------------------------------------------------------------------
# 1. Compare AIC / QAIC
# Note: For QuasiPoisson models, scale parameter > 1 means standard AIC is 
# unadjusted for overdispersion. We compute QAIC (Quasi-AIC).
# -----------------------------------------------------------------------------

dispersion_est <- summary(gam_with_eul)$dispersion

# QAIC function: -2*logLik / dispersion + 2*edf
get_qaic <- function(model, scale) {
  -2 * logLik(model)[1] / scale + 2 * sum(model$edf)
}

qaic_with <- get_qaic(gam_with_eul, dispersion_est)
qaic_no   <- get_qaic(gam_no_eul, dispersion_est)

# -----------------------------------------------------------------------------
# 2. Extract Key Metrics into Comparison Table
# -----------------------------------------------------------------------------

model_comp <- tibble(
  Model = c("With Eulachon", "Without Eulachon (Week Only)"),
  Formula = c(
    "s(week) + s(eulachon) + offset(log_scale)",
    "s(week) + offset(log_scale)"
  ),
  EDF = c(
    sum(gam_with_eul$edf),
    sum(gam_no_eul$edf)
  ),
  Deviance_Explained = c(
    summary(gam_with_eul)$dev.expl * 100,
    summary(gam_no_eul)$dev.expl * 100
  ),
  Adj_R2 = c(
    summary(gam_with_eul)$r.sq,
    summary(gam_no_eul)$r.sq
  ),
  Scale_Dispersion = c(
    summary(gam_with_eul)$scale,
    summary(gam_no_eul)$scale
  ),
  REML_Score = c(
    gam_with_eul$gcv.ubre,
    gam_no_eul$gcv.ubre
  ),
  QAIC = c(qaic_with, qaic_no)
) %>%
  mutate(
    Delta_QAIC = QAIC - min(QAIC),
    Deviance_Explained = sprintf("%.1f%%", Deviance_Explained),
    Adj_R2 = round(Adj_R2, 3),
    EDF = round(EDF, 2),
    QAIC = round(QAIC, 2),
    Delta_QAIC = round(Delta_QAIC, 2)
  )

# -----------------------------------------------------------------------------
# 3. Print Comparison Summary & Likelihood Ratio Test
# -----------------------------------------------------------------------------

cat("\n========================================================\n")
cat("                GAM MODEL FIT COMPARISON                \n")
cat("========================================================\n\n")

print(as.data.frame(model_comp))

cat("\n--- Formal Model Comparison (ANOVA / Likelihood Ratio Test) ---\n")
# F-test is appropriate for QuasiPoisson families to account for overdispersion
anova_res <- anova(gam_no_eul, gam_with_eul, test = "F")
print(anova_res)

# Save comparison summary to CSV
write.csv(model_comp, file.path(output_dir, "csl_predict_gam_w_wo_eulachon_model_comparison.csv"), row.names = FALSE)


# print(as.data.frame(model_comp))
#                         Model                                   Formula  EDF Deviance_Explained Adj_R2 Scale_Dispersion REML_Score QAIC Delta_QAIC
# 1                With Eulachon s(week) + s(eulachon) + offset(log_scale) 8.08              61.8%  0.722         126.1755   445.3218   NA         NA
# 2 Without Eulachon (Week Only)               s(week) + offset(log_scale) 4.74              58.8%  0.697         129.7445   444.3088   NA         NA


# Analysis of Deviance Table
# 
# Model 1: csl_nonpup_total_emb_lisa ~ s(week, k = 5) + offset(log_annual_scale)
# Model 2: csl_nonpup_total_emb_lisa ~ s(week, k = 5) + s(eulachon_reconstructed, 
#                                                         k = 5) + offset(log_annual_scale)
# Resid. Df Resid. Dev     Df Deviance      F  Pr(>F)  
# 1    148.88      15838                                 
# 2    144.66      14679 4.2187   1159.2 2.1777 0.07087 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1