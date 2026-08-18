#check total_bonn or csl_bonn to predict missing csl_emb data
#DO NOT USE
# QAIC (Annual Mean Scale): 76.56821 
# QAIC (total_bonn Scale):  93.28634 
# QAIC (csl_bonn Scale):  97.36923


library(tidyverse)
library(mgcv)

# -----------------------------------------------------------------------------
# 1. Prepare Dataset with total_bonn Scaling
# -----------------------------------------------------------------------------
week_bonn_prep <- week_all %>%
  filter(month %in% c(4, 5, 6), year >= 2011, year <= 2024) %>%
  mutate(
    year_factor = factor(year),
    
    # Log-transform total_bonn to serve as a proportional magnitude scalar
    # (Adding 1 inside log handles any potential zero values)
    log_total_bonn = log(pmax(1, total_bonn)),
    log_bonn = log(pmax(1, csl_bonn)),
    
    # Ensure Eulachon inputs are present
    eulachon_input = eulachon_ssb_4week_est
  )

# -----------------------------------------------------------------------------
# 2. Fit Unified GAM with Continuous total_bonn Predictor
# -----------------------------------------------------------------------------
# This single model handles all years uniformly:
#  - s(week) shapes the baseline spring curve
#  - s(eulachon_input) shifts the peak timing dynamically
#  - log_bonn scales the overall vertical height of the curve per year

csl_total_bonn_gam <- gam(
  csl_nonpup_total_emb ~ s(week, k = 5) + 
    s(eulachon_input, k = 5) + 
    log_total_bonn,
  data = week_bonn_prep %>% filter(!is.na(csl_nonpup_total_emb) & !is.na(eulachon_input) & !is.na(total_bonn)),
  family = quasipoisson(link = "log"),
  method = "REML"
)


csl_bonn_gam <- gam(
  csl_nonpup_total_emb ~ s(week, k = 5) + 
    s(eulachon_input, k = 5) + 
    log_bonn,
  data = week_bonn_prep %>% filter(!is.na(csl_nonpup_total_emb) & !is.na(eulachon_input) & !is.na(total_bonn)),
  family = quasipoisson(link = "log"),
  method = "REML"
)

# Function to calculate QAIC for quasi-Poisson GAMs
get_qaic <- function(model) {
  # 1. Fit an equivalent Poisson model to get log-likelihood
  poisson_mod <- update(model, family = poisson)
  
  # 2. Extract dispersion parameter (phi) from original quasi-Poisson model
  phi <- summary(model)$dispersion
  
  # 3. Extract log-likelihood and parameter count k
  ll <- as.numeric(logLik(poisson_mod))
  k  <- sum(model$edf) # Effective degrees of freedom
  
  # 4. Compute QAIC
  qaic <- (-2 * ll / phi) + (2 * k)
  return(qaic)
}

# Compare your models using QAIC
qaic_model1 <- get_qaic(csl_scaled_gam)
qaic_model2 <- get_qaic(csl_total_bonn_gam)
qaic_model3 <- get_qaic(csl_bonn_gam)

cat("QAIC (Annual Mean Scale):", qaic_model1, "\n")
cat("QAIC (total_bonn Scale): ", qaic_model2, "\n")
cat("QAIC (csl_bonn Scale): ", qaic_model3, "\n")

# Inspect model summary—check the t-stat / p-value on log_bonn!
summary(csl_bonn_gam)

# -----------------------------------------------------------------------------
# 3. Predict CSL Phenology Across All Years (Including 2019)
# -----------------------------------------------------------------------------
# Because log_bonn is present in all years (including 2019), 
# 2019 now gets its magnitude directly from 2019's total_bonn value!

bonn_preds_link <- predict(
  csl_bonn_gam,
  newdata = week_bonn_prep,
  type = "link",
  se.fit = TRUE
)

week_plot_bonn <- week_bonn_prep %>%
  mutate(
    fit_link     = bonn_preds_link$fit,
    se_link      = bonn_preds_link$se.fit,
    
    # Back-transform to response count scale
    csl_gam_pred = exp(fit_link),
    csl_lwr_95   = exp(fit_link - 1.96 * se_link),
    csl_upr_95   = exp(fit_link + 1.96 * se_link)
  )

# -----------------------------------------------------------------------------
# 4. Plot Evaluation Across All Years
# -----------------------------------------------------------------------------
ggplot(week_plot_bonn, aes(x = week)) +
  geom_ribbon(aes(ymin = csl_lwr_95, ymax = csl_upr_95), fill = "darkgreen", alpha = 0.2) +
  geom_line(aes(y = csl_gam_pred, color = "total_bonn Scaled GAM"), size = 1) +
  geom_point(aes(y = csl_nonpup_total_emb, color = "Observed Data"), size = 2, na.rm = TRUE) +
  facet_wrap(~ year, scales = "free_y") +
  scale_color_manual(values = c("Observed Data" = "black", "total_bonn Scaled GAM" = "darkgreen")) +
  labs(
    title = "CSL Phenology Scaled by Continuous `total_bonn` (2011–2024)",
    subtitle = "Single unified model using log(total_bonn) to anchor annual magnitude across all years",
    x = "Week Number (Apr–Jun)",
    y = "CSL Non-Pup Count",
    color = "Series"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

