suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

# -----------------------------------------------------------------------------
# 1. Define All SSL Candidate Variables
# -----------------------------------------------------------------------------
ssl_cands <- c(
  "ssl_model_eric",
  "ssl_west_pup_pred",
  "ssl_seak_pup_pred",
  "ssl_west_np_pred",
  "ssl_east_pup_pred",
  "ssl_east_np_pred",
  "ssl_cent_pup_pred",
  "ssl_cent_np_pred"
)

# Filter candidate list to only those present in ssl.dat
ssl_cands <- ssl_cands[ssl_cands %in% names(ssl.dat)]

if (length(ssl_cands) == 0) {
  stop("None of the specified SSL candidate columns were found in ssl.dat.")
}

# -----------------------------------------------------------------------------
# 2. Loop Across SSL Candidates & Fit Linear Interaction Models
# -----------------------------------------------------------------------------
results_list <- list()

for (i in seq_along(ssl_cands)) {
  ssl_var <- ssl_cands[i]
  
  # Select & Clean Input Data
  d0 <- ssl.dat %>%
    select(
      year,
      sar_ishark_resid,
      ssl_raw = all_of(ssl_var),
      sst_raw = sst_wgoa_coastwatch_junjulaug,
      f1_raw  = capelin_avg,
      f2_raw  = x13_stka_herr_matbiom
    ) %>%
    filter(if_all(everything(), ~ is.finite(.x))) %>%
    filter(year >= 1998, year <= 2021)
  
  if (nrow(d0) < 10) next
  
  # Scale Variables (Z-scores for standardized betas)
  d <- d0 %>%
    mutate(
      y        = as.vector(scale(sar_ishark_resid)),
      ssl      = as.vector(scale(ssl_raw)),
      sst      = as.vector(scale(sst_raw)),
      f_cap    = as.vector(scale(f1_raw)),
      f_herr   = as.vector(scale(f2_raw)),
      
      # Explicit Linear Interaction Terms
      ssl_x_sst    = ssl * sst,
      ssl_x_cap    = ssl * f_cap,
      ssl_x_herr   = ssl * f_herr
    )
  
  # Fit Linear Interaction Model
  mod <- lm(y ~ ssl + ssl_x_sst + ssl_x_cap + ssl_x_herr, data = d)
  
  sm <- summary(mod)
  td <- broom::tidy(mod)
  
  get_term <- function(term_name, field) {
    out <- td %>% filter(term == term_name) %>% pull(!!sym(field))
    if (length(out) == 0) NA_real_ else out[1]
  }
  
  # Extract Coefficients & Predictions
  results_list[[i]] <- tibble(
    ssl_variable    = ssl_var,
    n               = nrow(d),
    r2              = sm$r.squared,
    adj_r2          = sm$adj.r.squared,
    aic             = AIC(mod),
    
    # Main SSL Effect (Expect Negative)
    b_ssl           = get_term("ssl", "estimate"),
    p_ssl           = get_term("ssl", "p.value"),
    
    # SSL x SST Interaction (Negative = Warm SST increases SSL predation)
    b_ssl_x_sst     = get_term("ssl_x_sst", "estimate"),
    p_ssl_x_sst     = get_term("ssl_x_sst", "p.value"),
    
    # SSL x Capelin Interaction (Positive = Capelin buffers salmon)
    b_ssl_x_cap     = get_term("ssl_x_cap", "estimate"),
    p_ssl_x_cap     = get_term("ssl_x_cap", "p.value"),
    
    # SSL x Herring Interaction (Positive = Herring buffers salmon)
    b_ssl_x_herr    = get_term("ssl_x_herr", "estimate"),
    p_ssl_x_herr    = get_term("ssl_x_herr", "p.value")
  )
}

res_ssl_interactions <- bind_rows(results_list)

# -----------------------------------------------------------------------------
# 3. Rank & Display Results
# -----------------------------------------------------------------------------
ranked_ssl_models <- res_ssl_interactions %>%
  mutate(
    # Check if interactions match expected biological signs:
    # b_ssl_x_sst < 0 (warmth increases predation)
    # b_ssl_x_herr > 0 OR b_ssl_x_cap > 0 (forage buffers predation)
    expected_signs = (b_ssl_x_sst < 0) & (b_ssl_x_herr > 0 | b_ssl_x_cap > 0)
  ) %>%
  arrange(desc(r2))

cat("\n=== SSL LINEAR INTERACTION MODEL RANKINGS ===\n")
print(
  ranked_ssl_models %>% 
    select(ssl_variable, r2, adj_r2, aic, b_ssl, b_ssl_x_sst, b_ssl_x_herr, b_ssl_x_cap, expected_signs)
)

# Export Summary to CSV
write_csv(ranked_ssl_models, file.path(out_dir, "ssl_linear_interaction_sweep.csv"))