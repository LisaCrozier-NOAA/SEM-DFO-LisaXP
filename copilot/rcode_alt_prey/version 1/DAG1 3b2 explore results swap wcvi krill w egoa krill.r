suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

out_dir <- "copilot/outputs_9"

# -----------------------------------------------------------------------------
# 1. Load Shifted Positive Dataset (+5.0)
# -----------------------------------------------------------------------------
sem_raw <- read.csv(file.path(out_dir, "sem_altprey_data_complete_1998_2021.csv"), row.names = NULL)
names(sem_raw) <- tolower(names(sem_raw))

sem_complete_data <- sem_raw %>%
  mutate(across(where(is.numeric) & !matches("year"), ~ .x + 5.0))

# -----------------------------------------------------------------------------
# 2. Build Interaction Terms for Both Models
# -----------------------------------------------------------------------------
df_comp <- sem_complete_data

# Common Predators & Prey
ak_pred  <- "x11_ssl_seak_pup_pred"
herr     <- "x13_stka_herr_matbiom"
capelin  <- "x13_mid_il_capelin"
egoa_kr  <- "x12_egoa.krill"
wcvi_kr  <- "x12_dfa_biomasseuphshelfsum"

# Common interaction terms
df_comp$ak_int_1 <- df_comp[[ak_pred]] * df_comp[[herr]]
df_comp$ak_int_2 <- df_comp[[ak_pred]] * df_comp[[capelin]]

# Model 39A (Original with eGOA Krill)
df_comp$ak_int_3_egoa <- df_comp[[ak_pred]] * df_comp[[egoa_kr]]

# Model 39B (Variant with WCVI Krill)
df_comp$ak_int_3_wcvi <- df_comp[[ak_pred]] * df_comp[[wcvi_kr]]

# -----------------------------------------------------------------------------
# 3. Fit Both SEM Models
# -----------------------------------------------------------------------------
syntax_egoa <- "
  x07_dfa_cpue_intsprjunhw ~ x09_dfa_hakeage5plus
  x16_sar ~ x07_dfa_cpue_intsprjunhw + x11_ssl_seak_pup_pred + ak_int_1 + ak_int_2 + ak_int_3_egoa
"

syntax_wcvi <- "
  x07_dfa_cpue_intsprjunhw ~ x09_dfa_hakeage5plus
  x16_sar ~ x07_dfa_cpue_intsprjunhw + x11_ssl_seak_pup_pred + ak_int_1 + ak_int_2 + ak_int_3_wcvi
"

fit_egoa <- sem(syntax_egoa, data = df_comp, std.lv = TRUE, missing = "ML", warn = FALSE)
fit_wcvi <- sem(syntax_wcvi, data = df_comp, std.lv = TRUE, missing = "ML", warn = FALSE)

# -----------------------------------------------------------------------------
# 4. Extract & Compare Fit Statistics
# -----------------------------------------------------------------------------
extract_model_stats <- function(fit, label) {
  fm <- fitMeasures(fit, c("aic", "bic", "cfi", "rmsea", "pvalue"))
  r2 <- inspect(fit, "r2")
  
  tibble(
    model       = label,
    aic         = fm[["aic"]],
    bic         = fm[["bic"]],
    cfi         = fm[["cfi"]],
    pvalue      = fm[["pvalue"]],
    cpue_r2     = as.numeric(r2[["x07_dfa_cpue_intsprjunhw"]]),
    sar_r2      = as.numeric(r2[["x16_sar"]])
  )
}

comparison_summary <- bind_rows(
  extract_model_stats(fit_egoa, "Model 39A (eGOA Krill)"),
  extract_model_stats(fit_wcvi, "Model 39B (WCVI Krill)")
) %>%
  mutate(delta_aic = aic - min(aic))

cat("\n====================================================================\n")
cat(" MODEL COMPARISON: eGOA KRILL vs. WCVI KRILL                       \n")
cat("====================================================================\n")

print(comparison_summary %>% as_tibble(), n = Inf)

# -----------------------------------------------------------------------------
# 5. Extract & Compare Parameter Estimates
# -----------------------------------------------------------------------------
extract_params <- function(fit, label) {
  parameterEstimates(fit) %>%
    filter(op == "~") %>%
    select(lhs, term = rhs, est, se, z, pvalue) %>%
    mutate(
      model = label,
      resolved_species = case_when(
        term == "x09_dfa_hakeage5plus" ~ "Hake (Age 5+)",
        term == "x07_dfa_cpue_intsprjunhw" ~ "Early CPUE Link",
        term == "x11_ssl_seak_pup_pred" ~ "Steller Sea Lions (SEAK)",
        term == "ak_int_1" ~ "Sitka Herring Buffer",
        term == "ak_int_2" ~ "Capelin Buffer",
        term == "ak_int_3_egoa" ~ "eGOA Krill Buffer",
        term == "ak_int_3_wcvi" ~ "WCVI Krill Buffer",
        TRUE ~ term
      )
    )
}

comparison_params <- bind_rows(
  extract_params(fit_egoa, "Model 39A (eGOA Krill)"),
  extract_params(fit_wcvi, "Model 39B (WCVI Krill)")
)

cat("\n====================================================================\n")
cat(" PARAMETER ESTIMATES COMPARISON                                     \n")
cat("====================================================================\n")

print(
  comparison_params %>%
    select(model, stage = lhs, resolved_species, est, se, z, pvalue) %>%
    as_tibble(),
  n = Inf
)
