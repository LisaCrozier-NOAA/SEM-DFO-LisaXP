#Load packages---------
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

out_dir <- "copilot/outputs_6"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# index functions-----------------------------------------------------------------------------
# 1. Define Portable Index Functions with Fixed Benchmark Scaling

calc_i_ssl <- function(ssl, sst, herring, capelin, 
                       beta_sst  =  0.25, 
                       beta_herr = -2.00, 
                       beta_cap  = -1.00,
                       # Baseline reference scale (1998-2021)
                       ssl_ref  = list(mean = mean(ssl, na.rm=T), sd = sd(ssl, na.rm=T)),
                       sst_ref  = list(mean = mean(sst, na.rm=T), sd = sd(sst, na.rm=T)),
                       herr_ref = list(mean = mean(herring, na.rm=T), sd = sd(herring, na.rm=T)),
                       cap_ref  = list(mean = mean(capelin, na.rm=T), sd = sd(capelin, na.rm=T))) {
  
  ssl_z  <- (ssl - ssl_ref$mean) / ssl_ref$sd
  sst_z  <- (sst - sst_ref$mean) / sst_ref$sd
  herr_z <- (herring - herr_ref$mean) / herr_ref$sd
  cap_z  <- (capelin - cap_ref$mean) / cap_ref$sd
  
  raw_index <- ssl_z * (1.0 + beta_sst * sst_z + beta_herr * herr_z + beta_cap * cap_z)
  
  # Standardize final output
  as.vector(scale(raw_index))
}

calc_i_shark <- function(shark, temp_raw_Mt, temp_Ot,
                         alpha_ot     = 0.372,
                         Q10          = 2.0,
                         overlap_form = "logistic",
                         shark_ref    = list(mean = mean(shark, na.rm=T), sd = sd(shark, na.rm=T)),
                         temp_ref_Mt  = mean(temp_raw_Mt, na.rm=T)) {
  
  shark_z <- (shark - shark_ref$mean) / shark_ref$sd
  M_t     <- Q10^((temp_raw_Mt - temp_ref_Mt) / 10)
  
  if (overlap_form == "logistic") {
    O_t <- plogis(alpha_ot * temp_Ot) * 2
  } else if (overlap_form == "linear") {
    O_t <- pmax(0.1, 1.0 + alpha_ot * temp_Ot)
  } else {
    O_t <- rep(1.0, length(temp_Ot))
  }
  
  raw_index <- shark_z * M_t * O_t
  
  # Standardize final output
  as.vector(scale(raw_index))
}

# Load data-----------------------------------------------------------------------------
# 2. Load Raw Component Data & Compute Both Indices
# -----------------------------------------------------------------------------

# A. Shark & Climate Data
shark_dat <- read.csv("data_Lisa/shark_wide.csv") %>% 
  clean_names() %>%
  select(year, goa_pacific_sleeper_shark)

prey_clim_dat <- read.csv("data_Lisa/goa_prey_clim_raw_trends_avg.csv", row.names = 1) %>% 
  clean_names()

# B. Sea Lion Counts
ssl_dat <- read.csv("data_Lisa/ssl.dat.csv") %>% 
  clean_names() %>% 
  select(year, ssl_seak_pup_pred)

# C. Assemble Combined Index Input Dataframe
index_inputs <- shark_dat %>%
  left_join(prey_clim_dat %>% select(year, sst_wgoa_coastwatch_junjulaug, pdo_djf, 
                                     egoa_bio_stka_herr_matbiom, wgoa_bio_mid_il_capelin), by = "year") %>%
  left_join(ssl_dat, by = "year") %>%
  mutate(
    # Compute single-year shark index
    x15_ishark_z = calc_i_shark(
      shark       = goa_pacific_sleeper_shark,
      temp_raw_Mt = sst_wgoa_coastwatch_junjulaug,
      temp_Ot     = pdo_djf,
      alpha_ot    = 0.372,
      overlap_form = "linear"
    ),
    
    # Compute sea lion salmon hazard index
    x15_issl_z = calc_i_ssl(
      ssl     = ssl_seak_pup_pred,
      sst     = sst_wgoa_coastwatch_junjulaug,
      herring = egoa_bio_stka_herr_matbiom,
      capelin = wgoa_bio_mid_il_capelin,
      beta_sst  =  0.25,
      beta_herr = -2.00,
      beta_cap  = -1.00
    )
  )

# D. Load Guild Master Data & Join Computed Indices
guild.dfas1 <- read.csv("data_Lisa/guild.dfasAK.csv", row.names = NULL) %>% 
  clean_names() %>% 
  select(-any_of("x10_harbor_seal_cr_2yr_lead")) %>%
  left_join(index_inputs %>% select(year, x15_ishark_z, x15_issl_z), by = "year") %>%
  mutate(across(-year, ~ as.vector(scale(.)))) %>%
  mutate(
    # Compute single-year shark index
    x15_ishark_DFAsleeper_z = calc_i_shark(
      shark       = x15_dfa_sleeper_shark_bsai_pred_ak,
      temp_raw_Mt = x21_sst_wgoa_coastwatch_junjulaug,
      temp_Ot     = x21_pdo_djf,
      alpha_ot    = 0.372,
      overlap_form = "linear"
    ))


# Filter dataset to 1998-2021 window FIRST, then keep ONLY columns with all 24 years available
guild_dfas1_24yr <- guild.dfas1 %>%
  filter(year >= 1998, year <= 2021) %>%
  select(where(~ sum(!is.na(.)) >= 24))  # <--- THIS LINE ENFORCES 24-YEAR COMPLETENESS

write_csv(guild_dfas1_24yr, file.path(out_dir, "guild_dfas1_24yr.csv"))


# -----------------------------------------------------------------------------
# 3. Candidate Sweeps Setup
# -----------------------------------------------------------------------------

ak_prey_cands <- names(guild_dfas1_24yr)[grepl("^x21_|^x12_|^x13_|^x14_", names(guild_dfas1_24yr))]
ak_pred_cands <- names(guild_dfas1_24yr)[grepl("^x10_|^x15_", names(guild_dfas1_24yr))]

message("Found ", length(ak_prey_cands), " Prey candidates and ", 
        length(ak_pred_cands), " Predator candidates with complete 24-year series.")



# Helper function to fit lavaan SEM
fit_sem_two_topdown <- function(pred1_col, pred2_col, df_data) {
  df_model <- df_data %>%
    filter(year >= 1998, year <= 2021) %>%
    mutate(
      Pred1 = .data[[pred1_col]],
      Pred2 = .data[[pred2_col]]
    ) %>%
    filter(
      !is.na(x07_dfa_cpue_int_spr_jun_hw),
      !is.na(x09_dfa_hake_age5plus),
      !is.na(x16_sar),
      !is.na(Pred1),
      !is.na(Pred2)
    )
  
  if (nrow(df_model) < 24) return(NULL)
  
  sem_text <- '
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
    x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + Pred1 + Pred2
  '
  
  fit <- tryCatch({
    sem(sem_text, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  
  fm <- fitMeasures(fit, c("aic", "bic", "cfi", "rmsea"))
  pe <- parameterEstimates(fit)
  
  get_est <- function(lhs_v, rhs_v, field) {
    val <- pe %>% filter(lhs == lhs_v, op == "~", rhs == rhs_v) %>% pull(!!sym(field))
    if (length(val) == 0) NA_real_ else val[1]
  }
  
  tibble(
    term1_col     = pred1_col,
    term2_col     = pred2_col,
    n             = nrow(df_model),
    aic           = fm[["aic"]],
    bic           = fm[["bic"]],
    cfi           = fm[["cfi"]],
    rmsea         = fm[["rmsea"]],
    
    b_hake_cpue   = get_est("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "est"),
    p_hake_cpue   = get_est("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "pvalue"),
    
    b_cpue_sar    = get_est("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "est"),
    p_cpue_sar    = get_est("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "pvalue"),
    
    b_term1_sar   = get_est("x16_sar", "Pred1", "est"),
    p_term1_sar   = get_est("x16_sar", "Pred1", "pvalue"),
    
    b_term2_sar   = get_est("x16_sar", "Pred2", "est"),
    p_term2_sar   = get_est("x16_sar", "Pred2", "pvalue")
  )
}


# -----------------------------------------------------------------------------
# 4A. SWEEP 1: Two Predator Combinations
# -----------------------------------------------------------------------------
pred_pairs <- combn(ak_pred_cands, 2, simplify = FALSE)

res_pred_pairs <- purrr::map_dfr(pred_pairs, function(pair) {
  fit_sem_two_topdown(pair[1], pair[2], guild.dfas1)
}) %>%
  mutate(model_type = "2 Predators")

# -----------------------------------------------------------------------------
# 4B. SWEEP 2: One Predator + One Prey Combinations
# -----------------------------------------------------------------------------
prey_pred_grid <- expand.grid(
  pred_col = ak_pred_cands,
  prey_col = ak_prey_cands,
  stringsAsFactors = FALSE
)

res_pred_prey <- purrr::map2_dfr(
  prey_pred_grid$pred_col, 
  prey_pred_grid$prey_col, 
  function(p1, p2) {
    fit_sem_two_topdown(p1, p2, guild.dfas1)
  }
) %>%
  mutate(model_type = "1 Predator + 1 Prey")

# -----------------------------------------------------------------------------
# 5. Combine, Compare against Baseline, and Rank
# -----------------------------------------------------------------------------
df_base <- guild.dfas1 %>%
  filter(year >= 1998, year <= 2021) %>%
  filter(!is.na(x07_dfa_cpue_int_spr_jun_hw), !is.na(x09_dfa_hake_age5plus), !is.na(x16_sar))

fit_base <- sem('
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw
', data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE)

base_aic <- fitMeasures(fit_base, "aic")

all_competing_models <- bind_rows(res_pred_pairs, res_pred_prey) %>%
  mutate(
    delta_aic_vs_base    = aic - base_aic,
    contains_new_shark   = (term1_col == "x15_ishark_z" | term2_col == "x15_ishark_z"),
    contains_new_ssl     = (term1_col == "x15_issl_z"   | term2_col == "x15_issl_z"),
    contains_new_indices = contains_new_shark | contains_new_ssl
  ) %>%
  arrange(aic)

# -----------------------------------------------------------------------------
# 6. Display Model Summaries
# -----------------------------------------------------------------------------

cat("\n=== TOP 15 OVERALL SEM MODELS ===\n")
#all_competing_models<-read.csv(file.path("copilot/outputs_5", "sharkDFA_ssl_sem_pred_prey_sst_x16_linear.csv"))
write_csv(all_competing_models, file.path(out_dir, "sharkDFA_ssl_sem_pred_prey_sst_x16_linear.csv"))

top_eco_model<-  all_competing_models %>% 
  filter(!grepl("x21_",term1_col),!grepl("x21_",term2_col)) %>%
  select(model_type, term1_col, term2_col, aic, delta_aic_vs_base, 
         b_term1_sar, p_term1_sar, b_term2_sar, p_term2_sar, contains_new_indices) %>% 
  head(15)
top_eco_model[,1:3]
#             model_type                          term1_col                                term2_col
# 1          2 Predators x15_dfa_sleeper_shark_bsai_pred_ak                               x15_issl_z
# 2          2 Predators                         x15_issl_z                  x15_ishark_DFAsleeper_z
# 3  1 Predator + 1 Prey                         x15_issl_z                            x13_ammod_wai
# 4  1 Predator + 1 Prey            x15_ishark_DFAsleeper_z           x12_dfa_biomass_euph_shelf_sum
# 5  1 Predator + 1 Prey            x15_ishark_DFAsleeper_z x13_pollock_biomass_go_aage3plus_pred_ak
