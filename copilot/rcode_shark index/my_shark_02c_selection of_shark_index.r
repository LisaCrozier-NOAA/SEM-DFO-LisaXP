

#RESULTS
# shark_col                 importance_aicw param           level                         shark_temp_col                overlap_form overlap_slope shark_transform  rank importance_rel_top1
# 1 goa_pacific_sleeper_shark         0.948   shark_col       goa_pacific_sleeper_shark     NA                            NA                    NA   NA                  1             1      
# 2 NA                                0.905   shark_transform log1p                         NA                            NA                    NA   log1p               1             1      
# 3 NA                                0.811   shark_temp_col  sst_wgoa_coastwatch_junjulaug sst_wgoa_coastwatch_junjulaug NA                    NA   NA                  1             1      
# 4 NA                                0.574   overlap_form    linear                        NA                            linear                NA   NA                  1             1      
# 5 NA                                0.424   overlap_form    logistic                      NA                            logistic              NA   NA                  2             0.739  
# 6 NA                                0.308   overlap_slope   1                             NA                            NA                     1   NA                  1             1      
# 7 NA                                0.283   overlap_slope   0.8                           NA                            NA                     0.8 NA                  2             0.918  
# 8 NA                                0.227   overlap_slope   0.6                           NA                            NA                     0.6 NA                  3             0.736  
# 9 NA                                0.181   shark_temp_col  sst_egoa_coastwatch_junjulaug sst_egoa_coastwatch_junjulaug NA                    NA   NA                  2             0.223  
# 10 NA                                0.144   overlap_slope   0.4                           NA                            NA                     0.4 NA                  4             0.466  
# 11 NA                                0.0951  shark_transform raw                           NA                            NA                    NA   raw                 2             0.105  
# 12 goa_salmon_shark                  0.0520  shark_col       goa_salmon_shark              NA                            NA                    NA   NA                  2             0.0548 
# 13 NA                                0.0362  overlap_slope   0.2                           NA                            NA                     0.2 NA                  5             0.117  
# 14 NA                                0.00801 shark_temp_col  swln_temp_fall_176to226m      swln_temp_fall_176to226m      NA                    NA   NA                  3             0.00987
# 15 NA                                0.00240 overlap_form    constant                      NA                            constant              NA   NA                  3             0.00418
# 16 NA                                0.00240 overlap_slope   0                             NA                            NA                     0   NA                  6             0.00777
>


# shark_col              shark_transform shark_temp_col overlap_form   aic aic_base delta_aic_vs_base improved_fit p_shark_sar p_shark_base
# 1 goa_pacific_sleeper_s… log1p           sst_wgoa_coas… linear        92.9     106.            -12.8  TRUE          0.00000549        0.142
# 2 goa_pacific_sleeper_s… log1p           sst_wgoa_coas… logistic      93.0     106.            -12.7  TRUE          0.00000643        0.142
# 3 goa_pacific_sleeper_s… log1p           sst_wgoa_coas… linear        93.1     106.            -12.6  TRUE          0.00000665        0.142
# 4 goa_pacific_sleeper_s… log1p           sst_wgoa_coas… linear        93.3     106.            -12.4  TRUE          0.00000855        0.142
# 5 goa_pacific_sleeper_s… log1p           sst_wgoa_coas… logistic      94.3     106.            -11.4  TRUE          0.0000218         0.142
# 6 goa_pacific_sleeper_s… log1p           sst_egoa_coas… linear        96.2     106.             -9.51 TRUE          0.000115          0.142
# 7 goa_pacific_sleeper_s… log1p           sst_egoa_coas… linear        96.2     106.             -9.47 TRUE          0.000119          0.142
# 8 goa_pacific_sleeper_s… log1p           sst_wgoa_coas… linear        96.4     106.             -9.27 TRUE          0.000141          0.142
# 9 goa_pacific_sleeper_s… log1p           sst_wgoa_coas… logistic      96.5     106.             -9.22 TRUE          0.000147          0.142
# 10 goa_pacific_sleeper_s… log1p           sst_egoa_coas… linear        96.9     106.             -8.83 TRUE          0.000206          0.142



suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(janitor)
})

# ---------------------------
# Setup & Directory
# ---------------------------
out_dir <- "copilot/outputs_shark"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Assume sem.dat is loaded in global env
sem.dat<-read.csv(file.path("data_Lisa/sem_master_data.csv"))%>% 
  clean_names()

salmon_dat <- sem.dat %>% 
  select(year, contains("x07"), contains("x16_sar"), contains("sar_"), contains("x09_dfa_hake"))

shark_dat <- read.csv(file.path("data_Lisa/shark_wide.csv")) %>% 
  clean_names() %>%
  select(year, goa_pacific_sleeper_shark, goa_salmon_shark)

sst_dat <- read.csv(file.path("data_Lisa/goa_prey_clim_raw_trends_avg.csv"), row.names = 1) %>% 
  clean_names() %>%
  select(year, sst_egoa_coastwatch_junjulaug, sst_wgoa_coastwatch_junjulaug, swln_temp_fall_176to226m) %>%
  mutate(
    dT_sst_egoa = sst_egoa_coastwatch_junjulaug - mean(sst_egoa_coastwatch_junjulaug, na.rm = TRUE),
    dT_sst_wgoa = sst_wgoa_coastwatch_junjulaug - mean(sst_wgoa_coastwatch_junjulaug, na.rm = TRUE),
    dT_SewardLine_176to226m = swln_temp_fall_176to226m - mean(swln_temp_fall_176to226m, na.rm = TRUE),
    sc_sst_egoa = as.vector(scale(sst_egoa_coastwatch_junjulaug)),
    sc_sst_wgoa = as.vector(scale(sst_wgoa_coastwatch_junjulaug)),
    sc_SewardLine_176to226m = as.vector(scale(swln_temp_fall_176to226m))        
  )

col_salmon_juv   <- "x07_dfa_cpue_int_spr_jun_hw"
col_salmon_adult <- "x16_sar"

shark_temp_cands_raw <- c("sst_egoa_coastwatch_junjulaug", "sst_wgoa_coastwatch_junjulaug", "swln_temp_fall_176to226m")
shark_cands          <- names(shark_dat)[-1]

safe_scale <- function(x) {
  if (all(is.na(x))) return(x)
  as.vector(scale(x))
}

data_base <- salmon_dat %>%
  left_join(shark_dat, by = "year") %>%
  left_join(sst_dat, by = "year")

# ---------------------------
# 1. Integrated Index Function
# ---------------------------
fit_shark_variant_sem <- function(col_shark, shark_temp_col, Q10 = 2, overlap_form, overlap_slope, shark_transform) {
  
  if (!col_shark %in% names(data_base)) return(NULL)
  shark_vec <- data_base[[col_shark]]
  
  # Abundance Transformation (Strictly Non-Negative)
  if (shark_transform == "raw") {
    Shark_pos <- shark_vec
  } else if (shark_transform == "log1p") {
    minv <- suppressWarnings(min(shark_vec, na.rm = TRUE))
    if (!is.finite(minv)) return(NULL)
    Shark_pos <- log1p(shark_vec - minv)
  } else {
    return(NULL)
  }
  
  # Temperature & Metabolic Response
  if (!shark_temp_col %in% names(data_base)) return(NULL)
  T_shark_raw <- data_base[[shark_temp_col]]
  T_shark_z   <- safe_scale(T_shark_raw)
  
  Tref_shark  <- mean(T_shark_raw, na.rm = TRUE)
  M_t         <- Q10^((T_shark_raw - Tref_shark) / 10)
  
  # Spatial Overlap Form
  O_t <- if (overlap_form == "constant") {
    rep(1, nrow(data_base))
  } else if (overlap_form == "linear") {
    pmax(0.1, 1 + overlap_slope * T_shark_z)
  } else if (overlap_form == "logistic") {
    plogis(overlap_slope * T_shark_z) * 2
  } else {
    return(NULL)
  }
  
  # Integrated Index in Positive Space, THEN Z-scale
  fit_dat <- data_base %>%
    mutate(
      I_Shark_raw = Shark_pos * M_t * O_t,
      I_Shark     = safe_scale(I_Shark_raw)
    ) %>%
    filter(year >= 1998, year <= 2021) %>%
    rename(
      x07_dfa_cpue_int_spr_jun_hw = !!sym(col_salmon_juv),
      x16_sar                   = !!sym(col_salmon_adult)
    ) %>%
    filter(
      !is.na(x16_sar),
      !is.na(x07_dfa_cpue_int_spr_jun_hw),
      !is.na(x09_dfa_hake_age5plus),
      !is.na(I_Shark)
    )
  
  if (nrow(fit_dat) < 12) return(NULL)
  
  single_model_text <- '
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus 
    x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + I_Shark
  '
  
  fit_sem <- tryCatch({
    lavaan::sem(single_model_text, data = fit_dat, missing = "ML", warn = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(fit_sem) || !lavaan::lavInspect(fit_sem, "converged")) return(NULL)
  
  fit_measures <- lavaan::fitMeasures(fit_sem, c("aic", "bic", "cfi", "rmsea"))
  pe           <- lavaan::parameterEstimates(fit_sem)
  
  get_path <- function(lhs_val, rhs_val, field) {
    val <- pe %>% 
      filter(lhs == lhs_val, op == "~", rhs == rhs_val) %>% 
      pull(!!sym(field))
    if (length(val) == 0) NA_real_ else val[1]
  }
  
  tibble(
    shark_col        = col_shark,
    shark_temp_col   = shark_temp_col,
    Q10              = Q10,
    overlap_form     = overlap_form,
    overlap_slope    = overlap_slope,
    shark_transform  = shark_transform,
    n                = nrow(fit_dat),
    
    aic              = fit_measures[["aic"]],
    bic              = fit_measures[["bic"]],
    cfi              = fit_measures[["cfi"]],
    rmsea            = fit_measures[["rmsea"]],
    
    b_hake_cpue      = get_path("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "est"),
    p_hake_cpue      = get_path("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "pvalue"),
    b_cpue_sar       = get_path("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "est"),
    p_cpue_sar       = get_path("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "pvalue"),
    b_shark_sar      = get_path("x16_sar", "I_Shark", "est"),
    p_shark_sar      = get_path("x16_sar", "I_Shark", "pvalue")
  )
}

# ---------------------------
# 2. Baseline Model Function (Unadjusted Shark Input)
# ---------------------------
fit_baseline_sem <- function(col_shark, shark_transform) {
  if (!col_shark %in% names(data_base)) return(NULL)
  shark_vec <- data_base[[col_shark]]
  
  if (shark_transform == "raw") {
    Shark_pos <- shark_vec
  } else if (shark_transform == "log1p") {
    minv <- suppressWarnings(min(shark_vec, na.rm = TRUE))
    if (!is.finite(minv)) return(NULL)
    Shark_pos <- log1p(shark_vec - minv)
  } else {
    return(NULL)
  }
  
  fit_dat <- data_base %>%
    mutate(I_Shark = safe_scale(Shark_pos)) %>% # Simply scale raw/log count
    filter(year >= 1998, year <= 2021) %>%
    rename(
      x07_dfa_cpue_int_spr_jun_hw = !!sym(col_salmon_juv),
      x16_sar                   = !!sym(col_salmon_adult)
    ) %>%
    filter(
      !is.na(x16_sar),
      !is.na(x07_dfa_cpue_int_spr_jun_hw),
      !is.na(x09_dfa_hake_age5plus),
      !is.na(I_Shark)
    )
  
  if (nrow(fit_dat) < 12) return(NULL)
  
  single_model_text <- '
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus 
    x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + I_Shark
  '
  
  fit_sem <- tryCatch({
    lavaan::sem(single_model_text, data = fit_dat, missing = "ML", warn = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(fit_sem) || !lavaan::lavInspect(fit_sem, "converged")) return(NULL)
  
  fit_measures <- lavaan::fitMeasures(fit_sem, c("aic", "bic"))
  pe           <- lavaan::parameterEstimates(fit_sem)
  
  get_path <- function(lhs_val, rhs_val, field) {
    val <- pe %>% 
      filter(lhs == lhs_val, op == "~", rhs == rhs_val) %>% 
      pull(!!sym(field))
    if (length(val) == 0) NA_real_ else val[1]
  }
  
  tibble(
    shark_col       = col_shark,
    shark_transform = shark_transform,
    n               = nrow(fit_dat),
    aic_base        = fit_measures[["aic"]],
    bic_base        = fit_measures[["bic"]],
    b_shark_base    = get_path("x16_sar", "I_Shark", "est"),
    p_shark_base    = get_path("x16_sar", "I_Shark", "pvalue")
  )
}

# ---------------------------
# 3. Grid Sweep (Raw and Log1p ONLY)
# ---------------------------
grid_integrated <- tidyr::crossing(
  shark_col       = shark_cands,
  shark_temp_col  = shark_temp_cands_raw,
  Q10             = 2,
  overlap_form    = c("constant", "linear", "logistic"),
  overlap_slope   = c(0.2, 0.4, 0.6, 0.8,1),
  shark_transform = c("raw", "log1p") # EXCLUSIVELY raw and log1p
) %>%
  mutate(overlap_slope = if_else(overlap_form == "constant", 0, overlap_slope)) %>%
  distinct()

# Run Integrated Sweep
res_shark <- purrr::pmap_dfr(
  grid_integrated,
  function(shark_col, shark_temp_col, Q10, overlap_form, overlap_slope, shark_transform) {
    fit_shark_variant_sem(shark_col, shark_temp_col, Q10, overlap_form, overlap_slope, shark_transform)
  }
)

# Run Baseline Sweep
grid_baseline <- tidyr::crossing(
  shark_col       = shark_cands,
  shark_transform = c("raw", "log1p")
)

res_baseline <- purrr::pmap_dfr(
  grid_baseline,
  function(shark_col, shark_transform) {
    fit_baseline_sem(shark_col, shark_transform)
  }
)

# ---------------------------
# 4. Compare Integrated vs. Raw Baseline
# ---------------------------
comparison_ranked <- res_shark %>%
  left_join(res_baseline, by = c("shark_col", "shark_transform", "n")) %>%
  mutate(
    delta_aic_vs_base = aic - aic_base, # Negative = Integrated model is BETTER than raw baseline
    improved_fit      = delta_aic_vs_base < -2,
    sig_shark_int     = !is.na(p_shark_sar) & p_shark_sar < 0.05,
    sig_shark_base    = !is.na(p_shark_base) & p_shark_base < 0.05
  ) %>%
  arrange(delta_aic_vs_base)

# Print top comparison results to console
print(comparison_ranked %>% 
        select(
          shark_col, shark_transform, shark_temp_col, overlap_form, 
          aic, aic_base, delta_aic_vs_base, improved_fit, 
          p_shark_sar, p_shark_base
        ) %>% head(10))


comparison_ranked %>% 
  filter(shark_transform == "raw") %>% 
  select(
    shark_col, shark_transform, shark_temp_col, overlap_form, 
    aic, aic_base, delta_aic_vs_base, improved_fit, 
    p_shark_sar, p_shark_base
  ) %>% 
  head(10)

comparison_ranked %>% 
  group_by(shark_col, shark_transform) %>% 
  slice_min(order_by = delta_aic_vs_base, n = 1) %>% 
  ungroup() %>% 
  select(
    shark_col, shark_transform, shark_temp_col, overlap_form, 
    aic, aic_base, delta_aic_vs_base, p_shark_sar, p_shark_base
  )

# Save comparison results
write_csv(comparison_ranked, file.path(out_dir, "shark_index_vs_raw_baseline_comparison.csv"))



#IMPORTANCE------------
# ---------------------------
# Parameter Importance via AIC Weights
# ---------------------------
w <- res_shark %>%
  filter( !is.na(aic)) %>%
  mutate(
    delta_aic  = aic - min(aic, na.rm = TRUE),
    rel_like   = exp(-0.5 * delta_aic),
    aic_weight = rel_like / sum(rel_like, na.rm = TRUE)
  ) 

imp <- bind_rows(
  w %>% group_by( shark_col)       %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_col", level=shark_col),
  w %>% group_by( shark_temp_col)  %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_temp_col", level=shark_temp_col),
  w %>% group_by( overlap_form)    %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="overlap_form", level=overlap_form),
  w %>% group_by( overlap_slope)   %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="overlap_slope", level=as.character(overlap_slope)),
  w %>% group_by( shark_transform) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_transform", level=shark_transform)
) %>%
  group_by( param) %>%
  mutate(
    rank = rank(-importance_aicw, ties.method="min"),
    importance_rel_top1 = importance_aicw / max(importance_aicw)
  ) %>%
  ungroup()%>% 
  arrange(desc(importance_aicw))

(imp)

imp %>% 
  filter(param == "overlap_form") %>% 
  arrange(desc(importance_aicw))

imp %>% 
  filter(param == "overlap_slope") %>% 
  arrange(desc(importance_aicw))

write_csv(imp, file.path(out_dir, "shark_param_importance_sem_x16.csv"))


#PLOT INDEX ----------
library(tidyverse)
library(ggplot2)

# Set global theme for publication-ready plots
theme_set(theme_bw(base_size = 12))

# ---------------------------
# 1. Prepare Plot Data matching the SEM frame
# ---------------------------
top_shark_col <- "goa_pacific_sleeper_shark"
top_temp_col  <- "sst_wgoa_coastwatch_junjulaug"
top_slope     <- 0.8  # Top model slope

plot_dat <- data_base %>%
  # Filter to exact SEM model years and rename salmon variables
  filter(year >= 1998, year <= 2021) %>%
  rename(
    x07_dfa_cpue_int_spr_jun_hw = !!sym(col_salmon_juv),
    x16_sar                   = !!sym(col_salmon_adult)
  ) %>%
  filter(
    !is.na(x16_sar),
    !is.na(x07_dfa_cpue_int_spr_jun_hw),
    !is.na(x09_dfa_hake_age5plus),
    !is.na(!!sym(top_shark_col)),
    !is.na(!!sym(top_temp_col))
  ) %>%
  mutate(
    # A) Raw log1p count (Z-scaled)
    N_log1p = safe_scale(log1p(!!sym(top_shark_col))),
    
    # B) Integrated Index (WGOA SST + Logistic Overlap)
    T_raw   = !!sym(top_temp_col),
    T_z     = safe_scale(T_raw),
    M_t     = 2^((T_raw - mean(T_raw, na.rm = TRUE)) / 10),
    O_t     = plogis(top_slope * T_z) * 2,
    
    I_Shark_raw = log1p(!!sym(top_shark_col)) * M_t * O_t,
    I_Shark_z   = safe_scale(I_Shark_raw)
  )

# ---------------------------
# PLOT 1: Time Series (Raw log1p vs Integrated Index)
# ---------------------------
p1_dat <- plot_dat %>%
  select(year, `Log1p Raw Shark` = N_log1p, `Climate-Integrated Index` = I_Shark_z) %>%
  pivot_longer(-year, names_to = "Metric", values_to = "Z_Score")

p1 <- ggplot(p1_dat, aes(x = year, y = Z_Score, color = Metric, linetype = Metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Log1p Raw Shark" = "gray40", "Climate-Integrated Index" = "firebrick")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
  labs(
    title = "Sleeper Shark: Raw Counts vs. Climate-Integrated Index",
    subtitle = "Heatwave years (e.g., 2014-2016) amplify the integrated index relative to raw counts",
    x = "Year",
    y = "Standardized Metric (Z-score)",
    color = "Index Type",
    linetype = "Index Type"
  ) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(p1)
ggsave(file.path(out_dir, "plot1_shark_index_timeseries.png"), plot = p1, width = 9, height = 5)

# ---------------------------
# PLOT 2: No-Shark Target Residuals vs. Shark Predation
# ---------------------------
# Fit the base SEM path without sharks (SAR ~ Juv_CPUE)
no_shark_lm <- lm(x16_sar ~ x07_dfa_cpue_int_spr_jun_hw, data = plot_dat)

res_dat <- plot_dat %>%
  mutate(
    unexplained_sar_residual = residuals(no_shark_lm),
    I_Shark_z = I_Shark_z
  )

p2 <- ggplot(res_dat, aes(x = year)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Blue bars: Unexplained SAR Residuals
  geom_col(aes(y = unexplained_sar_residual), fill = "steelblue", alpha = 0.7, width = 0.6) +
  # Red line: Inverted Integrated Shark Index (-I_Shark)
  geom_line(aes(y = -I_Shark_z), color = "firebrick", linewidth = 1.1) +
  geom_point(aes(y = -I_Shark_z), color = "firebrick", size = 2) +
  scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
  labs(
    title = "Target Residuals vs. Integrated Shark Predation",
    subtitle = "Blue Bars = Unexplained SAR Residuals (No-Shark Model) | Red Line = Inverted Shark Index (-I_Shark)",
    x = "Year",
    y = "Unexplained SAR Anomaly / Inverted Shark Index"
  ) +
  theme(panel.grid.minor = element_blank())

print(p2)
ggsave(file.path(out_dir, "plot2_no_shark_target_residuals.png"), plot = p2, width = 9, height = 5)
