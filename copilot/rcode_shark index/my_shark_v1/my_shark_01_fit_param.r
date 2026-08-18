suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})

# ---------------------------
# Setup & Directory
# ---------------------------
out_dir <- "copilot/outputs_4"

if (!exists("out_dir")) out_dir <- "copilot/outputs_4"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Data Ingestion
# ---------------------------
sem_master_data <- read.csv("data_Lisa/sem_master_data.csv", row.names = NULL) %>% clean_names()
ak_dat          <- read.csv("data_Lisa/data_all_tested_columns_annual.csv")

# Join DFA Hake variable into main dataset
ak_dat <- left_join(
  sem_master_data %>% select(year, x09_dfa_hake_age5plus),
  ak_dat, 
  by = "year"
)

# ---------------------------
# Variable Definitions
# ---------------------------
col_salmon_juv   <- "x07_dfa_cpue_int_spr_jun_hw"
col_salmon_adult <- "x16_sar"

# Helper for standardizing candidate search
find_col <- function(df, patterns) {
  names(df)[grep(paste(patterns, collapse = "|"), names(df), ignore.case = TRUE)]
}

# Shark temperature candidates
shark_temp_cands <- c(
  "sst_wgoa_coastwatch_junjulaug",
  "bts_wgoa_temp_sum_195to205m",
  "swln_temp_spr_176to226m",
  "ll_wgoa_temp_sum_246to255m",
  "pdo_djf",
  "enso_dj"
)
shark_temp_cands <- shark_temp_cands[shark_temp_cands %in% names(ak_dat)]

if (length(shark_temp_cands) == 0) {
  stop("None of the chosen shark temperature candidates were found in ak_dat.")
}

# Shark abundance candidates
shark_cands <- find_col(ak_dat, c("sleeper_shark", "salmon_shark"))
if (length(shark_cands) == 0) {
  shark_cands <- find_col(sem_master_data, c("sleeper_shark", "salmon_shark"))
}
if (length(shark_cands) == 0) {
  stop("No shark abundance candidates found.")
}

# Safe scaling utility function
safe_scale <- function(x) {
  if (all(is.na(x))) return(x)
  as.vector(scale(x))
}

# ---------------------------
# Build Base Table (Cleaned of SSL & Forage)
# ---------------------------
data_base <- tibble(year = sort(unique(c(ak_dat$year, sem_master_data$year)))) %>%
  left_join(
    ak_dat %>% select(year, any_of(shark_temp_cands), any_of(shark_cands)),
    by = "year"
  ) %>%
  left_join(
    sem_master_data %>% select(
      year, 
      x09_dfa_hake_age5plus, 
      all_of(col_salmon_juv), 
      all_of(col_salmon_adult)
    ),
    by = "year"
  )

# ---------------------------
# SEM Model Fitting Function
# ---------------------------
fit_shark_variant_sem <- function(col_shark, shark_temp_col, Q10, overlap_form, overlap_slope, shark_transform) {
  
  # 1. Shark Abundance Transformation
  if (!col_shark %in% names(data_base)) return(NULL)
  shark_vec <- data_base[[col_shark]]
  
  if (shark_transform == "z") {
    Shark_use <- safe_scale(shark_vec)
  } else if (shark_transform == "log1p_z") {
    minv <- suppressWarnings(min(shark_vec, na.rm = TRUE))
    if (!is.finite(minv)) return(NULL)
    Shark_use <- safe_scale(log1p(shark_vec - minv))
  } else if (shark_transform == "z_roll2") {
    z0   <- safe_scale(shark_vec)
    zlag <- dplyr::lag(z0, 1)
    Shark_use <- rowMeans(cbind(z0, zlag), na.rm = TRUE)
  } else {
    return(NULL)
  }
  
  # 2. Shark Temperature Driver & Metabolic Response
  if (!shark_temp_col %in% names(data_base)) return(NULL)
  T_shark_raw <- data_base[[shark_temp_col]]
  T_shark     <- safe_scale(T_shark_raw)
  
  Tref_shark  <- mean(T_shark_raw, na.rm = TRUE)
  M_t         <- Q10^((T_shark_raw - Tref_shark) / 10)
  
  # 3. Spatial Overlap Form
  O_t <- if (overlap_form == "constant") {
    rep(1, nrow(data_base))
  } else if (overlap_form == "linear") {
    pmax(0.1, 1 + overlap_slope * T_shark)
  } else if (overlap_form == "logistic") {
    plogis(overlap_slope * T_shark) * 2
  } else {
    return(NULL)
  }
  
  # 4. Construct Integrated Index & Prep Model Frame
  fit_dat <- data_base %>%
    mutate(I_Shark = Shark_use * M_t * O_t) %>%
    filter(year >= 1998, year <= 2021) %>%
    rename(
      x07_dfa_cpue_int_spr_jun_hw = !!sym(col_salmon_juv),
      x16_sar                  = !!sym(col_salmon_adult)
    ) %>%
    filter(
      !is.na(x16_sar),
      !is.na(x07_dfa_cpue_int_spr_jun_hw),
      !is.na(x09_dfa_hake_age5plus),
      !is.na(I_Shark)
    )
  
  if (nrow(fit_dat) < 12) return(NULL)
  
  # 5. Lavaan SEM Model Syntax
  single_model_text <- '
    # Structural Model
    x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus 
    x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + I_Shark
  '
  
  # 6. Fit SEM with Error Handling
  fit_sem <- tryCatch({
    lavaan::sem(single_model_text, data = fit_dat, missing = "ML", warn = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(fit_sem) || !lavaan::lavInspect(fit_sem, "converged")) return(NULL)
  
  # 7. Extract Model Metrics and Parameters
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
    
    # Structural Path Coefficients & P-Values
    b_hake_cpue      = get_path("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "est"),
    p_hake_cpue      = get_path("x07_dfa_cpue_int_spr_jun_hw", "x09_dfa_hake_age5plus", "pvalue"),
    b_cpue_sar       = get_path("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "est"),
    p_cpue_sar       = get_path("x16_sar", "x07_dfa_cpue_int_spr_jun_hw", "pvalue"),
    b_shark_sar      = get_path("x16_sar", "I_Shark", "est"),
    p_shark_sar      = get_path("x16_sar", "I_Shark", "pvalue")
  )
}

# ---------------------------
# Parameter Grid Sweep
# ---------------------------
grid <- tidyr::crossing(
  shark_col       = shark_cands,
  shark_temp_col  = shark_temp_cands,
  Q10             = c(1.5, 2.0, 2.5, 3.0),
  overlap_form    = c("constant", "linear", "logistic"),
  overlap_slope   = c(0.2, 0.4, 0.6, 0.8),
  shark_transform = c("z", "log1p_z", "z_roll2")
) %>%
  mutate(overlap_slope = if_else(overlap_form == "constant", 0, overlap_slope)) %>%
  distinct()

# Run iterations
res_shark <- purrr::pmap_dfr(
  grid,
  function(shark_col, shark_temp_col, Q10, overlap_form, overlap_slope, shark_transform) {
    fit_shark_variant_sem(
      col_shark       = shark_col,
      shark_temp_col  = shark_temp_col,
      Q10             = Q10,
      overlap_form    = overlap_form,
      overlap_slope   = overlap_slope,
      shark_transform = shark_transform
    )
  }
)

if (nrow(res_shark) == 0) stop("No shark SEM variants converged.")

# ---------------------------
# Output & Ranking
# ---------------------------
res_shark <- res_shark %>%
  mutate(
    model_group = case_when(n == 19 ~ "short", n == 24 ~ "long", TRUE ~ NA_character_),
    sig_shark   = !is.na(p_shark_sar) & p_shark_sar < 0.05
  )

write_csv(res_shark, file.path(out_dir, "shark_sem_sweep_all.csv"))

# Rank models by structural path significance and lowest AIC
ranked <- res_shark %>%
  mutate(
    sig_hake  = !is.na(p_hake_cpue) & p_hake_cpue < 0.05,
    sig_cpue  = !is.na(p_cpue_sar)  & p_cpue_sar  < 0.05,
    sig_shark = !is.na(p_shark_sar) & p_shark_sar < 0.05
  ) %>%
  arrange(
    desc(sig_shark),     # Prioritize models where the shark effect is significant
    desc(sig_cpue),      # Next prioritize juvenile salmon CPUE -> SAR path
    aic,                 # Lowest AIC
    p_shark_sar,         # Strongest p-value for shark term
    desc(abs(b_shark_sar))
  )

head(ranked)

write_csv(ranked, file.path(out_dir, "shark_sem_sweep_ranked.csv"))

# ---------------------------
# Parameter Importance via AIC Weights
# ---------------------------
w <- res_shark %>%
  filter(!is.na(model_group), !is.na(aic)) %>%
  group_by(model_group) %>%
  mutate(
    delta_aic  = aic - min(aic, na.rm = TRUE),
    rel_like   = exp(-0.5 * delta_aic),
    aic_weight = rel_like / sum(rel_like, na.rm = TRUE)
  ) %>%
  ungroup()

imp <- bind_rows(
  w %>% group_by(model_group, shark_col)       %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_col", level=shark_col),
  w %>% group_by(model_group, shark_temp_col)  %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_temp_col", level=shark_temp_col),
  w %>% group_by(model_group, Q10)             %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="Q10", level=as.character(Q10)),
  w %>% group_by(model_group, overlap_form)    %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="overlap_form", level=overlap_form),
  w %>% group_by(model_group, overlap_slope)   %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="overlap_slope", level=as.character(overlap_slope)),
  w %>% group_by(model_group, shark_transform) %>% summarise(importance_aicw = sum(aic_weight), .groups="drop") %>% mutate(param="shark_transform", level=shark_transform)
) %>%
  group_by(model_group, param) %>%
  mutate(
    rank = rank(-importance_aicw, ties.method="min"),
    importance_rel_top1 = importance_aicw / max(importance_aicw)
  ) %>%
  ungroup()

imp

write_csv(imp, file.path(out_dir, "shark_param_importance_sem.csv"))

message("Finished! SEM sweep outputs saved in: ", out_dir)
