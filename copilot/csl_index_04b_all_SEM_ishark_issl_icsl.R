

#Load packages---------
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})



out_dir <- "copilot/outputs_7"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#1) setup the data ----

guild.dfasAK <- read.csv(file.path("data_Lisa/guild.dfasAK.csv")) %>% clean_names()

ssl.dat <- read.csv("data_Lisa/ssl.dat.csv") %>% clean_names() %>%
  select(year,ssl_seak_pup_pred,sst_wgoa_coastwatch_junjulaug) %>%
  left_join(guild.dfasAK %>%
              select(year, x13_stka_herr_matbiom,x13_egoa_herring,x13_mid_il_capelin),
            join_by("year"))

#Shark function------
ishark_fxn <- function(N_shark, 
                       T_mt,
                       Tref_mt=mean(T_mt, na.rm = TRUE), 
                       T_ot, 
                       Q10 = 2, 
                       overlap_slope = 1, 
                       transform_shark = "log1p", 
                       scale_final = TRUE) {
  
  # 1. Input checks
  if (length(N_shark) != length(T_mt) || length(N_shark) != length(T_ot)) {
    stop("Input vectors (N_shark, T_mt, T_ot) must be of equal length.")
  }
  
  # 2. Transform Abundance (Keep non-negative)
  if (transform_shark == "log1p") {
    minv <- suppressWarnings(min(N_shark, na.rm = TRUE))
    N_pos <- log1p(N_shark - minv)
  } else if (transform_shark == "raw") {
    N_pos <- N_shark
  } else {
    stop("transform_shark must be either 'log1p' or 'raw'.")
  }
  
  # 3. Calculate Metabolic Multiplier (Mt) using raw T_mt differential
  M_t     <- Q10^((T_mt - Tref_mt) / 10)
  
  # 4. Calculate Spatial Overlap Multiplier (Ot) using Z-scaled T_ot for Logistic curve
  # Safe Z-scaling function
  scale_vector <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }
  
  T_ot_z <- scale_vector(T_ot)
  O_t    <- plogis(overlap_slope * T_ot_z) * 2
  
  # 5. Compute Integrated Predation Index in Positive Space
  I_Shark_raw <- N_pos * M_t * O_t
  
  # 6. Standardize or Return Raw Index
  if (scale_final) {
    I_Shark_out <- scale_vector(I_Shark_raw)
  } else {
    I_Shark_out <- I_Shark_raw
  }
  
  return(I_Shark_out)
}


shark_raw<-read.csv(file.path("data_Lisa/shark_wide.csv")) %>% 
  clean_names() %>%
  select(year,goa_pacific_sleeper_shark) %>%
  full_join(guild.dfasAK %>%
              select(year,x21_sst_wgoa_coastwatch_junjulaug),
            by="year") %>%
  mutate(ishark=ishark_fxn(
    N_shark=goa_pacific_sleeper_shark,
    T_mt=x21_sst_wgoa_coastwatch_junjulaug,
    T_ot=x21_sst_wgoa_coastwatch_junjulaug))
head(shark_raw)

fit_dat<-guild.dfasAK %>% full_join(shark_raw %>% select(year,ishark),by="year") %>%
  select(year,x07_dfa_cpue_int_spr_jun_hw,x09_dfa_hake_age5plus,x16_sar,ishark)


df_top_ssl <- ssl.dat %>%
  mutate(
    # Core predictors
    ssl = as.vector(scale(ssl_seak_pup_pred)),
    sst = as.vector(scale(sst_wgoa_coastwatch_junjulaug)),
    f_cap = as.vector(scale(x13_mid_il_capelin)),
    f_herr = as.vector(scale(x13_stka_herr_matbiom)),
    
    # Linear interaction products
    ssl_x_sst  = ssl * sst,
    ssl_x_cap  = ssl * f_cap,
    ssl_x_herr = ssl * f_herr
  )%>%
  full_join(fit_dat,by="year")

#Define ssl index----------
df_top_ssl <- df_top_ssl %>%
  mutate(
    # Calculate Simplified Hazard Index (No empirical regression tuning!)
    # Negated so higher = higher predation hazard
    I_SSL_simple = -1 * ( 1.0 * ssl - 0.3 * (ssl * sst) + 2.0 * (ssl * f_herr) + 1.0 * (ssl * f_cap) ),
    
    # Standardize the resulting composite index
    issl = as.vector(scale(I_SSL_simple))
  )

guild.dfasAK_ishark_issl<-full_join(guild.dfasAK,
                                    df_top_ssl %>% select(year,ishark,issl),
                                    by="year")
write.csv(guild.dfasAK_ishark_issl,   file.path("copilot/outputs_6", "guild.dfasAK_ishark_issl.csv"), row.names = FALSE)


guild.dfasAK_ishark_issl<-read.csv(   file.path("copilot/outputs_6", "guild.dfasAK_ishark_issl.csv"), row.names = FALSE)





#Add csl index------

csl_path   <- file.path("copilot/outputs_6", "chinook_numeric_csl_exposure_risk_1998_2024.csv")

csl_raw <- read.csv(csl_path) %>% clean_names()

# ---------------------------
# 1) Build lag-2 CSL index (idxE10)
# For SAR year t, use CSL/chinook at t+2 adult return year
# ---------------------------
safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - m) / s)
}

csl_idx <- csl_raw %>%
  transmute(
    year,
    eul_z       = safe_z(eulachon_during_chinook),
    shad_z      = safe_z(shad_during_chinook),
    x15_icsl_risk_eulachon75_shad25 = -(0.75 * eul_z + 0.25 * shad_z),
    x15_icsl_cslchin_ratio = safe_z(log1p(csl_during_chinook) - log1p(chinook_estuary_total))
  ) %>%
  mutate(year = year - 2)  # align t+2 CSL to SAR year t


# ---------------------------
# 2) Assemble guild table as in original workflow------
# ---------------------------
guild_raw<-guild.dfasAK_ishark_issl

guild <- guild_raw %>%
  left_join(csl_idx, by = "year") %>%
  mutate(across(-year, ~ as.vector(scale(.))))

guild_dfas1_24yr <- guild %>%
  filter(year >= 1998, year <= 2021) %>%
  select(where(~ sum(!is.na(.)) >= 24))%>%
  rename(X15_ishark=ishark, X15_issl=issl) %>%
  select(-x15_shark_catch_go_a_pred_ak,-x10_harbor_seal_cr_2yr_lead)

sort(names(guild_dfas1_24yr))


# ---------------------------
# 3) Candidate sets (same logic + new CSL index)--------
# ---------------------------
ak_prey_cands <- names(guild_dfas1_24yr)[grepl("^x21_|^x12_|^x13_|^x14_", names(guild_dfas1_24yr))]
ak_pred_cands <- names(guild_dfas1_24yr)[grepl("^x10_|^x15_", names(guild_dfas1_24yr))]
all_cands <-c(ak_prey_cands,ak_pred_cands)

message("Prey candidates: ", length(ak_prey_cands))
message("Pred candidates: ", length(ak_pred_cands))

# response and structural base terms
base_needed <- c("year","x07_dfa_cpue_int_spr_jun_hw","x09_dfa_hake_age5plus","x16_sar")
if (!all(base_needed %in% names(guild_dfas1_24yr))) {
  stop("Missing required SEM columns in guild_dfas1_24yr.")
}

# ---------------------------
# 4) General SEM fitting helper for k extra terms
# ---------------------------
fit_sem_k <- function(extra_terms, df_data, model_type_label) {
  keep_cols <- unique(c(base_needed, extra_terms))
  df_model <- df_data %>%
    filter(year >= 1998, year <= 2021) %>%
    select(all_of(keep_cols)) %>%
    drop_na()
  
  if (nrow(df_model) < 24) return(NULL)
  
  rhs <- paste(c("x07_dfa_cpue_int_spr_jun_hw", extra_terms), collapse = " + ")
  sem_text <- paste0(
    "x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus\n",
    "x16_sar ~ ", rhs, "\n"
  )
  
  fit <- tryCatch(
    sem(sem_text, data = df_model, std.lv = TRUE, missing = "ML", warn = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  
  fm <- fitMeasures(fit, c("aic", "bic", "cfi", "rmsea"))
  r2_sar <- inspect(fit, "r2")[["x16_sar"]]
  pe <- parameterEstimates(fit)
  
  # Extract coeffs for added terms
  get_coef <- function(term) {
    row <- pe %>% filter(lhs == "x16_sar", op == "~", rhs == term)
    if (nrow(row) == 0) return(c(est = NA_real_, p = NA_real_))
    c(est = row$est[1], p = row$pvalue[1])
  }
  coef_list <- lapply(extra_terms, get_coef)
  
  tibble(
    model_type = model_type_label,
    terms = paste(extra_terms, collapse = " + "),
    n = nrow(df_model),
    aic = fm[["aic"]],
    bic = fm[["bic"]],
    cfi = fm[["cfi"]],
    rmsea = fm[["rmsea"]],
    sar_r2 = as.numeric(r2_sar),
    term1 = extra_terms[1],
    term2 = ifelse(length(extra_terms) >= 2, extra_terms[2], NA_character_),
    term3 = ifelse(length(extra_terms) >= 3, extra_terms[3], NA_character_),
    b1 = coef_list[[1]]["est"], p1 = coef_list[[1]]["p"],
    b2 = ifelse(length(extra_terms) >= 2, coef_list[[2]]["est"], NA_real_),
    p2 = ifelse(length(extra_terms) >= 2, coef_list[[2]]["p"], NA_real_),
    b3 = ifelse(length(extra_terms) >= 3, coef_list[[3]]["est"], NA_real_),
    p3 = ifelse(length(extra_terms) >= 3, coef_list[[3]]["p"], NA_real_)
  )
}

# ---------------------------
# 5) Baseline model for comparison
# ---------------------------
df_base <- guild_dfas1_24yr %>%
  filter(year >= 1998, year <= 2021) %>%
  select(year, x07_dfa_cpue_int_spr_jun_hw, x09_dfa_hake_age5plus, x16_sar) %>%
  drop_na()

fit_base <- sem('
  x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
  x16_sar ~ x07_dfa_cpue_int_spr_jun_hw
', data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE)

base_aic <- fitMeasures(fit_base, "aic")
base_r2  <- inspect(fit_base, "r2")[["x16_sar"]]

# ---------------------------
# 6) Build model sets--------
#    Simple global sweep:
#    - all 2-variable combos from all_cands
#    - all 3-variable combos from all_cands
# ---------------------------
# all 2-term and 3-term combinations
comb2 <- combn(all_cands, 2, simplify = FALSE)
comb3 <- combn(all_cands, 3, simplify = FALSE)

res_2 <- purrr::map_dfr(comb2, ~fit_sem_k(.x, guild_dfas1_24yr, "2var_global"))
res_3 <- purrr::map_dfr(comb3, ~fit_sem_k(.x, guild_dfas1_24yr, "3var_global"))

all_models <- bind_rows(res_2,res_3)

if (nrow(all_models) == 0) stop("No converged models found in 2- or 3-variable sweeps.")

# annotate model composition / index usage
all_models <- all_models %>%
  mutate(
    n_terms = if_else(is.na(term3), 2L, 3L),
    has_ishark = str_detect(terms, "ishark"),
    has_issl   = str_detect(terms, "issl"),
    has_icsl   = str_detect(terms, "icsl"),
    has_any_index = has_ishark | has_issl | has_icsl,
    index_count = as.integer(has_ishark) + as.integer(has_issl) + as.integer(has_icsl),
    model_label = paste0("M", row_number())
  ) %>%
  mutate(
    delta_aic_vs_base = aic - base_aic,
    delta_r2_vs_base  = sar_r2 - base_r2
  )

# ---------------------------
# 7) AIC weights + variable importance
# ---------------------------
all_models_aic <- all_models <- all_models %>%
  mutate(
    delta_aic_global = aic - min(aic, na.rm = TRUE),
    aic_weight_raw = exp(-0.5 * delta_aic_global),
    aic_weight = aic_weight_raw / sum(aic_weight_raw, na.rm = TRUE)
  ) %>%
  select(-aic_weight_raw)%>% 
  arrange(aic)

# long model-term table for importance
model_terms_long <- all_models %>%
  select(model_label, model_type, n_terms, terms, aic, aic_weight, term1, term2, term3,
         has_ishark, has_issl, has_icsl, has_any_index, index_count) %>%
  pivot_longer(cols = c(term1, term2, term3), names_to = "slot", values_to = "variable") %>%
  filter(!is.na(variable))

var_importance <- model_terms_long %>%
  group_by(variable) %>%
  summarise(
    importance = sum(aic_weight, na.rm = TRUE),
    n_models = n(),
    in_any_index_model_weight = sum(aic_weight[has_any_index], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(importance))

# optional variable tags
var_importance <- var_importance %>%
  mutate(
    var_group = case_when(
      str_detect(variable, "ishark|issl|icsl") ~ "Index",
      str_detect(variable, "^x10_|^x15_") ~ "Predator",
      str_detect(variable, "^x21_") ~ "Climate",
      str_detect(variable, "^x12_|^x13_") ~ "Prey",
      str_detect(variable, "^x14_") ~ "Competitor",
      TRUE ~ "Other"
    )
  )

# ---------------------------
# 8) Output tables
# ---------------------------
# sorted master outputs
# all_models_aic <- all_models %>% arrange(aic)
 all_models_r2  <- all_models %>% arrange(desc(sar_r2), aic)

write.csv(all_models,      file.path(out_dir, "sem_global_2var_models_full.csv"), row.names = FALSE)
# write.csv(all_models_aic,  file.path(out_dir, "sem_global_2var_models_sorted_aic.csv"), row.names = FALSE)
# write.csv(all_models_r2,   file.path(out_dir, "sem_global_2var_models_sorted_r2.csv"), row.names = FALSE)
write.csv(var_importance,  file.path(out_dir, "sem_global_2var_variable_importance_aicweights.csv"), row.names = FALSE)

# top subsets
top15_aic <- all_models_aic %>% slice_head(n = 15)
top15_r2  <- all_models_r2  %>% slice_head(n = 15)

write.csv(top15_aic, file.path(out_dir, "sem_global_top15_aic.csv"), row.names = FALSE)
write.csv(top15_r2,  file.path(out_dir, "sem_global_top15_r2.csv"), row.names = FALSE)

# ---------------------------
# 9) Variable importance plot
# ---------------------------
importance_plot_df <- var_importance %>%
  mutate(variable = fct_reorder(variable, importance))

p_imp <- ggplot(importance_plot_df, aes(x = variable, y = importance, fill = var_group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~var_group, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(
    title = "Variable importance from global 2- and 3-term SEM sweep",
    subtitle = "Importance = summed AIC weights across models containing variable",
    x = "Variable",
    y = "AIC-weighted importance",
    fill = NULL
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_imp
ggsave(
  file.path(out_dir, "sem_global_variable_importance_aicweights.png"),
  p_imp, width = 10, height = 12, dpi = 170
)

# ---------------------------
# 10) Console summary
# ---------------------------
cat("\n=== Baseline ===\n")
cat("Base AIC:", round(base_aic, 3), " | Base SAR R2:", round(base_r2, 3), "\n")

cat("\n=== Top 15 by AIC ===\n")
print(top15_aic[1:5,] %>%
        select(model_label, model_type, n_terms, terms, aic, delta_aic_vs_base, sar_r2, delta_r2_vs_base,
               has_ishark, has_issl, has_icsl, has_any_index, index_count))

cat("\n=== Top 15 by SAR R2 ===\n")
print(top15_r2 %>%
        select(model_label, model_type, n_terms, terms, sar_r2, delta_r2_vs_base, aic, delta_aic_vs_base,
               has_ishark, has_issl, has_icsl, has_any_index, index_count))
top15_r2[1:5,4:6]
cat("\n=== Top 20 variable importance ===\n")
print(var_importance %>% slice_head(n = 20))

top15_aic[1:5,"terms"]
top15_r2[1:5,"terms"]

message("\nWrote:")
message(" - ", file.path(out_dir, "sem_global_2var_3var_models_full.csv"))
message(" - ", file.path(out_dir, "sem_global_2var_3var_models_sorted_aic.csv"))
message(" - ", file.path(out_dir, "sem_global_2var_3var_models_sorted_r2.csv"))
message(" - ", file.path(out_dir, "sem_global_top15_aic.csv"))
message(" - ", file.path(out_dir, "sem_global_top15_r2.csv"))
message(" - ", file.path(out_dir, "sem_global_variable_importance_aicweights.csv"))
message(" - ", file.path(out_dir, "sem_global_variable_importance_aicweights.png"))



#better plot?------------
# ---------------------------
# 9) Plot variable importance
# (builds off your prior bar chart style, simplified to final importance)--------
# ---------------------------
cartoon_colors <- c(
  "Foundations"    = "#A3C1AD",
  "Data/Assembly"  = "#D9E4EC",
  "Functional"     = "#5D8AA8",
  "Simulation"     = "#E67E22",
  "Dimension/DFA"  = "#7F8C8D",
  "Selection/AIC"  = "#F39C12",
  "Latent/AK"      = "#2E4053",
  "Predator"       = "#2E4053",
  "Prey/Climate"   = "#D35400",
  "Other"          = "#7F8C8D"
)

plot_df <- var_importance %>%
  mutate(variable = fct_reorder(variable, importance))

p_imp <- ggplot(plot_df, aes(x = variable, y = importance, fill = group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = cartoon_colors, drop = FALSE) +
  theme_minimal() +
  labs(
    x = "Indicator",
    y = "AIC-weighted importance",
    fill = NULL,
    title = "Final variable importance across 2- and 3-variable SEM sweeps",
    subtitle = "Importance = sum of AIC weights across models containing variable"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "sem_variable_importance_aicweights_with_cslidx.png"), p_imp,
       width = 10, height = 12, dpi = 170)

# # ---------------------------
# # repeat? 10) Console summary---------
# # ---------------------------
# cat("\n=== Baseline ===\n")
# cat("Base AIC:", round(base_aic, 3), " | Base SAR R2:", round(base_r2, 3), "\n")
# 
# cat("\n=== Top 15 models by AIC ===\n")
# print(top15 %>% select(model_type, terms, aic, delta_aic_vs_base, sar_r2, delta_r2_vs_base, aic_weight))
# 
# cat("\n=== Top 20 variable importance ===\n")
# print(var_importance %>% slice_max(importance, n = 20))
# 
# message("\nWrote:")
# message(" - ", file.path(out_dir, "sem_sweep_2var_3var_with_cslidx_models.csv"))
# message(" - ", file.path(out_dir, "sem_top15_models_with_cslidx.csv"))
# message(" - ", file.path(out_dir, "sem_variable_importance_aicweights_with_cslidx.csv"))
# message(" - ", file.path(out_dir, "sem_variable_importance_aicweights_with_cslidx.png"))