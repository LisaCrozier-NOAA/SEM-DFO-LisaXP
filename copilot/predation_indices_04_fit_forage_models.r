#results: fixed forage fish to  stka_herr_matbiom + mid_il_capelin


suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

if (!exists("out_dir")) out_dir <- "copilot/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# expects these objects from prior scripts:
# ak_raw, sem, shark, clim_tr(optional), guess_year_col(), find_col(), safe_scale()
# plus fixed salmon columns:
col_salmon_juv   <- "x07_dfa_cpue_int_spr_jun_hw"
col_salmon_adult <- "x16_sar"

# ---------------------------
# 1) Build annual tables
# ---------------------------
year_col_ak <- guess_year_col(ak_raw)
year_col_sem <- guess_year_col(sem)
year_col_shark <- guess_year_col(shark)

ak_yr <- ak_raw %>%
  rename(year = all_of(year_col_ak)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

sem_yr <- sem %>%
  rename(year = all_of(year_col_sem)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

shark_yr <- shark %>%
  rename(year = all_of(year_col_shark)) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

stopifnot(col_salmon_juv %in% names(sem_yr), col_salmon_adult %in% names(sem_yr))

# ---------------------------
# 2) Lock selected drivers
# ---------------------------
col_sst <- "sst_wgoa_coastwatch_junjulaug"
col_ssl <- "ssl_west_pup_pred"
use_climate <- FALSE  # fixed to T_sst

if (!col_sst %in% names(ak_yr)) stop("Locked SST column not found in ak_yr.")
if (!col_ssl %in% names(ak_yr)) stop("Locked SSL column not found in ak_yr.")

# shark candidates remain open
shark_cands <- find_col(shark_yr, c("shark","biomass","cpue","salmon_shark"))
if (length(shark_cands) == 0) shark_cands <- find_col(sem_yr, c("shark","salmon_shark"))
if (length(shark_cands) == 0) stop("No shark candidates found.")

# ---------------------------
# 3) Forage candidates (selective)
# ---------------------------
# user-vetted forage discovery
forage_ak  <- find_col(ak_raw, c("capelin", "sand_lance", "sand lance", "ammod", "herring", "herr"))
forage_sem <- find_col(sem %>% select(!contains("x05")), c("capelin", "sand_lance", "sand lance", "ammod", "herring", "forage"))

forage_ak_use  <- intersect(unique(forage_ak), names(ak_yr))
forage_sem_use <- intersect(unique(forage_sem), names(sem_yr))
forage_pool <- unique(c(forage_ak_use, forage_sem_use))

if (length(forage_pool) < 2) stop("Need at least 2 forage candidates in annual tables.")

# optional: keep only columns with enough non-NA in 1998-2021
valid_years <- 1998:2021
tmp_dat <- tibble(year = valid_years) %>%
  left_join(ak_yr %>% select(year, any_of(forage_pool)), by="year") %>%
  left_join(sem_yr %>% select(year, any_of(forage_pool)), by="year", suffix = c("_ak","_sem"))

# simplify by using original names presence from either source table:
enough_data <- function(colname) {
  v <- c(
    if (colname %in% names(ak_yr))  ak_yr %>% filter(year %in% valid_years) %>% pull(all_of(colname)) else NA_real_,
    if (colname %in% names(sem_yr)) sem_yr %>% filter(year %in% valid_years) %>% pull(all_of(colname)) else NA_real_
  )
  sum(!is.na(v)) >= 12
}
forage_pool <- forage_pool[purrr::map_lgl(forage_pool, enough_data)]

# make forage sets:
#  - singles
#  - pairs
#  - triples
make_sets <- function(x, k) combn(x, k, simplify = FALSE)
forage_sets <- c(
  make_sets(forage_pool, 1),
  if (length(forage_pool) >= 2) make_sets(forage_pool, 2) else list(),
  if (length(forage_pool) >= 3) make_sets(forage_pool, 3) else list()
)

# label
forage_set_tbl <- tibble(
  forage_set_id = seq_along(forage_sets),
  forage_cols = forage_sets,
  forage_label = purrr::map_chr(forage_sets, ~paste(.x, collapse = " + "))
)

write_csv(forage_set_tbl %>% select(forage_set_id, forage_label),
          file.path(out_dir, "phase2_forage_sets_tested.csv"))

# ---------------------------
# 4) Fit function with locked SST/SSL and variable forage set + shark
# ---------------------------
fit_locked_combo <- function(col_shark, forage_cols,
                             k0 = 0.0, kT = 1.0, kF = 1.0, Q10 = 2.0) {
  
  forage_cols_ak  <- intersect(forage_cols, names(ak_yr))
  forage_cols_sem <- intersect(forage_cols, names(sem_yr))
  forage_cols_all <- unique(c(forage_cols_ak, forage_cols_sem))
  if (length(forage_cols_all) == 0) return(NULL)
  
  data_selected <- tibble(year = sort(unique(c(ak_yr$year, sem_yr$year, shark_yr$year)))) %>%
    left_join(
      ak_yr %>% select(year, SST_raw = all_of(col_sst), SSL_raw = all_of(col_ssl), any_of(forage_cols_ak)),
      by = "year"
    ) %>%
    left_join(
      sem_yr %>% select(year, all_of(col_salmon_juv), all_of(col_salmon_adult), any_of(forage_cols_sem)),
      by = "year"
    ) %>%
    left_join(
      shark_yr %>% select(year, Shark_raw = all_of(col_shark)),
      by = "year"
    )
  
  forage_scaled <- data_selected %>%
    select(any_of(forage_cols_all)) %>%
    mutate(across(everything(), safe_scale))
  
  data_selected <- data_selected %>%
    mutate(F_raw = rowMeans(as.matrix(forage_scaled), na.rm = TRUE))
  
  T_ref <- mean(data_selected$SST_raw, na.rm = TRUE)
  
  data_selected <- data_selected %>%
    mutate(
      T_sst   = safe_scale(SST_raw),
      T_use   = T_sst,
      SSL_z   = safe_scale(SSL_raw),
      Shark_z = safe_scale(Shark_raw),
      F_z     = safe_scale(F_raw)
    )
  
  # switching function
  p_switch <- plogis(k0 + kT * data_selected$T_use - kF * data_selected$F_z)
  I_SSL <- data_selected$SSL_z * p_switch
  
  # shark index
  O_t <- pmax(0.2, 1 + 0.3 * data_selected$T_use)
  M_t <- Q10^((data_selected$SST_raw - T_ref)/10)
  I_Shark <- data_selected$Shark_z * M_t * O_t
  
  fit_dat <- data_selected %>%
    mutate(I_SSL = I_SSL, I_Shark = I_Shark) %>%
    filter(year >= 1998, year <= 2021) %>%
    filter(
      !is.na(.data[[col_salmon_adult]]),
      !is.na(.data[[col_salmon_juv]]),
      !is.na(I_SSL),
      !is.na(I_Shark)
    )
  
  if (nrow(fit_dat) < 12) return(NULL)
  
  mod <- lm(
    reformulate(c(col_salmon_juv, "I_SSL", "I_Shark"), response = col_salmon_adult),
    data = fit_dat
  )
  td <- broom::tidy(mod)
  sm <- summary(mod)
  
  get_term <- function(term, field) {
    out <- td %>% filter(.data$term == term) %>% pull(!!sym(field))
    if (length(out)==0) NA_real_ else out[1]
  }
  
  tibble(
    n = nrow(fit_dat),
    aic = AIC(mod),
    r2 = sm$r.squared,
    b_cpue = get_term(col_salmon_juv, "estimate"),
    b_ssl = get_term("I_SSL", "estimate"),
    b_shark = get_term("I_Shark", "estimate"),
    p_cpue = get_term(col_salmon_juv, "p.value"),
    p_ssl = get_term("I_SSL", "p.value"),
    p_shark = get_term("I_Shark", "p.value")
  )
}

# ---------------------------
# 5) Run grid: forage set x shark
# ---------------------------
grid <- tidyr::crossing(
  forage_set_id = forage_set_tbl$forage_set_id,
  col_shark = shark_cands
) %>%
  left_join(forage_set_tbl %>% select(forage_set_id, forage_cols, forage_label), by = "forage_set_id")

res <- purrr::pmap_dfr(
  grid,
  function(forage_set_id, col_shark, forage_cols, forage_label) {
    fit <- fit_locked_combo(col_shark = col_shark, forage_cols = forage_cols,
                            k0 = 0.0, kT = 1.0, kF = 1.0, Q10 = 2.0)
    if (is.null(fit)) return(tibble())
    bind_cols(
      tibble(
        sst_col = col_sst,
        forcing = "T_sst",
        ssl_col = col_ssl,
        shark_col = col_shark,
        forage_set_id = forage_set_id,
        forage_label = forage_label
      ),
      fit
    )
  }
)

if (nrow(res)==0) stop("No valid models fit in forage-selection sweep.")

res_ranked <- res %>% arrange(aic, desc(r2))
write_csv(res, file.path(out_dir, "phase2_forage_selection_all.csv"))
write_csv(res_ranked, file.path(out_dir, "phase2_forage_selection_ranked.csv"))

# ---------------------------
# 6) Importance of forage sets + individual forage variables
# ---------------------------
res_w <- res %>%
  mutate(model_group = case_when(n==19 ~ "short", n==24 ~ "long", TRUE ~ NA_character_)) %>%
  filter(!is.na(model_group)) %>%
  group_by(model_group) %>%
  mutate(
    delta_aic = aic - min(aic, na.rm=TRUE),
    rel_like = exp(-0.5 * delta_aic),
    aic_weight = rel_like / sum(rel_like, na.rm=TRUE)
  ) %>%
  ungroup()

# set-level importance
imp_sets <- res_w %>%
  group_by(model_group, forage_set_id, forage_label) %>%
  summarise(
    importance_aicw = sum(aic_weight, na.rm=TRUE),
    n_models = n(),
    min_aic = min(aic, na.rm=TRUE),
    mean_r2 = mean(r2, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  group_by(model_group) %>%
  mutate(rank = rank(-importance_aicw, ties.method="min"),
         importance_rel_top1 = importance_aicw / max(importance_aicw, na.rm=TRUE)) %>%
  ungroup() %>%
  arrange(model_group, rank)

write_csv(imp_sets, file.path(out_dir, "phase2_forage_set_importance_short_long.csv"))

# variable-level importance (sum weights across models whose forage set includes var)
expand_vars <- forage_set_tbl %>%
  select(forage_set_id, forage_cols) %>%
  tidyr::unnest_longer(forage_cols, values_to = "forage_var")

imp_vars <- res_w %>%
  select(model_group, forage_set_id, aic_weight) %>%
  left_join(expand_vars, by = "forage_set_id", relationship = "many-to-many") %>%
  group_by(model_group, forage_var) %>%
  summarise(
    importance_aicw = sum(aic_weight, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  group_by(model_group) %>%
  mutate(rank = rank(-importance_aicw, ties.method="min"),
         importance_rel_top1 = importance_aicw / max(importance_aicw, na.rm=TRUE)) %>%
  ungroup() %>%
  arrange(model_group, rank)

write_csv(imp_vars, file.path(out_dir, "phase2_forage_variable_importance_short_long.csv"))

# quick plots
p_sets <- ggplot(imp_sets %>% group_by(model_group) %>% slice_head(n=15),
                 aes(x=reorder(forage_label, importance_aicw), y=importance_aicw, fill=model_group)) +
  geom_col() + coord_flip() + theme_bw() +
  labs(title="Top forage sets by AIC-weighted importance", x="Forage set", y="Summed AIC weight")
ggsave(file.path(out_dir, "phase2_forage_set_importance_top15.png"), p_sets, width=12, height=8, dpi=150)

p_sets$data$forage_label[1:10]
# [1] "stka_herr_matbiom + capelin_w_go_a + est_abund_herring_recruits"      
# [2] "stka_herr_matbiom + mid_il_capelin + est_abund_herring_recruits"      
# [3] "stka_herr_matbiom + sitka_herring_e_go_a + est_abund_herring_recruits"
# [4] "stka_herr_matbiom + crg_herr_matbiom + est_abund_herring_recruits"    
# [5] "stka_herr_matbiom + capelin_w_go_a"                                   
# [6] "mid_il_capelin + sitka_herring_e_go_a + est_abund_herring_recruits"   
# [7] "stka_herr_matbiom + crg_herr_matbiom + capelin_w_go_a"                
# [8] "stka_herr_matbiom + capelin_w_go_a + sitka_herring_e_go_a"            
# [9] "stka_herr_matbiom + mid_il_capelin"                                   
# [10] "capelin_w_go_a + sitka_herring_e_go_a + est_abund_herring_recruits" 


p_vars <- ggplot(imp_vars, aes(x=reorder(forage_var, importance_aicw), y=importance_aicw, fill=model_group)) +
  geom_col(position="dodge") + coord_flip() + theme_bw() +
  labs(title="Forage variable importance (AIC-weighted)", x="Forage variable", y="Summed AIC weight")
ggsave(file.path(out_dir, "phase2_forage_variable_importance.png"), p_vars, width=10, height=7, dpi=150)

message("Done.\n",
        "- phase2_forage_selection_all.csv\n",
        "- phase2_forage_selection_ranked.csv\n",
        "- phase2_forage_set_importance_short_long.csv\n",
        "- phase2_forage_variable_importance_short_long.csv\n")