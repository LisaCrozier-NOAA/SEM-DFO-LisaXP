suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

if (!exists("out_dir")) out_dir <- "copilot/outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Expected objects in env:
# ak_yr, sem_yr, shark_yr, safe_scale

# -------------------------
# Final fixed settings
# -------------------------
col_salmon_juv   <- "x07_dfa_cpue_int_spr_jun_hw"
col_salmon_adult <- "x16_sar"

# SSL side (fixed)
col_sst_ssl <- "sst_wgoa_coastwatch_junjulaug"
col_ssl <- "ssl_west_pup_pred"
forage_fixed <- c("stka_herr_matbiom", "mid_il_capelin")

k0 <- 0.0
kT <- 1.5
kF <- 2.0

# Shark side (fixed except species)
shark_models <- c("pacific_sleeper_shark_goa", "salmon_shark_goa")
shark_temp_col <- "pdo_djf"
Q10 <- 2.0
overlap_slope <- 0.4
shark_transform <- "z"   # simplest
overlap_form <- "linear"

# -------------------------
# Checks
# -------------------------
stopifnot(
  col_sst_ssl %in% names(ak_yr),
  col_ssl %in% names(ak_yr),
  shark_temp_col %in% names(ak_yr),
  col_salmon_juv %in% names(sem_yr),
  col_salmon_adult %in% names(sem_yr)
)

missing_sharks <- setdiff(shark_models, names(shark_yr))
if (length(missing_sharks) > 0) {
  stop("Missing shark columns in shark_yr: ", paste(missing_sharks, collapse = ", "))
}

forage_ak_use <- intersect(forage_fixed, names(ak_yr))
forage_sem_use <- intersect(forage_fixed, names(sem_yr))
forage_cols_all <- unique(c(forage_ak_use, forage_sem_use))
if (length(forage_cols_all) == 0) stop("No fixed forage columns found in ak_yr/sem_yr.")

# -------------------------
# Build base data + fixed SSL pieces
# -------------------------
data_base <- tibble(year = sort(unique(c(ak_yr$year, sem_yr$year, shark_yr$year)))) %>%
  left_join(
    ak_yr %>% select(year, SST_ssl = all_of(col_sst_ssl), SSL_raw = all_of(col_ssl),
                     SharkTemp_raw = all_of(shark_temp_col), any_of(forage_ak_use)),
    by = "year"
  ) %>%
  left_join(
    sem_yr %>% select(year, all_of(col_salmon_juv), all_of(col_salmon_adult), any_of(forage_sem_use)),
    by = "year"
  )

forage_scaled <- data_base %>%
  select(any_of(forage_cols_all)) %>%
  mutate(across(everything(), safe_scale))

data_base <- data_base %>%
  mutate(
    F_raw = rowMeans(as.matrix(forage_scaled), na.rm = TRUE),
    T_ssl = safe_scale(SST_ssl),
    SSL_z = safe_scale(SSL_raw),
    F_z   = safe_scale(F_raw)
  )

p_switch_ssl <- plogis(k0 + kT * data_base$T_ssl - kF * data_base$F_z)
I_SSL_fixed <- data_base$SSL_z * p_switch_ssl

# -------------------------
# Fit function for one shark species
# -------------------------
fit_one_final <- function(shark_col) {
  d <- data_base %>%
    left_join(shark_yr %>% select(year, Shark_raw = all_of(shark_col)), by = "year")
  
  Shark_use <- safe_scale(d$Shark_raw)
  
  T_shark_raw <- d$SharkTemp_raw
  T_shark_z <- safe_scale(T_shark_raw)
  Tref <- mean(T_shark_raw, na.rm = TRUE)
  
  M_t <- Q10^((T_shark_raw - Tref)/10)
  O_t <- pmax(0.1, 1 + overlap_slope * T_shark_z)
  I_Shark <- Shark_use * M_t * O_t
  
  fit_dat <- d %>%
    mutate(
      p_switch_ssl = p_switch_ssl,
      I_SSL = I_SSL_fixed,
      Shark_z = Shark_use,
      T_shark_z = T_shark_z,
      M_shark = M_t,
      O_shark = O_t,
      I_Shark = I_Shark
    ) %>%
    filter(year >= 1998, year <= 2021) %>%
    filter(
      !is.na(.data[[col_salmon_adult]]),
      !is.na(.data[[col_salmon_juv]]),
      !is.na(I_SSL),
      !is.na(I_Shark)
    )
  
  mod <- lm(
    reformulate(c(col_salmon_juv, "I_SSL", "I_Shark"), response = col_salmon_adult),
    data = fit_dat
  )
  mod_base <- lm(
    reformulate(c(col_salmon_juv, "I_SSL"), response = col_salmon_adult),
    data = fit_dat
  )
  
  list(
    shark_col = shark_col,
    fit_dat = fit_dat,
    mod = mod,
    mod_base = mod_base
  )
}

# -------------------------
# Run final models
# -------------------------
fits <- lapply(shark_models, fit_one_final)

# coefficient summary
coef_tbl <- purrr::map_dfr(fits, function(x) {
  td <- broom::tidy(x$mod)
  td %>%
    mutate(shark_col = x$shark_col,
           n = nrow(x$fit_dat),
           aic = AIC(x$mod),
           aic_base = AIC(x$mod_base),
           delta_aic_shark = AIC(x$mod) - AIC(x$mod_base),
           r2 = summary(x$mod)$r.squared)
})

write_csv(coef_tbl, file.path(out_dir, "final_models_coefficients_goa_sharks.csv"))

# model-level summary
model_tbl <- coef_tbl %>%
  group_by(shark_col, n, aic, aic_base, delta_aic_shark, r2) %>%
  summarise(
    p_cpue = p.value[term == col_salmon_juv][1],
    p_ssl = p.value[term == "I_SSL"][1],
    p_shark = p.value[term == "I_Shark"][1],
    b_cpue = estimate[term == col_salmon_juv][1],
    b_ssl = estimate[term == "I_SSL"][1],
    b_shark = estimate[term == "I_Shark"][1],
    .groups = "drop"
  ) %>%
  mutate(
    sig_cpue = p_cpue < 0.05,
    sig_ssl = p_ssl < 0.05,
    sig_shark = p_shark < 0.05
  )

write_csv(model_tbl, file.path(out_dir, "final_models_summary_goa_sharks.csv"))
print(model_tbl)

# -------------------------
# Graph helpers
# -------------------------
plot_partial <- function(fit_obj, focal = c("I_SSL","I_Shark")) {
  focal <- match.arg(focal)
  dat <- fit_obj$fit_dat
  mod <- fit_obj$mod
  
  med_cpue <- median(dat[[col_salmon_juv]], na.rm = TRUE)
  med_ssl <- median(dat$I_SSL, na.rm = TRUE)
  med_shark <- median(dat$I_Shark, na.rm = TRUE)
  
  if (focal == "I_SSL") {
    nd <- tibble(
      !!col_salmon_juv := med_cpue,
      I_SSL = seq(min(dat$I_SSL, na.rm = TRUE), max(dat$I_SSL, na.rm = TRUE), length.out = 200),
      I_Shark = med_shark
    )
    nd$pred <- predict(mod, newdata = nd)
    
    ggplot(dat, aes(x = I_SSL, y = .data[[col_salmon_adult]])) +
      geom_point(alpha = 0.8) +
      geom_line(data = nd, aes(x = I_SSL, y = pred), color = "blue3", linewidth = 1.1) +
      theme_bw() +
      labs(
        title = paste0("SAR response vs I_SSL (", fit_obj$shark_col, ")"),
        subtitle = "Partial prediction: CPUE and I_Shark held at medians",
        x = "I_SSL",
        y = col_salmon_adult
      )
  } else {
    nd <- tibble(
      !!col_salmon_juv := med_cpue,
      I_SSL = med_ssl,
      I_Shark = seq(min(dat$I_Shark, na.rm = TRUE), max(dat$I_Shark, na.rm = TRUE), length.out = 200)
    )
    nd$pred <- predict(mod, newdata = nd)
    
    ggplot(dat, aes(x = I_Shark, y = .data[[col_salmon_adult]])) +
      geom_point(alpha = 0.8, color = "firebrick3") +
      geom_line(data = nd, aes(x = I_Shark, y = pred), color = "firebrick4", linewidth = 1.1) +
      theme_bw() +
      labs(
        title = paste0("SAR response vs I_Shark (", fit_obj$shark_col, ")"),
        subtitle = "Partial prediction: CPUE and I_SSL held at medians",
        x = "I_Shark",
        y = col_salmon_adult
      )
  }
}

plot_switch <- function(fit_obj) {
  dat <- fit_obj$fit_dat
  ggplot(dat, aes(x = F_z, y = p_switch_ssl, color = T_ssl)) +
    geom_point(size = 2, alpha = 0.85) +
    scale_color_viridis_c() +
    theme_bw() +
    labs(
      title = paste0("SSL prey switching (", fit_obj$shark_col, ")"),
      subtitle = paste0("p_switch = logistic(", k0, " + ", kT, "*T_ssl - ", kF, "*F_z)"),
      x = "F_z (forage index, scaled)",
      y = "p_switch_ssl",
      color = "T_ssl"
    )
}

plot_shark_overlap <- function(fit_obj) {
  dat <- fit_obj$fit_dat
  ggplot(dat, aes(x = T_shark_z, y = O_shark)) +
    geom_point(alpha = 0.8, color = "darkgreen") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    theme_bw() +
    labs(
      title = paste0("Shark overlap function O(T) (", fit_obj$shark_col, ")"),
      subtitle = paste0("O = max(0.1, 1 + ", overlap_slope, "*T_shark_z)"),
      x = "T_shark_z (PDO standardized)",
      y = "O_shark"
    )
}

# -------------------------
# Save plots
# -------------------------
for (f in fits) {
  nm <- f$shark_col
  
  p1 <- plot_partial(f, "I_SSL")
  p2 <- plot_partial(f, "I_Shark")
  p3 <- plot_switch(f)
  p4 <- plot_shark_overlap(f)
  
  ggsave(file.path(out_dir, paste0("final_", nm, "_sar_vs_I_SSL.png")), p1, width = 7.5, height = 5, dpi = 160)
  ggsave(file.path(out_dir, paste0("final_", nm, "_sar_vs_I_Shark.png")), p2, width = 7.5, height = 5, dpi = 160)
  ggsave(file.path(out_dir, paste0("final_", nm, "_ssl_switch.png")), p3, width = 7.5, height = 5, dpi = 160)
  ggsave(file.path(out_dir, paste0("final_", nm, "_shark_overlap.png")), p4, width = 7.5, height = 5, dpi = 160)
}

message("Done. Wrote final model tables + plots to ", out_dir)