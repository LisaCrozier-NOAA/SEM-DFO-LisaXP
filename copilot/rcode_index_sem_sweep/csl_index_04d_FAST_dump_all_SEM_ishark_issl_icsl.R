# ==============================================================================
# Fixed Subset Processor: Handles 2-Var & 3-Var Dataframes & Corrects Index Plotting
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

out_dir <- "copilot/outputs_7"

# 1. Load Pre-Calculated Models & Clean Joined Column Names safely (.x / .y)
# Reads fast_2var_all_models.csv (or sem_global_2var_models_full.csv if available)
model_file <- file.path(out_dir, "fast_2var_all_models.csv")
if (!file.exists(model_file)) {
  model_file <- file.path(out_dir, "sem_global_2var_models_full.csv")
}

all_models <- read.csv(model_file, stringsAsFactors = FALSE) %>%
  mutate(
    across(any_of(c("term1", "term2", "term3")), ~ str_remove(., "\\.[xy]$")),
    terms = str_remove_all(terms, "\\.[xy]")
  )

# 2. Corrected Directional Expectation & Group Classification
get_expected_sign <- function(v) {
  if (is.na(v) || v == "") return(NA_character_)
  case_when(
    str_detect(v, "(?i)ishark|issl|icsl") ~ "Negative",
    str_detect(v, "(?i)^x10_|^x15_")       ~ "Negative",
    str_detect(v, "(?i)^x21_")              ~ "Negative",
    str_detect(v, "(?i)^x12_|^x13_")        ~ "Positive",
    str_detect(v, "(?i)^x14_")              ~ "Negative",
    TRUE                                   ~ "Unknown"
  )
}

# Ensure Index terms are checked FIRST before general x15_ predator matches
get_var_group <- function(v) {
  case_when(
    str_detect(v, "(?i)ishark|issl|icsl") ~ "Index",
    str_detect(v, "(?i)^x10_|^x15_")       ~ "Predator",
    str_detect(v, "(?i)^x21_")              ~ "Climate",
    str_detect(v, "(?i)^x12_|^x13_")        ~ "Prey",
    str_detect(v, "(?i)^x14_")              ~ "Competitor",
    TRUE                                   ~ "Other"
  )
}

check_sign <- function(b, exp) {
  if (is.na(b) || is.na(exp) || exp == "Unknown") return(NA)
  if (exp == "Positive") return(b > 0)
  if (exp == "Negative") return(b < 0)
  return(NA)
}

# 3. Master Processing Function
process_subset <- function(df_sub, plotname, mytitle) {
  cat("\n====================================================================\n")
  cat(" PROCESSING SUBSET:", toupper(plotname), "\n")
  cat(" Title:", mytitle, "\n")
  cat("====================================================================\n")
  
  # Recalculate Weights & Sign Compliance
  df_sub <- df_sub %>%
    mutate(
      delta_aic_global = aic - min(aic, na.rm = TRUE),
      aic_weight       = exp(-0.5 * delta_aic_global) / sum(exp(-0.5 * delta_aic_global), na.rm = TRUE),
      exp1             = sapply(term1, get_expected_sign),
      exp2             = sapply(term2, get_expected_sign),
      pass1            = mapply(check_sign, b1, exp1),
      pass2            = mapply(check_sign, b2, exp2),
      prop_signs_met   = (replace_na(pass1, FALSE) + replace_na(pass2, FALSE)) / 2
    ) %>%
    arrange(aic)
  
  # Long Format Evaluation
  terms_long <- bind_rows(
    df_sub %>% select(aic_weight, variable = term1, estimate = b1, exp_sign = exp1, sign_met = pass1),
    df_sub %>% select(aic_weight, variable = term2, estimate = b2, exp_sign = exp2, sign_met = pass2)
  ) %>%
    filter(!is.na(variable)) %>%
    mutate(var_group = sapply(variable, get_var_group))
  
  # Calculate Variable Importance
  var_imp <- terms_long %>%
    group_by(variable, var_group, exp_sign) %>%
    summarise(
      importance            = sum(aic_weight, na.rm = TRUE),
      times_in_models       = n(),
      prop_sign_met_aic_wtd = sum(aic_weight[sign_met == TRUE], na.rm = TRUE) / sum(aic_weight, na.rm = TRUE),
      .groups               = "drop"
    ) %>%
    arrange(desc(importance))
  
  cat("\n--- TOP 10 VARIABLES IN THIS SUBSET ---\n")
  print(
    var_imp %>% 
      slice_head(n = 10) %>%
      transmute(
        Variable = variable,
        Group = var_group,
        Expected = exp_sign,
        `AIC Importance` = round(importance, 4),
        `Sign Match %` = paste0(round(prop_sign_met_aic_wtd * 100, 1), "%")
      )
  )
  
  # Generate Corrected Plot
  plot_df <- var_imp %>% mutate(variable = fct_reorder(variable, importance))
  
  p_imp <- ggplot(plot_df, aes(x = variable, y = importance, fill = var_group)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~var_group, scales = "free_y", ncol = 1) +
    scale_fill_manual(
      values = c(
        "Index"      = "#2563eb", # Bright Blue
        "Predator"   = "#dc2626", # Red
        "Prey"       = "#059669", # Green
        "Climate"    = "#d97706", # Orange
        "Competitor" = "#7c3aed"  # Purple
      ),
      drop = FALSE
    ) +
    theme_minimal(base_size = 11) +
    labs(
      title = mytitle,
      subtitle = "Importance = Summed AIC weights | Corrected Index Faceting",
      x = "Variable",
      y = "AIC-Weighted Importance",
      fill = NULL
    ) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  out_filename <- file.path(out_dir, paste0("sem_var_importance_", plotname, "_fixed.png"))
  ggsave(out_filename, p_imp, width = 10, height = 11, dpi = 170)
  
  return(list(var_importance = var_imp, plot = p_imp))
}

# ==============================================================================
# RUN & DISPLAY ALL PLOTS
# ==============================================================================

sub1 <- process_subset(all_models, "alldrivers_2var", "2 Variables: All AK Drivers")
sub2 <- process_subset(all_models %>% filter(!grepl("x21", terms)), "noclimate_2var", "2 Variables: No Climate")
sub3 <- process_subset(all_models %>% filter(!grepl("x21", terms) & !grepl("icsl", terms)), "noclimate_noratio_2var", "2 Variables: No Climate & No CSL Ratio")

# Display plots directly to viewer
print(sub1$plot)
print(sub2$plot)
print(sub3$plot)