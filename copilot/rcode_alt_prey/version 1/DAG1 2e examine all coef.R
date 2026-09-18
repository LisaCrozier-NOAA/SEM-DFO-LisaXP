suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(ggplot2)
})

out_dir <- "copilot/outputs_9"

# -----------------------------------------------------------------------------
# 1. Load Data & Lookup Tables
# -----------------------------------------------------------------------------
sem_complete_data <- read.csv(file.path(out_dir, "sem_altprey_data_complete_1998_2021.csv"), row.names = NULL)
names(sem_complete_data) <- tolower(names(sem_complete_data))

all_two_stage_models <- read.csv(file.path(out_dir, "two_stage_factorial_sem_sweep.csv"))

ncc_candidates <- read.csv(file.path(out_dir, "all_pred_dfa_altprey.csv")) %>% 
  clean_names() %>% mutate(across(everything(), tolower)) %>%
  filter(region == "ncc" & grepl("^(x08_|x09_|x11_)", pred_data_col))

ak_candidates <- read.csv(file.path(out_dir, "all_pred_dfa_altprey.csv")) %>% 
  clean_names() %>% mutate(across(everything(), tolower)) %>%
  filter(region == "ak" | grepl("^x10_", pred_data_col))

# -----------------------------------------------------------------------------
# 2. Comprehensive Parameter Extraction (Main Effects + Interactions)
# -----------------------------------------------------------------------------
extract_full_model_parameters <- function(top_models_df, raw_data, ncc_lookup, ak_lookup) {
  param_records <- list()
  
  for (i in seq_len(nrow(top_models_df))) {
    row <- top_models_df[i, ]
    ncc_pred <- row$ncc_predator
    ak_pred  <- row$ak_predator
    
    df_m <- raw_data
    
    # Re-build NCC interaction terms
    ncc_row <- ncc_lookup %>% filter(pred_data_col == ncc_pred) %>% slice(1)
    ncc_preys <- c(ncc_row$altprey1_data_col, ncc_row$altprey2_data_col, ncc_row$altprey3_data_col)
    ncc_preys <- ncc_preys[!is.na(ncc_preys) & ncc_preys != "" & ncc_preys != "na"]
    ncc_ints <- c()
    for (k in seq_along(ncc_preys)) {
      pc <- ncc_preys[k]
      if (pc %in% names(df_m)) {
        in_name <- paste0("ncc_int_", k)
        df_m[[in_name]] <- df_m[[ncc_pred]] * df_m[[pc]]
        ncc_ints <- c(ncc_ints, in_name)
      }
    }
    
    # Re-build AK interaction terms
    ak_row <- ak_lookup %>% filter(pred_data_col == ak_pred) %>% slice(1)
    ak_preys <- c(ak_row$altprey1_data_col, ak_row$altprey2_data_col, ak_row$altprey3_data_col)
    ak_preys <- ak_preys[!is.na(ak_preys) & ak_preys != "" & ak_preys != "na"]
    ak_ints <- c()
    for (k in seq_along(ak_preys)) {
      pc <- ak_preys[k]
      if (!pc %in% names(df_m)) {
        bp <- str_remove(pc, "_2yrlead$")
        if (bp %in% names(df_m)) df_m[[pc]] <- dplyr::lead(df_m[[bp]], 2)
      }
      if (pc %in% names(df_m)) {
        in_name <- paste0("ak_int_", k)
        df_m[[in_name]] <- df_m[[ak_pred]] * df_m[[pc]]
        ak_ints <- c(ak_ints, in_name)
      }
    }
    
    # Active terms in pruned model
    active_ncc <- unlist(strsplit(row$active_ncc_ints, ",\\s*"))
    active_ncc <- active_ncc[active_ncc %in% ncc_ints]
    
    active_ak <- unlist(strsplit(row$active_ak_ints, ",\\s*"))
    active_ak <- active_ak[active_ak %in% ak_ints]
    
    # Build syntax
    rhs_ncc <- paste(c(ncc_pred, active_ncc), collapse = " + ")
    rhs_ak  <- paste(c("x07_dfa_cpue_intsprjunhw", ak_pred, active_ak), collapse = " + ")
    
    sem_syntax <- paste0(
      "x07_dfa_cpue_intsprjunhw ~ ", rhs_ncc, "\n",
      "x16_sar ~ ", rhs_ak
    )
    
    fit <- tryCatch(
      sem(sem_syntax, data = df_m, std.lv = TRUE, missing = "ML", warn = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(fit) && lavInspect(fit, "converged")) {
      pe <- parameterEstimates(fit) %>%
        filter(op == "~") %>%
        select(lhs, term = rhs, est, se, z, pvalue) %>%
        mutate(
          model_rank   = i,
          ncc_predator = ncc_pred,
          ak_predator  = ak_pred,
          aic          = row$aic,
          term_type    = case_when(
            term %in% c(ncc_pred, ak_pred) ~ "Main Predator Effect",
            term == "x07_dfa_cpue_intsprjunhw" ~ "Early Survival Link (CPUE)",
            TRUE ~ "Alternate Prey Interaction"
          ),
          effect_direction = case_when(
            est < 0 & term_type == "Main Predator Effect" ~ "Top-Down Predation (-)",
            est > 0 & term_type == "Main Predator Effect" ~ "Co-variance / Bottom-Up (+)",
            est < 0 & term_type == "Alternate Prey Interaction" ~ "Buffering (-)",
            est > 0 & term_type == "Alternate Prey Interaction" ~ "Apparent Competition (+)",
            TRUE ~ "Positive Link"
          )
        )
      
      param_records[[i]] <- pe
    }
  }
  
  bind_rows(param_records)
}

# -----------------------------------------------------------------------------
# 3. Execute & Export
# -----------------------------------------------------------------------------
top_valid_models <- all_two_stage_models %>%
  filter(model_type == "Pruned", pvalue >= 0.05, cfi >= 0.90) %>%
  arrange(aic) %>%
  slice_head(n = 15)

full_param_table <- extract_full_model_parameters(
  top_valid_models, 
  sem_complete_data, 
  ncc_candidates, 
  ak_candidates
) %>%
  as_tibble()

write.csv(full_param_table, file.path(out_dir, "full_model_parameters_main_and_interactions.csv"), row.names = FALSE)

# Display Summary Table
cat("\n====================================================================\n")
cat(" MAIN PREDATOR EFFECTS & INTERACTION TERMS FOR TOP VALID MODELS      \n")
cat("====================================================================\n")

print(
  full_param_table %>%
    filter(term_type != "Early Survival Link (CPUE)") %>%
    select(model_rank, ak_predator, term, term_type, est, pvalue, effect_direction),
  n = Inf
)


print(
  full_param_table %>%
    filter(term_type != "Early Survival Link (CPUE)") %>%
    filter(term_type == "Main Predator Effect") %>%
  #  filter(est < 0) %>%
#    filter(pvalue > 0.05 & est > 0) %>%
#      filter(pvalue < 0.05 & est > 0) %>%
  select(model_rank, ak_predator, term, term_type, est, pvalue, effect_direction),
  n = Inf
)

head(full_param_table)

print(
  full_param_table %>%
    filter(term_type != "Early Survival Link (CPUE)") %>%
  #  filter(model_rank < 5) %>%
    filter(term == "ak_int_1") %>%
    #    filter(term_type == "Main Predator Effect") %>%
    #  filter(est < 0) %>%
    #    filter(pvalue > 0.05 & est > 0) %>%
    #      filter(pvalue < 0.05 & est > 0) %>%
    select(model_rank, ak_predator, term, term_type, est, pvalue, effect_direction),
  n = Inf
)

#exclude models 7 and 11 because they had strong positive predator main effects
# model_rank ak_predator                  term                         term_type              est     pvalue effect_direction           
# 1          7 x11_ssl_seak_pup_pred        x11_ssl_seak_pup_pred        Main Predator Effect 0.641 0.00000167 Co-variance / Bottom-Up (+)
# 2         11 x15_pacificcodbiomass_predak x15_pacificcodbiomass_predak Main Predator Effect 0.448 0.0230     Co-variance / Bottom-Up (+)

print(
  full_param_table %>%
    filter(term_type != "Early Survival Link (CPUE)") %>%
    filter(!model_rank %in% c(7,11)) %>%
   # filter(term == "ak_int_1") %>%
        filter(term_type == "Main Predator Effect") %>%
    #  filter(est < 0) %>%
    #    filter(pvalue > 0.05 & est > 0) %>%
    #      filter(pvalue < 0.05 & est > 0) %>%
    select(model_rank, ak_predator, term, term_type, est, pvalue, effect_direction),
  n = Inf
)
