# topmod_x16<-read.csv(file.path(out_dir,"fair_ranked_models_x16_sar.csv")) %>% 
#   select(model_type, term1_col, term2_col, n_years, aic, delta_aic_vs_base, 
#          b_term1_sar, p_term1_sar, b_term2_sar, p_term2_sar, contains_new_shark) %>% 
#   head(5)
# 
# print(topmod_x16)
#   model_type                               term1_col                         term2_col n_years      aic delta_aic_vs_base b_term1_sar  p_term1_sar b_term2_sar  p_term2_sar contains_new_shark
# 1 1 Predator + 1 Prey                    x15_shark_enso_roll2                       x21_pdo_djf      24 105.7353         -16.53250  -0.5002842 1.533218e-05  -0.5755494 2.666867e-06               TRUE
# 2 1 Predator + 1 Prey                    x15_shark_enso_roll2 x21_sst_wgoa_coastwatch_junjulaug      24 107.0390         -15.22886  -0.3832021 8.481156e-04  -0.5959722 4.994543e-06               TRUE
# 3 1 Predator + 1 Prey                    x15_shark_enso_roll2    x12_dfa_biomass_euph_shelf_sum      24 108.7214         -13.54647  -0.5634725 1.220921e-05  -0.5353229 5.844474e-05               TRUE

#=================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
})

# -----------------------------------------------------------------------------
# 1. Extract Top Shark Index & Clean Names
# -----------------------------------------------------------------------------
      out_dir <- "copilot/outputs_4"
      
      # 'all_sims' contains the scenario output from the top model (enso_dj | z_roll2)
      
      all_sims<-read.csv(file.path(out_dir, "shark_predictions.csv"))
      
      top_shark_index <- all_sims %>%
        filter(Scenario == "enso_dj | z_roll2") %>%
        select(year, x15_shark_enso_roll2 = I_Shark)
      
      # Load, clean, scale guild dataset
      guild.dfas1 <- read.csv("data_Lisa/guild.dfasAK.csv", row.names = NULL) %>% 
        clean_names() %>% 
        select(-x10_harbor_seal_cr_2yr_lead) %>%
        left_join(top_shark_index, by = "year") %>%
        mutate(across(-year, ~ as.vector(scale(.))))

# -----------------------------------------------------------------------------
# 3. Identify Candidate Variable Sets
# -----------------------------------------------------------------------------
      # Temperature (x21_) and Prey candidates (x12_, x13_, x14_)
          ak_prey_cands <- names(guild.dfas1)[grepl("^x21_|^x12_|^x13_|^x14_", names(guild.dfas1))]
      
      # Predator candidates (x10_, x15_ including your new x15_shark_enso_roll2)
          ak_pred_cands <- names(guild.dfas1)[grepl("^x10_|^x15_", names(guild.dfas1))]
      
          message("Found ", length(ak_prey_cands), " Prey candidates and ", 
              length(ak_pred_cands), " Predator candidates (including new shark index).")

      # Helper function to fit lavaan SEM safely and pull path statistics
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
                
                if (nrow(df_model) < 12) return(NULL)
                
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
# 4A. SWEEP 1: Two Predator Combinations (Predator vs. Predator)
# -----------------------------------------------------------------------------
          # Generate unique pairs of predator candidates
          pred_pairs <- combn(ak_pred_cands, 2, simplify = FALSE)
          
          res_pred_pairs <- purrr::map_dfr(pred_pairs, function(pair) {
            fit_sem_two_topdown(pair[1], pair[2], guild.dfas1)
          }) %>%
            mutate(model_type = "2 Predators")

# -----------------------------------------------------------------------------
# 4B. SWEEP 2: One Predator + One Prey Combinations (Predator vs. Prey)
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
          # Fit baseline model (CPUE -> SAR only, no AK terms)
          df_base <- guild.dfas1 %>%
            filter(year >= 1998, year <= 2021) %>%
            filter(!is.na(x07_dfa_cpue_int_spr_jun_hw), !is.na(x09_dfa_hake_age5plus), !is.na(x16_sar))
          
          fit_base <- sem('
            x07_dfa_cpue_int_spr_jun_hw ~ x09_dfa_hake_age5plus
            x16_sar ~ x07_dfa_cpue_int_spr_jun_hw
          ', data = df_base, std.lv = TRUE, missing = "ML", warn = FALSE)
          
          base_aic <- fitMeasures(fit_base, "aic")
          
          # Combine all competitive model sweeps
          all_competing_models <- bind_rows(res_pred_pairs, res_pred_prey) %>%
            mutate(
              delta_aic_vs_base = aic - base_aic,
              contains_new_shark = (term1_col == "x15_shark_enso_roll2" | term2_col == "x15_shark_enso_roll2")
            ) %>%
            arrange(aic)

# -----------------------------------------------------------------------------
# 6. Display Model Selection Summary
# -----------------------------------------------------------------------------
          # Top 15 Overall Models
          cat("\n=== TOP 15 OVERALL SEM MODELS ===\n")
          print(
            all_competing_models %>% 
              select(model_type, term1_col, term2_col, aic, delta_aic_vs_base, 
                     b_term1_sar, p_term1_sar, b_term2_sar, p_term2_sar, contains_new_shark) %>% 
              head(15)
          )
          
          # Frequency of New Shark Index in Top Models
          shark_in_top20 <- sum(all_competing_models$contains_new_shark[1:20])
          cat("\nThe new shark index (x15_shark_enso_roll2) appears in", shark_in_top20, "out of the top 20 models.\n")
          
          
          topmod<- all_competing_models %>% 
            select(model_type, term1_col, term2_col, aic, delta_aic_vs_base, 
                   b_term1_sar, p_term1_sar, b_term2_sar, p_term2_sar, contains_new_shark) %>% 
            head(15)
          
          topmod[1,]
          