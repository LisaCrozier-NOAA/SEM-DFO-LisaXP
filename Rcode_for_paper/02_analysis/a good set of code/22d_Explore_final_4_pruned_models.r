

library(lavaan)

# save_rdata    <- file.path(proj_dir, "Routput_for_paper/models_DAG1B_latent_stages.RData")
# save(dag1b_candidate_models, df_master, file = save_rdata)
 load(file.path(proj_dir, "Routput_for_paper/models_DAG1B_latent_stages.RData"),verbose=TRUE)

top4_keys <- c("Model_1", "Model_3", "Model_4", "Model_5")

# Run sem_output_summary_fxn on each final pruned model
for (m in top4_keys) {
  cat("\n====================================================================\n")
  cat("MODEL SUMMARY:", m, "\n")
  cat("====================================================================\n")
  fit_obj <- dag1b_candidate_models[[m]]$Stage3_Final_Pruned_Regression$fit
  
  # if (exists("sem_output_summary_fxn")) {
  #   sem_output_summary_fxn(fit_obj)
  # } else {
    print(summary(fit_obj, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE))
  # }
}


fit_obj <- dag1b_candidate_models[["Model_1"]]$Stage3_Final_Pruned_Regression$fit

 summary(fit_obj, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
 sem_output_summary_fxn(fit_obj)
 