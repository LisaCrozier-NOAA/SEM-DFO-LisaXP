
# Does size matter only in a "bad ocean"?

dat <- read.csv(file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv"), stringsAsFactors = FALSE);head(dat)

hist(dat$X16_DFA_sar_UCsprCh_smoltyr,breaks=200)

fit<-lm(X16_DFA_sar_UCsprCh_smoltyr~X06_DFA_IGF_mu_smoltyr,data=dat);summary(fit)
fit<-lm(X07_cpue_IntSprJunHW_smoltyr~X06_DFA_IGF_mu_smoltyr,data=dat);summary(fit)
fit<-lm(X07_cpue_IntSprJunHW_smoltyr~X06_DFA_IGF_mu_smoltyr,data=dat %>% filter(X16_DFA_sar_UCsprCh_smoltyr<0));summary(fit)
fit<-lm(X07_cpue_IntSprJunHW_smoltyr~X06_DFA_IGF_mu_smoltyr,data=dat %>% filter(X16_DFA_sar_UCsprCh_smoltyr< -0.5));summary(fit)


#Script: "Rcode_for_paper/02_analysis/01x_growth_test_dag1B.R"


#see if growth-sar is usually significant, and if aic is lower for a or b
library(tidyverse)

# ------------------------------------------------------------------------------
# 1. Define Paths to Doug's Full Sweep Outputs
# ------------------------------------------------------------------------------
dir_1a <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1A_long"
dir_1b <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1B_long"

res_1a <- read.csv(file.path(dir_1a, "SEMresultsByClus.csv"), stringsAsFactors = FALSE)
res_1b <- read.csv(file.path(dir_1b, "SEMresultsByClus.csv"), stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 2. Overall AIC Comparison: DAG 1A vs DAG 1B
# ------------------------------------------------------------------------------
cat("====================================================================\n")
cat("1. AIC COMPARISON: DAG 1A vs DAG 1B (Long, 24yr)\n")
cat("====================================================================\n")

cat(sprintf("Top Model AIC - DAG 1A: %.2f | DAG 1B: %.2f  --> (DAG 1B is %.2f AIC units lower)\n", 
            min(res_1a$AIC, na.rm = TRUE), 
            min(res_1b$AIC, na.rm = TRUE),
            min(res_1a$AIC, na.rm = TRUE) - min(res_1b$AIC, na.rm = TRUE)))

# Compare top 10 models for both
comp_top10 <- tibble(
  Rank = 1:10,
  `DAG 1A ModNum` = res_1a$modNum[1:10],
  `DAG 1A AIC`    = round(res_1a$AIC[1:10], 2),
  `DAG 1A AGFI` = round(res_1a$AGFI[1:10], 3),
  `DAG 1A p_Chi2` = round(res_1a$p_Chi2[1:10], 3),
  `DAG 1B ModNum` = res_1b$modNum[1:10],
  `DAG 1B AIC`    = round(res_1b$AIC[1:10], 2),
  `DAG 1B AGFI` = round(res_1b$AGFI[1:10], 3),
  `DAG 1B p_Chi2` = round(res_1b$p_Chi2[1:10], 3),
  `Delta AIC (1A - 1B)` = round(res_1a$AIC[1:10] - res_1b$AIC[1:10], 2)
)
print(as.data.frame(comp_top10))

# ------------------------------------------------------------------------------
# 3. Filter "Acceptable" Models for DAG 1B
# ------------------------------------------------------------------------------
min_aic_1b <- min(res_1b$AIC, na.rm = TRUE)

# Standard SEM acceptability criteria:
# Delta AIC <= 2.0 (or <= 4.0), p_Chi2 >= 0.05, CFI >= 0.90, SRMR <= 0.08
acc_1b_delta2 <- res_1b %>%
  mutate(deltaAIC = AIC - min_aic_1b) %>%
  filter(deltaAIC <= 2.0, p_Chi2 >= 0.05, CFI >= 0.90, SRMR <= 0.08)

acc_1b_delta4 <- res_1b %>%
  mutate(deltaAIC = AIC - min_aic_1b) %>%
  filter(deltaAIC <= 4.0, p_Chi2 >= 0.05, CFI >= 0.90, SRMR <= 0.08)

cat("\n====================================================================\n")
cat("2. ACCEPTABLE DAG 1B MODELS COUNT\n")
cat("====================================================================\n")
cat(sprintf("Number of DAG 1B models within Delta AIC <= 2.0 (and good fit): %d\n", nrow(acc_1b_delta2)))
cat(sprintf("Number of DAG 1B models within Delta AIC <= 4.0 (and good fit): %d\n", nrow(acc_1b_delta4)))

# ------------------------------------------------------------------------------
# 4. Check Significance & Magnitude of Growth -> SAR in Acceptable 1B Models
# ------------------------------------------------------------------------------
pe_1b_path <- file.path(dir_1b, "parameterEstimates.csv")

if (file.exists(pe_1b_path)) {
  pe_1b <- read.csv(pe_1b_path, stringsAsFactors = FALSE)
  
  # Growth -> SAR path in DAG 1B is coded as lhs == "logSAR", op == "~", rhs == "Growth"
  growth_sar_all <- pe_1b %>%
    filter(op == "~" & lhs %in% c("logSAR", "SAR") & rhs == "Growth")
  
  # Join with acceptable model sets
  growth_sar_delta2 <- growth_sar_all %>%
    filter(modNum %in% acc_1b_delta2$modNum)
  
  growth_sar_delta4 <- growth_sar_all %>%
    filter(modNum %in% acc_1b_delta4$modNum)
  
  cat("\n====================================================================\n")
  cat("3. 'Growth -> SAR' LINK EVALUATION IN ACCEPTABLE 1B MODELS\n")
  cat("====================================================================\n")
  
  summarize_growth_path <- function(df_path, label) {
    n_total <- nrow(df_path)
    if (n_total == 0) return(cat("No models to summarize for", label, "\n"))
    
    n_sig     <- sum(df_path$pvalue < 0.05, na.rm = TRUE)
    n_pos     <- sum(df_path$est > 0, na.rm = TRUE)
    n_pos_sig <- sum(df_path$est > 0 & df_path$pvalue < 0.05, na.rm = TRUE)
    
    cat(sprintf("\n--- %s (N = %d models) ---\n", label, n_total))
    cat(sprintf("Significant (p < 0.05): %d / %d (%.1f%%)\n", n_sig, n_total, (n_sig / n_total) * 100))
    cat(sprintf("Positive sign (+):      %d / %d (%.1f%%)\n", n_pos, n_total, (n_pos / n_total) * 100))
    cat(sprintf("Positive AND Sig:      %d / %d (%.1f%%)\n", n_pos_sig, n_total, (n_pos_sig / n_total) * 100))
    cat(sprintf("Mean Beta (std):       %.3f (Range: [%.3f, %.3f])\n", 
                mean(df_path$std.all, na.rm = TRUE), 
                min(df_path$std.all, na.rm = TRUE), 
                max(df_path$std.all, na.rm = TRUE)))
    cat(sprintf("Mean p-value:          %.4f\n", mean(df_path$pvalue, na.rm = TRUE)))
  }
  
  summarize_growth_path(growth_sar_delta2, "Delta AIC <= 2.0 (Strict)")
  summarize_growth_path(growth_sar_delta4, "Delta AIC <= 4.0 (Broad)")
  
} else {
  cat("\nparameterEstimates.csv not found in DAG1B directory; refitting loop will be needed.\n")
}
