#nice graph of coefficients!!!

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

out_dir <- "copilot/outputs_9"

# Load coefficients table and lookup table
coef_data <- read.csv(file.path(out_dir, "two_stage_interaction_coefficients.csv"))
pred_lookup <- read.csv(file.path(out_dir, "all_pred_dfa_altprey.csv")) %>% 
  clean_names() %>%
  mutate(across(everything(), tolower))

# Map species names back to int terms
map_prey_species <- function(df, lookup) {
  df %>%
    left_join(lookup, by = c("ak_predator" = "pred_data_col")) %>%
    mutate(
      species_name = case_when(
        term == "ak_int_1" ~ altprey1_data_col,
        term == "ak_int_2" ~ altprey2_data_col,
        term == "ak_int_3" ~ altprey3_data_col,
        term == "ncc_int_1" ~ altprey1_data_col,
        term == "ncc_int_2" ~ altprey2_data_col,
        term == "ncc_int_3" ~ altprey3_data_col,
        TRUE ~ term
      ),
      clean_species = case_when(
        grepl("pollock", species_name) ~ "Pollock (Age 1+)",
        grepl("hake", species_name)    ~ "Hake (Age 5+)",
        grepl("herr|sardine", species_name) ~ "Herring / Sardine",
        grepl("capelin", species_name) ~ "Capelin",
        grepl("krill", species_name)   ~ "Krill",
        TRUE ~ species_name
      )
    )
}

mapped_coefs <- map_prey_species(coef_data, pred_lookup) %>%
  filter(!is.na(est))

# Export Mapped Table
write.csv(mapped_coefs, file.path(out_dir, "mapped_interaction_coefficients_by_species.csv"), row.names = FALSE)

# Generate Visualization
p <- ggplot(mapped_coefs, aes(x = reorder(paste(ak_predator, clean_species, sep = " x "), est), y = est, fill = sign_check)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("Buffering (Negative)" = "#2b5c8f", "Positive (Check Overfit)" = "#d95f02")) +
  labs(
    title = "Alternate Prey Interaction Coefficients Across Top Two-Stage SEMs",
    subtitle = "Negative = True Prey Buffering | Positive = Apparent Competition / Overfitting",
    x = "Predator x Alternate Prey Interaction",
    y = "Standardized Regression Coefficient (Beta)",
    fill = "Effect Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p
ggsave(file.path(out_dir, "interaction_coefficients_by_species.png"), plot = p, width = 11, height = 8)

cat("Successfully mapped species and saved coefficient plot to copilot/outputs_9/interaction_coefficients_by_species.png!\n")