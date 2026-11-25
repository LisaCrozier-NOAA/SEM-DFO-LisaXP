# ==============================================================================
# Script: 02_analysis/03b_export_peff_display_tables.R
# Purpose: Build ready-to-print Quarto display CSVs from saved Peff sweep CSVs
# ==============================================================================

library(tidyverse)

proj_dir       <- file.path(getwd(), "Rcode_for_paper")
data_out_dir   <- file.path(proj_dir, "Outputs", "data")
crosswalk_path <- file.path(proj_dir, "metadata", "master_name_crosswalk.csv")

# 1. Load saved sweep summaries (no model re-running needed)
res_24 <- read_csv(file.path(data_out_dir, "obj3_data_peff_sweep_all_models_24yr.csv"), show_col_types = FALSE)
res_19 <- read_csv(file.path(data_out_dir, "obj3_data_peff_sweep_all_models_19yr.csv"), show_col_types = FALSE)

# 2. Load predator name crosswalk
lookup_dict <- character(0)
if (file.exists(crosswalk_path)) {
  cw <- read_csv(crosswalk_path, show_col_types = FALSE)
  lookup_dict <- c(setNames(cw$LisaName, cw$dfa_cols), setNames(cw$LisaName, cw$DFAname))
}

shorten_pred_label <- function(x) {
  x <- gsub("^(X08_|X09_|X10_|X11_|X15_)", "", x)
  x <- gsub("DFA_", "", x)
  x <- gsub("_predAK", "", x)
  x <- gsub("Biomass", "", x)
  ifelse(nchar(x) > 24, paste0(substr(x, 1, 22), ".."), x)
}

# 3. Format & export clean display CSVs
export_peff_display_csv <- function(sweep_df, out_csv_path) {
  sweep_df %>%
    mutate(
      raw_mod_variant = mod_variant,
      Predator_Clean = case_when(
        predator %in% names(lookup_dict) ~ shorten_pred_label(lookup_dict[predator]),
        TRUE ~ shorten_pred_label(predator)
      ),
      # Replace underscores with spaces so LaTeX never needs escaping
      Predator_Clean = gsub("_", " ", Predator_Clean),
      mod_prey_clean = gsub("_", " ", mod_prey),
      
      # Compact mechanism labels
      mod_variant = gsub("Joint Temp \\+ AltPrey", "Joint Temp+Prey", mod_variant),
      mod_variant = gsub("AltPrey Buffered", "AltPrey Only", mod_variant),
      mod_variant = gsub("Temp Overlap Only", "Temp Only", mod_variant),
      mod_variant = gsub("Nominal Biomass \\(Base\\)", "Base (Nominal)", mod_variant),
      
      base_str = paste0(sprintf("%.2f", base_beta), " (", sprintf("%.2f", base_p), ")"),
      mod_str  = paste0(sprintf("%.2f", mod_beta), " (", sprintf("%.2f", mod_p), ")"),
      daic_str = sprintf("%.2f", delta_aic),
      
      # Green = significant negative top-down or delta_aic < 0; Red = positive & p < 0.20
      Base_Display = case_when(
        (raw_mod_variant == "Nominal Biomass (Base)" | mod_variant == "Base (Nominal)") &
          base_p <= 0.05 & base_beta < 0 ~ paste0("\\textcolor[HTML]{1B5E20}{\\textbf{", base_str, "}}"),
        base_beta > 0 & base_p < 0.20   ~ paste0("\\textcolor[HTML]{B71C1C}{\\textbf{", base_str, "}}"),
        TRUE                            ~ base_str
      ),
      Mod_Display = if_else(
        (raw_mod_variant != "Nominal Biomass (Base)" & mod_variant != "Base (Nominal)") &
          mod_p <= 0.05 & mod_beta < 0,
        paste0("\\textcolor[HTML]{1B5E20}{\\textbf{", mod_str, "}}"),
        mod_str
      ),
      Delta_Display = if_else(
        delta_aic < 0,
        paste0("\\textcolor[HTML]{1B5E20}{\\textbf{", daic_str, "}}"),
        daic_str
      ),
      Base_AIC_Display = sprintf("%.2f", base_aic),
      Mod_AIC_Display  = sprintf("%.2f", mod_aic)
    ) %>%
    arrange(mod_aic) %>%
    select(
      Predator_Clean, Base_Display, Base_AIC_Display,
      mod_variant, mod_prey_clean, Mod_Display, Mod_AIC_Display, Delta_Display
    ) %>%
    write_csv(out_csv_path)
}

export_peff_display_csv(res_24, file.path(data_out_dir, "obj3_table_display_24yr.csv"))
export_peff_display_csv(res_19, file.path(data_out_dir, "obj3_table_display_19yr.csv"))
