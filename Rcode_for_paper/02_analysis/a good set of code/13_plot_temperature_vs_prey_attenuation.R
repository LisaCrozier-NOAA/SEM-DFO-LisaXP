

# ==============================================================================
# Script: 13_plot_temperature_vs_prey_attenuation.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Load saved sweep result CSVs, reconstruct $P_{\text{eff}}$, and plot
#          2-panel figures (Panel A: Temp/Spatial, Panel B: Prey 1 & 2) 
#          ordering Spatial/Temp models first.
# Output: Rcode_for_paper/Routput_for_paper/plots/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------------------------
# STEP 0: PATHS & INPUT FILE SELECTION
# ------------------------------------------------------------------------------
years_included<-"19yr"
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")
plot_dir     <- file.path(master_out, "plots",years_included)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Select the Dataset & Matching Sweep CSV
# --- TOGGLE SWITCH HERE ---
#years_included<-"24yr"
dataset_file <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
sweep_file   <- file.path(data_out_dir, paste0("Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_",years_included,".csv"))

# (For 19-year run, uncomment below:)
# dataset_file <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv")
# sweep_file   <- file.path(data_out_dir, "Goal3_5Guild_Joint_Sweep_Results_minmax_combo_size_2_19yr_noNA_dougsmoothed.csv")

if (!file.exists(dataset_file)) stop("Missing dataset file: ", dataset_file)
if (!file.exists(sweep_file))   stop("Missing sweep CSV file: ", sweep_file)

blacklist <- c("canaryRockfish", "Harbor_seal_CR",
               "X09_DFA_ChinAbundSnakeFall_adultyr","X09_DFA_HakeAge5Plus_adultyr",
               "sharkCatchGoA","X15_spinyDogfishGoA_predAK_adultyr")
blacklist_pattern <- paste(blacklist, collapse = "|")

raw_sem_data <- read.csv(dataset_file, stringsAsFactors = FALSE)
sweep_df     <- read.csv(sweep_file, stringsAsFactors = FALSE) %>%
  filter(abs(delta_aic)>6) %>%
  filter(!grepl(blacklist_pattern, predator))



cat(sprintf("Loaded dataset with %d rows and sweep results with %d models.\n", nrow(raw_sem_data), nrow(sweep_df)))

# ------------------------------------------------------------------------------
# STEP 1: HELPER FOR MIN-MAX SCALING & EFFECTIVE PREDATOR RECONSTRUCTION
# ------------------------------------------------------------------------------
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# Scale all numeric indicators into [0, 1]
minmax_data <- raw_sem_data %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), min_max_scale))

# Identify spatial/temperature columns based on naming pattern
is_spatial_col <- function(col_name) {
  grepl("^(X01_habCompInd|X12_DFA_biomassEuph|X21_sst_egoa)", col_name)
}
which(grepl("^(X08_DFA_DC_corm_3_WS)", names(raw_sem_data)))
which(grepl("^(X08_DFA_DC_corm_3_WS)", guild_pred_lookup$pred_data_col))

# ------------------------------------------------------------------------------
# STEP 2: ORDER MODELS (SPATIAL/TEMP ATTENUATION FIRST)
# ------------------------------------------------------------------------------
sweep_df <- sweep_df %>%
  mutate(
    # Identify if selected_prey contains a spatial/temp indicator
    used_spatial = sapply(selected_prey, function(sp) {
      if (is.na(sp) || sp == "None (Unattenuated)") return(FALSE)
      gates <- unlist(strsplit(sp, " \\+ "))
      any(sapply(gates, is_spatial_col))
    })
  ) %>%
  # Sort: Spatial/Temp models FIRST, then by best Delta AIC
  arrange(desc(used_spatial), delta_aic) %>%
  as_tibble() %>%
  distinct()

cat("\n==============================================================================\n")
cat("      MODELS LOADED & ORDERED (SPATIAL/TEMP ATTENUATION FIRST)               \n")
cat("==============================================================================\n\n")
print(sweep_df %>% select(predator, selected_prey, fitted_weights, used_spatial, beta_peff, delta_aic), n = 50)

# ------------------------------------------------------------------------------
# STEP 3: 2-PANEL PLOTTING FUNCTION
# ------------------------------------------------------------------------------
plot_saved_model_row <- function(model_row, df_model) {
  pred_col    <- model_row$predator
  gates_str   <- model_row$selected_prey
  weights_str <- model_row$fitted_weights
  
  if (!pred_col %in% names(df_model)) return(NULL)
  
  Years  <- df_model$Year
  P_norm <- df_model[[pred_col]]
  
  # Parse gates and weights
  has_gates <- !is.na(gates_str) && gates_str != "None (Unattenuated)" && gates_str != "None"
  
  if (has_gates) {
    gates   <- unlist(strsplit(gates_str, " \\+ "))
    weights <- as.numeric(unlist(strsplit(weights_str, "; ")))
  } else {
    gates   <- c()
    weights <- c()
  }
  
  # Identify Spatial/Temp Gate vs. Prey Gates
  sp_gate_idx   <- if (has_gates) which(sapply(gates, is_spatial_col)) else integer(0)
  prey_gate_idx <- if (has_gates) setdiff(seq_along(gates), sp_gate_idx) else integer(0)
  
  has_sp   <- length(sp_gate_idx) > 0
  has_prey <- length(prey_gate_idx) > 0
  
  # Calculate Reconstructed P_eff
  P_eff <- P_norm
  if (has_gates) {
    for (g in seq_along(gates)) {
      gate_var <- gates[g]
      w_val    <- weights[g]
      if (gate_var %in% names(df_model)) {
        A_norm <- df_model[[gate_var]]
        P_eff  <- P_eff * (1 - w_val * A_norm)
      }
    }
  }
  
  # Calculate Partial Spatial P_eff (for Left Panel)
  P_sp_eff <- P_norm
  if (has_sp) {
    sp_var <- gates[sp_gate_idx[1]]
    sp_w   <- weights[sp_gate_idx[1]]
    if (sp_var %in% names(df_model)) {
      P_sp_eff <- P_norm * (1 - sp_w * df_model[[sp_var]])
    }
  }
  
  # Layout Setup: 1 row, 2 columns
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0, 0, 3, 0))
  
  # ----------------------------------------------------------------------------
  # PANEL A: TEMPERATURE / SPATIAL OVERLAP ATTENUATION
  # ----------------------------------------------------------------------------
  if (has_sp) {
    sp_var <- gates[sp_gate_idx[1]]
    sp_w   <- weights[sp_gate_idx[1]]
    S_norm <- df_model[[sp_var]]
    
    clean_sp_name <- gsub("^X[0-9]{2}_", "", sp_var)
    
    plot(Years, P_norm, type = "b", pch = 16, col = "gray40", ylim = c(0, 1.15),
         main = "A. Spatial Overlap / Temp Gate", xlab = "Year", ylab = "Index Scale [0, 1]")
    lines(Years, P_sp_eff, type = "b", pch = 17, col = "firebrick3", lty = 2, lwd = 2)
    lines(Years, S_norm, type = "l", col = "pink", lty = 3, lwd = 3)
    
    legend("topleft", 
           legend = c("Raw Predator (P_norm)", 
                      paste0("Temp-Gated (w=", sp_w, ")"), 
                      paste0("Spatial Index: ", clean_sp_name)),
           col = c("gray40", "firebrick3", "pink"), 
           pch = c(16, 17, NA), lty = c(1, 2, 3), lwd = c(1, 2, 3), bty = "n", cex = 0.75)
  } else {
    plot(Years, P_norm, type = "b", pch = 16, col = "gray50", ylim = c(0, 1.15),
         main = "A. Spatial Overlap / Temp Gate (Not Selected)", xlab = "Year", ylab = "Index Scale [0, 1]")
    lines(Years, P_eff, type = "b", pch = 17, col = "firebrick3", lty = 2, lwd = 1.5)
    legend("topleft", 
           legend = c("Raw Predator (P_norm)", "Final Joint P_eff"), 
           col = c("gray50", "firebrick3"), pch = c(16, 17), lty = c(1, 2), bty = "n", cex = 0.8)
  }
  
  # ----------------------------------------------------------------------------
  # PANEL B: ALTERNATE PREY ATTENUATION (PREY 1 & PREY 2)
  # ----------------------------------------------------------------------------
  plot(Years, P_norm, type = "b", pch = 16, col = "gray40", ylim = c(0, 1.15),
       main = "B. Alternate Prey Attenuation", xlab = "Year", ylab = "Index Scale [0, 1]")
  lines(Years, P_eff, type = "b", pch = 17, col = "darkred", lty = 2, lwd = 2.5)
  
  leg_txt <- c("Raw Predator (P_norm)", "Final Joint Effective (P_eff)")
  leg_col <- c("gray40", "darkred")
  leg_pch <- c(16, 17)
  leg_lty <- c(1, 2)
  
  prey_colors <- c("steelblue3", "darkorange2")
  
  if (has_prey) {
    for (p_i in seq_along(prey_gate_idx)) {
      idx    <- prey_gate_idx[p_i]
      p_name <- gates[idx]
      
      if (p_name %in% names(df_model)) {
        A_norm  <- df_model[[p_name]]
        p_color <- prey_colors[p_i]
        
        lines(Years, A_norm, type = "l", col = p_color, lty = 3, lwd = 1.8)
        
        clean_p_name <- gsub("^X[0-9]{2}_", "", p_name)
        leg_txt <- c(leg_txt, paste0("Prey ", p_i, ": ", clean_p_name))
        leg_col <- c(leg_col, p_color)
        leg_pch <- c(leg_pch, NA)
        leg_lty <- c(leg_lty, 3)
      }
    }
  }
  
  legend("topright", legend = leg_txt, col = leg_col, pch = leg_pch, lty = leg_lty, 
         lwd = c(1, 2.5, rep(1.8, length(prey_gate_idx))), bty = "n", cex = 0.72)
  
  # Overall Title
  title(main = paste0("Predator: ", pred_col, "\nSelected Gates: ", gates_str, 
                      " | ΔAIC = ", round(model_row$delta_aic, 2), " | R² = ", round(model_row$final_r2, 3)), 
        outer = TRUE, cex.main = 1.1)
}

# ------------------------------------------------------------------------------
# STEP 4: GENERATE AND SAVE PLOTS
# ------------------------------------------------------------------------------
cat("\nGenerating comparative 2-panel figures for all sweep models...\n")

for (i in 1:nrow(sweep_df)) {
  m_row <- sweep_df[i, ]
  pred  <- m_row$predator
  
  plot_out_file <- file.path(plot_dir, paste0("Temp_vs_Prey_Plot_", pred,"_",years_included, ".png"))
  
  png(filename = plot_out_file, width = 10, height = 5, units = "in", res = 300)
  plot_saved_model_row(m_row, minmax_data)
  dev.off()
}

cat("\nPlotting complete! All figures saved to: ", plot_dir, "\n")
