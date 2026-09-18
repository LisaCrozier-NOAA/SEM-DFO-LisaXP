#plot sst top model marginal effects


library(tidyverse)
library(viridis)
library(patchwork)

# ------------------------------------------------------------------------------
# STEP 1: DEFINE PARAMETERS FOR PAIR ID 95
# ------------------------------------------------------------------------------
# NCC Parameters
b_p_ncc   <- 0.319   # X10_Californian_s_l_smoltYear_WS
b_t_ncc   <- -0.136  # X01_habCompInd
b_int_ncc <- -0.153  # ncc_sst_int

# AK Parameters
b_cpue    <- -2.633  # CPUE -> SAR Stage Link
b_p_ak    <- 1.280   # X15_PacificCodBiomass_predAK_smoltYear
b_t_ak    <- -1.003  # X21_sst_egoa_junjulaug
b_int_ak  <- 6.011   # ak_sst_int

# ------------------------------------------------------------------------------
# STEP 2: CREATE CONTINUOUS GRID ACROSS OBSERVED DOMAIN (-3 TO +3 SD)
# ------------------------------------------------------------------------------
grid_seq <- seq(-2.5, 2.5, length.out = 200)

# Grid 1: NCC Surface
ncc_grid <- expand_grid(
  P_ncc = grid_seq,
  T_ncc = grid_seq
) %>%
  mutate(
    CPUE = b_p_ncc * P_ncc + b_t_ncc * T_ncc + b_int_ncc * (P_ncc * T_ncc)
  )

# Grid 2: AK Surface (Direct Predator & Thermal Interaction Effects on SAR)
ak_grid <- expand_grid(
  P_ak = grid_seq,
  T_ak = grid_seq
) %>%
  mutate(
    SAR_Direct = b_p_ak * P_ak + b_t_ak * T_ak + b_int_ak * (P_ak * T_ak)
  )

# ------------------------------------------------------------------------------
# STEP 3: BUILD CONTOUR MAPS
# ------------------------------------------------------------------------------
# Plot A: NCC Contour Surface
p_ncc <- ggplot(ncc_grid, aes(x = P_ncc, y = T_ncc, z = CPUE)) +
  geom_contour_filled(bins = 15) +
  geom_contour(color = "white", alpha = 0.3) +
  scale_fill_viridis_d(option = "magma", name = "CPUE") +
  labs(
    title = "A) NCC Stage: CPUE Surface",
    subtitle = "Target ~ CA Sea Lion (P) + HabCompInd (T) + (P * T)",
    x = "CA Sea Lion Abundance (X10_Californian_s_l_smoltYear_WS)",
    y = "Habitat Compression Index (X01_habCompInd)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Plot B: AK Contour Surface
p_ak <- ggplot(ak_grid, aes(x = P_ak, y = T_ak, z = SAR_Direct)) +
  geom_contour_filled(bins = 15) +
  geom_contour(color = "white", alpha = 0.3) +
  scale_fill_viridis_d(option = "viridis", name = "SAR (Direct)") +
  labs(
    title = "B) AK Stage: Direct SAR Surface",
    subtitle = "Target ~ Pacific Cod (P) + EGoA SST (T) + (P * T)",
    x = "Pacific Cod Biomass (X15_PacificCodBiomass_predAK_smoltYear)",
    y = "Summer EGoA SST (X21_sst_egoa_junjulaug)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Combine and Export Plots
combined_plot <- p_ncc + p_ak + plot_layout(ncol = 2)

ggsave("Rcode_for_paper/Routput_for_paper/figures/Goal3_Pair95_Interaction_Contours.png", 
       combined_plot, width = 12, height = 5.5, dpi = 300)

print(combined_plot)
