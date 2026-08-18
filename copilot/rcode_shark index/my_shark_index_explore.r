#load data---------
salmon_dat<-read.csv(file.path("data_Lisa/sem_master_data.csv"))%>% 
  clean_names()  %>% 
  select(year,contains("x07"),contains("x16_sar"),contains("sar_"),contains("x09_dfa_hake"))
names(salmon_dat)


shark_dat<-read.csv(file.path("data_Lisa/shark_wide.csv"))%>% 
  clean_names() %>%
  select(year,goa_pacific_sleeper_shark,goa_salmon_shark)  
head(shark_dat)


sst<-read.csv(file.path("data_Lisa/ak_yr.csv"),row.names=NULL) %>%
  select(year,contains("sst"),contains("176to226m"),contains("pdo"),contains("enso"),contains("195to205m"),contains("246to255")) %>%
  filter(year>1997);names(sst)

#too many NAs in the other deep water temps, swln_temp_spr_176to226m had no NAs 

summary(sst)

sst0<-sst %>%
  select(year,contains("coastwatch"),contains("swln_temp_spr_176to226m"),contains("pdo"),contains("enso"))
summary(sst0)

ssl_dat<-read.csv(file.path("copilot/outputs_2/ssl.dat.csv"))%>% 
  clean_names() %>%
  select(!contains("_avg"))
names(ssl_dat)

goa_prey<-read.csv(file.path("data_Lisa/goa_prey_clim_raw_trends_avg.csv"),row.names = 1) %>% 
  clean_names() %>%
  select(1:12)
names(goa_prey)

names(shark_dat)
names(sst0)

pred_index_dat <- left_join(sst0 %>% select(year,sst_wgoa_coastwatch_junjulaug,pdo_djf),
                            shark_dat %>% select(year,goa_pacific_sleeper_shark),by="year") %>%
  left_join(ssl_dat %>% select(year,ssl_seak_pup_pred),by="year") %>%
  left_join(goa_prey %>% select(year,egoa_bio_stka_herr_matbiom,wgoa_bio_mid_il_capelin),by="year") %>%
  left_join(salmon_dat, by="year")
head(pred_index_dat)

#Raw data plots------------
library(tidyverse)

# 1. Join shark data with sst0 dataset
shark_base <- shark_dat %>%
  inner_join(sst0, by = "year")

# 2. Pivot long & classify categories
shark_inputs_long <- shark_base %>%
  pivot_longer(-year, names_to = "variable", values_to = "value") %>%
  mutate(
    category = case_when(
      variable %in% c("sst_egoa_coastwatch_junjulaug", "sst_wgoa_coastwatch_junjulaug", "swln_temp_spr_176to226m") ~ "1. Raw Temperature for Mt (°C)",
      variable %in% c("pdo_djf", "enso_dj") ~ "2. Synoptic Climate for Ot (Index)",
      variable %in% c("goa_pacific_sleeper_shark", "goa_salmon_shark") ~ "3. Shark Abundance Indices (Scaled)",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(variable) %>%
  mutate(
    # Scale ONLY the shark abundance variables for comparative plotting
    plot_value = if_else(category == "3. Shark Abundance Indices (Scaled)", as.vector(scale(value)), value)
  ) %>%
  ungroup()

# 3. Faceted Time-Series Plot
shark_data_plot<-
ggplot(shark_inputs_long, aes(x = year, y = plot_value, color = variable, group = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 1.5, alpha = 0.85) +
  facet_wrap(~ category, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(shark_base$year, na.rm = TRUE), max(shark_base$year, na.rm = TRUE), by = 2)) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Shark Index Input Series Dynamics by Role",
    subtitle = "Raw Temps (°C) for Mt | Unscaled PDO/ENSO for Ot | Scaled Shark Abundance",
    x = "Year",
    y = "Value",
    color = "Variable Name"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(ncol = 3))


ggsave("shark_data_plot.png",shark_data_plot)


# Phase 2: Permutation Engine ----------
#for $I_{\text{Shark}}$We will build every valid combination across:2 Shark Species: goa_pacific_sleeper_shark, goa_salmon_shark3 $M_t$ Raw Temperatures: sst_egoa_coastwatch_junjulaug, sst_wgoa_coastwatch_junjulaug, swln_temp_spr_176to226m2 $O_t$ Climate Overlap Indices: pdo_djf, enso_dj (unscaled)Total Combinations = $2 \times 3 \times 2 = 12$ permutations.Important Scaling Note: The function expects shark_scaled. Inside the grid generation, each shark species counts vector is $z$-score scaled before being passed into shark_index_fxn.v2().

# 1. Define Candidate Lists
shark_vars <- c("goa_pacific_sleeper_shark", "goa_salmon_shark")
mt_temp_vars <- c("sst_egoa_coastwatch_junjulaug", "sst_wgoa_coastwatch_junjulaug", "swln_temp_spr_176to226m")
ot_temp_vars <- c("pdo_djf", "enso_dj")

# 2. Build Permutation Grid (12 combinations)
shark_perm_grid <- expand.grid(
  shark_var = shark_vars,
  mt_var    = mt_temp_vars,
  ot_var    = ot_temp_vars,
  stringsAsFactors = FALSE
)

# 3. Compute I_Shark predictions matrix across all years
shark_perm_matrix <- map_dfc(1:nrow(shark_perm_grid), function(i) {
  row <- shark_perm_grid[i, ]
  
  # Scale shark input for the function
  s_scaled <- as.vector(scale(shark_base[[row$shark_var]]))
  t_raw    <- shark_base[[row$mt_var]]
  t_ref    <- mean(t_raw, na.rm = TRUE) # Baseline average temperature for Mt
  t_ot     <- shark_base[[row$ot_var]]  # Unscaled PDO/ENSO
  
  pred <- shark_index_fxn.v2(
    shark_scaled = s_scaled,
    temp_raw_Mt  = t_raw,
    temp_ref_Mt  = t_ref,
    temp_Ot      = t_ot,
    Q10          = 2,
    overlap_form = "logistic",
    overlap_slope = 1
  )
  
  tibble(!!paste0("shark_perm_", i) := pred)
})

# 4. Attach metadata into a long data frame for plotting
shark_perm_long <- shark_perm_grid %>%
  mutate(perm_id = paste0("shark_perm_", row_number())) %>%
  inner_join(
    shark_perm_matrix %>% 
      mutate(year = shark_base$year) %>% 
      pivot_longer(cols = starts_with("shark_perm_"), names_to = "perm_id", values_to = "i_shark_pred"),
    by = "perm_id"
  )


#Phase 3: Visualizing I_Shark Predictions Across Time & Drivers------
#Graph 3A: Chronological Timeline Dodged by Overlap Index ($O_t$)
#This plot shows the trajectory of $I_{\text{Shark}}$ predictions across time, 
#with points grouped by Shark Species (Color), $M_t$ Temperature Driver (Faceted Rows), and $O_t$ Climate Overlap Choice (Dodged Bands).

shark_index_drivers_plot<-
ggplot(shark_perm_long, aes(x = factor(year), y = i_shark_pred, color = shark_var, shape = ot_var)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  # Dodged points so Ot options don't cover each other
  geom_point(
    aes(group = interaction(shark_var, ot_var)),
    position = position_dodge(width = 0.6),
    size = 2.5,
    alpha = 0.9
  ) +
  facet_wrap(~ mt_var, ncol = 1, scales = "free_y") +
  scale_color_manual(
    values = c("goa_pacific_sleeper_shark" = "#2b5c8f", "goa_salmon_shark" = "#d95f02"),
    name = "Shark Species"
  ) +
  scale_shape_manual(
    values = c("pdo_djf" = 16, "enso_dj" = 17),
    name = "Ot Climate Driver"
  ) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Predicted I_Shark Index Dynamics Across Thermal & Overlap Formulations",
    subtitle = "Faceted by Mt Temperature Depth/Location | Color: Shark Species | Shape: Ot Driver",
    x = "Year",
    y = "Predicted Raw I_Shark Index"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 10)
  )

ggsave("shark_index_drivers.png",shark_index_drivers_plot)

#Graph 3B: Sensitivity of $I_{\text{Shark}}$ to Deep Water vs. Surface Water $M_t$
#To directly prove whether substituting the Deep Water Spring Temperature (swln_temp_spr_176to226m) for Surface Temperature changes the metabolic multiplier enough to matter, 
#this scatter plot contrasts Surface $I_{\text{Shark}}$ predictions directly against Deep $I_{\text{Shark}}$ predictions for matching years:

# Pivot deep vs surface predictions side-by-side
deep_vs_surface <- shark_perm_long %>%
  mutate(temp_type = if_else(mt_var == "swln_temp_spr_176to226m", "Deep_Water", "Surface_Water")) %>%
  group_by(year, shark_var, ot_var, temp_type) %>%
  summarise(mean_i_shark = mean(i_shark_pred, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = temp_type, values_from = mean_i_shark)

MT_temp_depth<-
ggplot(deep_vs_surface, aes(x = Surface_Water, y = Deep_Water, color = shark_var, shape = ot_var)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 3, alpha = 0.85) +
  geom_text(aes(label = year), vjust = -0.7, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("goa_pacific_sleeper_shark" = "#2b5c8f", "goa_salmon_shark" = "#d95f02")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Metabolic Driver Comparison: Surface Water vs. Deep Water (176–226m)",
    subtitle = "1:1 dashed line indicates identical index predictions regardless of depth",
    x = "I_Shark Prediction using Surface Water Temp (Mt)",
    y = "I_Shark Prediction using Deep Water Temp (Mt)",
    color = "Shark Species",
    shape = "Ot Driver"
  ) +
  theme(legend.position = "bottom")


ggsave("MT_temp_depth.png",MT_temp_depth)


#Compare ssl and shark indices--------------
library(tidyverse)

# 1. Build your final clean baseline dataset for SEM
sem_input_data <- pred_index_dat %>%
  mutate(
    # Selected SSL Index (Using Southeast Alaska Pup counts & standard climate drivers)
    I_SSL_final = ssl_index_fxn(
      ssl_scaled     = as.vector(scale(ssl_seak_pup_pred)),
      sst_scaled     = as.vector(scale(sst_wgoa_coastwatch_junjulaug)),
      herring_scaled = as.vector(scale(egoa_bio_stka_herr_matbiom)),
      capelin_scaled = as.vector(scale(wgoa_bio_mid_il_capelin))
    ),
    
    # Selected Shark Index (Pacific Sleeper Shark + Surface Temp + PDO)
    I_Shark_final = shark_index_fxn.v2(
      shark_scaled  = as.vector(scale(goa_pacific_sleeper_shark)),
      temp_raw_Mt   = sst_wgoa_coastwatch_junjulaug,
      temp_ref_Mt   = mean(sst_wgoa_coastwatch_junjulaug, na.rm = TRUE),
      temp_Ot       = pdo_djf,
      Q10           = 2,
      overlap_form  = "logistic",
      overlap_slope = 1
    )
  ) %>%
  # Scale final outputs for SEM stability
  mutate(
    I_SSL_sem   = as.vector(scale(I_SSL_final)),
    I_Shark_sem = as.vector(scale(I_Shark_final))
  )

head(sem_input_data)


# 2. Plot the final two SEM indicator trajectories together
ggplot(sem_input_data, aes(x = year)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(y = I_SSL_sem, color = "I_SSL (Sea Lion Index)"), linewidth = 1.1) +
  geom_line(aes(y = I_Shark_sem, color = "I_Shark (Sleeper Shark Index)"), linewidth = 1.1) +
  geom_point(aes(y = I_SSL_sem, color = "I_SSL (Sea Lion Index)"), size = 2) +
  geom_point(aes(y = I_Shark_sem, color = "I_Shark (Sleeper Shark Index)"), size = 2) +
  scale_color_manual(values = c("I_SSL (Sea Lion Index)" = "#2b5c8f", "I_Shark (Sleeper Shark Index)" = "#d95f02")) +
  scale_x_continuous(breaks = seq(min(sem_input_data$year), max(sem_input_data$year), by = 2)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Final Standardized SEM Index Inputs (1998–2021)",
    subtitle = "Comparing I_SSL and I_Shark trajectories ready for structural equation modeling",
    x = "Year",
    y = "Standardized Value (z-score)",
    color = "Index"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold")
  )
