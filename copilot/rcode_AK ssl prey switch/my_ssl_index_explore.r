

source("function/ssl_index_fxn.r")

ssl_index_fxn<-function (ssl_scaled, sst_scaled, herring_scaled, capelin_scaled,
                         scalar_ssl = 1,
                         scalar_sst = -0.3,
                         scalar_herring = 2,
                         scalar_capelin = 1) {
  
  # Uncalibrated Fixed-Weight Hazard Sum
  # I_SSL_raw = -1 * ( 1.0 * ssl_scaled 
  #                    - 0.3 * (ssl_scaled * sst_scaled) 
  #                    + 2.0 * (ssl_scaled * herring_scaled) 
  #                    + 1.0 * (ssl_scaled * capelin_scaled) )
  
  I_SSL_raw = -1 * ( scalar_ssl * ssl_scaled 
                     + scalar_sst * (ssl_scaled * sst_scaled) 
                     + scalar_herring * (ssl_scaled * herring_scaled) 
                     + scalar_capelin * (ssl_scaled * capelin_scaled) )
  
  # Standardize final index for SEM stability
  I_SSL_simple = as.vector(scale(I_SSL_raw))
  
  return(I_SSL_raw)
}

#Test fxn-----------

# Standard deviation sequence on x-axis (-2 to +2 SD)
sd_seq <- seq(-2, 2, length.out = 100)

# 2. Build marginal datasets for each of the 4 factors
# SSL Scaled varying
df_ssl <- data.frame(
  x_sd = sd_seq,
  factor = "1. SSL Scaled",
  ssl_state = "SSL = Variable",
  prediction = ssl_index_fxn(sd_seq, 0, 0, 0)
)

# SST Scaled varying
df_sst <- expand.grid(x_sd = sd_seq, ssl_state = c("SSL = 0 (Mean)", "SSL = +1 SD", "SSL = -1 SD")) %>%
  mutate(
    factor = "2. SST Scaled",
    ssl_val = case_when(ssl_state == "SSL = 0 (Mean)" ~ 0, 
                        ssl_state == "SSL = +1 SD" ~ 1,
                        TRUE ~ -1),
    prediction = ssl_index_fxn(ssl_val, x_sd, 0, 0)
  )

# Herring Scaled varying
df_herring <- expand.grid(x_sd = sd_seq, ssl_state = c("SSL = 0 (Mean)", "SSL = +1 SD", "SSL = -1 SD")) %>%
  mutate(
    factor = "3. Herring Scaled",
    ssl_val = case_when(ssl_state == "SSL = 0 (Mean)" ~ 0, ssl_state == "SSL = +1 SD" ~ 1, TRUE ~ -1),
    prediction = ssl_index_fxn(ssl_val, 0, x_sd, 0) # herring varies
  )

# Capelin Scaled varying
df_capelin <- expand.grid(x_sd = sd_seq, ssl_state = c("SSL = 0 (Mean)", "SSL = +1 SD", "SSL = -1 SD")) %>%
  mutate(
    factor = "4. Capelin Scaled",
    ssl_val = case_when(ssl_state == "SSL = 0 (Mean)" ~ 0, ssl_state == "SSL = +1 SD" ~ 1, TRUE ~ -1),
    prediction = ssl_index_fxn(ssl_val, 0, 0, x_sd)
  )

# Combine into single long data frame
df_all <- bind_rows(df_ssl, df_sst, df_herring, df_capelin)

# 3. Create the 4-panel prediction plot
ggplot(df_all, aes(x = x_sd, y = prediction, color = ssl_state, linetype = ssl_state)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ factor, scales = "free_y") +
  scale_color_manual(values = c(
    "SSL = 0 (Mean)" = "#1f77b4",
    "SSL = +1 SD"    = "#d62728",
    "SSL = -1 SD"    = "#2ca02c",
    "SSL = Variable" = "#1f77b4"
  )) +
  scale_linetype_manual(values = c(
    "SSL = 0 (Mean)" = "solid",
    "SSL = +1 SD"    = "dashed",
    "SSL = -1 SD"    = "dotted",
    "SSL = Variable" = "solid"
  )) +
  labs(
    title = "Marginal Response Curves across Variable Standard Deviations (-2 to +2 SD)",
    x = "Factor Value (Standard Deviations)",
    y = "Predicted Raw Index Value",
    color = "Sea Lion Baseline",
    linetype = "Sea Lion Baseline"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  )

#Test correlation structure in the actual database---------
safe_scale <- function(x) {
  if (all(is.na(x))) return(x)
  as.vector(scale(x))
}


salmon_dat<-read.csv(file.path("data_Lisa/sem_master_data.csv"))%>% 
  clean_names()  %>% 
  select(year,contains("x07"),contains("x16_sar"),contains("sar_"),contains("x09_dfa_hake"))
names(salmon_dat)


sst_dat<-read.csv(file.path("data_Lisa/goa_prey_clim_raw_trends_avg.csv"),row.names = 1) %>% 
  clean_names() %>%
  select(c(1,19:32)) %>%
  mutate(sst_wgoa_coastwatch_raw=sst_wgoa_coastwatch_junjulaug,
         sst_egoa_coastwatch_raw=sst_wgoa_coastwatch_junjulaug)
names(sst_dat)

shark_dat<-read.csv(file.path("data_Lisa/shark_wide.csv"))%>% 
  clean_names() %>%
  select(year,contains("goa")) %>%
  filter(year>=1998, year<=2021) %>%
  mutate(goa_pacific_sleeper_shark=safe_scale(goa_pacific_sleeper_shark))
names(shark_dat)

ssl_dat<-read.csv(file.path("copilot/outputs_2/ssl.dat.csv"))%>% 
  clean_names() %>%
  select(!contains("_avg"))
names(ssl_dat)

goa_prey<-read.csv(file.path("data_Lisa/goa_prey_clim_raw_trends_avg.csv"),row.names = 1) %>% 
  clean_names() %>%
  select(1:12)
round(apply(goa_prey[,-1],2,mean,na.rm=T),2)
round(apply(goa_prey[,-1],2,sd,na.rm=T),2)


data_base <- ssl_dat  %>%
  left_join(goa_prey, by = "year") %>%
  left_join(sst_dat %>% select(year,sst_wgoa_coastwatch_raw), by = "year") %>%
  left_join(shark_dat %>% select(year,goa_pacific_sleeper_shark,goa_salmon_shark), by = "year") 
names(data_base)

data_base <- data_base %>%
  mutate(i_ssl = ssl_index_fxn(ssl_scaled=ssl_seak_pup_pred,
                               sst_scaled=sst_wgoa_coastwatch_junjulaug ,
                               herring_scaled =  egoa_bio_stka_herr_matbiom,
                               capelin_scaled = wgoa_bio_mid_il_capelin)) %>%
  mutate(i_shark = shark_index_fxn.v2(shark_scaled=goa_pacific_sleeper_shark,
                                      temp_raw_Mt=sst_wgoa_coastwatch_raw,
                                      temp_ref_Mt=mean(sst_wgoa_coastwatch_raw,na.rm=T),
                                      temp_Ot = enso_dj
                                      )
         )
