#Shift focus to alternate prey for CR ssl
output_dir <- "copilot/outputs_csl_cr"

#pre-2011 Eulachon data
path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Data/New folder (3)/Eulachon_data_Gustavson_2010_status_review.xlsx"
library(readxl)
eulachon_landings<-read_xlsx(path=path,sheet=1,skip=1) %>%  select(Year,Total)%>%rename(eulachon_cr_landings_pounds=Total)
head(eulachon_landings)

      #same data but includes more years, don't need at the moment:
      # eulachon_count<-read_xlsx(path=path,sheet=2,skip=1) %>%  filter(Year>=1998) %>% rename(Cr_trib_pounds=`Total landings\r\n(pounds)`)
      # head(eulachon_count)
      # 
      # eulachon_dps<-read_xlsx(path=path,sheet=3,skip=1) %>%  filter(Year>=1998) %>%  select(1:2) %>%rename(Columbia_River_mt=`Columbia\r\nRiver`)
      # head(eulachon_dps)

plot()
#Other spp
library(readxl)


csl_weekly_1998_2024<-read.csv(file.path(output_dir, "csl_reconstructed_weekly_1998_2024.csv"), row.names = NULL)
eulachon_weekly_1998_2024<- read.csv(file.path(output_dir, "eulachon_reconstructed_weekly_1998_2024.csv"), row.names = NULL)

bonn_daily<- read_excel("data_Lisa/Adult_BONpassage_1976_2024.xlsx",sheet=1);head(bonn_daily)
shad_daily<- read.csv("data_Lisa/SHAD_bonn_19382026.csv",row.names = NULL);head(shad_daily)


# Load jake's data-----------------------------------------------------------------------------
jake_path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Jake Marshall/"
load(paste0(jake_path,"lre_dat_yearly.RData"),verbose=T)
head(lre_dat_yearly)
names(lre_dat_yearly)
write.csv(lre_dat_yearly,"data_Lisa/lre_dat_yearly.csv",row.names=FALSE)


week_all<-read.csv(paste0(jake_path,"Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.all.results.csv"),row.names=NULL)
head(week_all)

week_main<-read.csv(paste0(jake_path,"Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.main.results.csv"),row.names=NULL)
head(week_main)

year_main<-read.csv(paste0(jake_path,"Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.annual.main.results.csv"),row.names=NULL)
head(year_main)

year_main 

#find missing data in the week_main file---------
library(tidyverse)

# 1. Read in weekly data
week_main <- read.csv(
  paste0(jake_path, "Survival Meta analysis final 10012025/Jake Marshall -- Final Product/Jake.weekly.main.results.csv"),
  row.names = NULL
)

# Quick check of structure
head(week_main)

# 2. Filter for Spring weeks in 2011-2024 and isolate missing csl_nonpup_total_emb
spring_csl_missing <- week_main %>%
  # Filter for season and target years
  filter(
    month %in% c(4, 5, 6),
    #tolower(season) == "spring",
    year >= 2011, year <= 2024
  ) %>%
  # Select key metadata + target sealion count column
  select(year, week, month, date, season, csl_nonpup_total_emb) %>%
  # Identify rows where csl_nonpup_total_emb is NA
  filter(is.na(csl_nonpup_total_emb))

spring_csl_missing

#SSL -- only 4 surveys total at south jetty
dat<-week_all  %>% select(year,week,month,date,season,csl_nonpup_total_emb,csl_nobs_emb,ssl_nonpup_total_sj,ssl_nobs_sj)
dat %>% filter(!is.na(ssl_nonpup_total_sj)) %>% select(year,week,month,date,season,ssl_nonpup_total_sj,ssl_nobs_sj)
dat %>% filter(!is.na(ssl_nobs_sj)) %>% select(year,week,month,date,season,ssl_nonpup_total_sj,ssl_nobs_sj)

dat<-lre_dat_yearly  %>% select(year,csl_nonpup_total_emb,csl_nobs_emb,ssl_nonpup_total_sj,ssl_nobs_sj)
dat1<-dat %>% filter(!is.na(ssl_nonpup_total_sj)) %>% select(year,ssl_nonpup_total_sj,ssl_nobs_sj) %>% filter(year>=1997)
print(dat1,n=Inf)
# year ssl_nonpup_total_sj ssl_nobs_sj
# 1  1997                 135           4
# 2  1998                 301           2
# 3  1999                 371           6
# 4  2000                 284           6
# 5  2001                 590           2
# 6  2002                 281           3
# 7  2003                 707           1
# 8  2006                 705           1
# 9  2008                 386           1
# 10  2013                 399           1
# 11  2015                 421           1
# 12  2017                 358           1
# 13  2021                 217           2
# -----------------------------------------------------------------------------
# 3. Diagnostic Summaries
# -----------------------------------------------------------------------------

# A. List of exact missing year-week combinations
cat("--- MISSING SPRING WEEKS FOR CSL_NONPUP_TOTAL_EMB (2011-2024) ---/n")
cat("Total missing spring week-year entries:", nrow(spring_csl_missing), "/n/n")
print(spring_csl_missing, n = 50)

# B. Summary count of missing weeks by year
missing_by_year <- week_main %>%
  filter(tolower(season) == "spring", month %in% c(4,5,6), year >= 2011, year < 2024) %>%
  group_by(year) %>%
  summarise(
    total_spring_weeks = n(),
    missing_csl_weeks  = sum(is.na(csl_nonpup_total_emb)),
    available_csl_weeks = sum(!is.na(csl_nonpup_total_emb))
  )

cat("/n--- SPRING CSL COVERAGE SUMMARY BY YEAR ---/n")
print(missing_by_year, n = 20)

# C. Grid showing missing (NA) vs available (OK) spring weeks
spring_week_grid <- week_main %>%
  filter(tolower(season) == "spring", year >= 2011, year < 2024) %>%
  mutate(status = if_else(is.na(csl_nonpup_total_emb), "MISSING", "OK")) %>%
  pivot_wider(
    id_cols = week,
    names_from = year,
    values_from = status,
    values_fill = "NO_RECORD"
  ) %>%
  arrange(week)

cat("/n--- WEEK x YEAR AVAILABILITY GRID ---/n")
print(spring_week_grid, n = 30)

#linearly 

#process yearly abundance data------
lre_dat<-lre_dat_yearly %>% select(year,Shad_bonn_pass,eulachon_ssb_est,hake_cpue, 
                                   hs_nonpup_total_ds,csl_nonpup_total_emb,ssl_nonpup_total_sj,total_bonn) %>%
  filter(year>=1998,year<=2021)

head(lre_dat)
#the real problem is lack of eulachon data from 1998-2011
#let's interpolate the sea lion data, but drop hs because of 2015-2021 lack of data
#ssl at sk is also missing 4 years in a row 2009-2012. but assume constant?
#what are we going to do about eulachon? set it at its min?

lre_dat<-lre_dat %>% select(-hs_nonpup_total_ds) 

library(tidyverse)
library(zoo)

# Data prep pipeline for lre_dat
lre_dat_clean <- lre_dat %>%
  # 1. Ensure dataset is sorted chronologically
  arrange(year) %>%
  
  # 2. Set NAs in the first 4 years of total_bonn to 0
  mutate(
    total_bonn = if_else(row_number() <= 4 & is.na(total_bonn), 0, total_bonn)
  ) %>%
  
  # 3. Linearly interpolate missing data in all columns EXCEPT eulachon_ssb_est
  mutate(
    across(
      .cols = -c(year, eulachon_ssb_est),
      .fns  = ~ zoo::na.approx(., na.rm = FALSE)
    )
  )

# Verify summary of missing values
cat("--- Missing Value Summary After Prep ---/n")
lre_dat_clean %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "na_count") %>%
  print(n = 30)

dat<-lre_dat_clean


#fit spline to eulachon-------
library(tidyverse)
library(mgcv) # For GAM smooth splines

# 1. Fit GAM Smooth Spline to observed Eulachon data
# Using a penalized thin-plate regression spline (s())
eulachon_gam <- gam(
  eulachon_ssb_est ~ s(year, k = 8), 
  data = lre_dat_clean,
  method = "REML"
)

# 2. Generate smooth predictions across all years (including pre-2011)
eulachon_pred <- lre_dat_clean %>%
  select(year, eulachon_ssb_est) %>%
  mutate(
    # Fit prediction line
    fit_ssb = as.vector(predict(eulachon_gam, newdata = ., type = "response")),
    # Get standard errors for 95% confidence intervals
    se_ssb  = predict(eulachon_gam, newdata = ., type = "se.fit"),
    lwr_95  = fit_ssb - 1.96 * se_ssb,
    upr_95  = fit_ssb + 1.96 * se_ssb,
    # Enforce non-negativity constraint (SSB cannot be < 0)
    lwr_95  = pmax(0, lwr_95),
    fit_ssb = pmax(0, fit_ssb)
  )

# 3. Plot Observed vs. Spline Reconstructed Eulachon Trajectory
ggplot(eulachon_pred, aes(x = year)) +
  # 95% Confidence Interval Ribbon
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), fill = "steelblue", alpha = 0.2) +
  # Spline Fitted Line
  geom_line(aes(y = fit_ssb, color = "Spline Fit"), size = 1.2) +
  # Actual Observed Points
  geom_point(aes(y = eulachon_ssb_est, color = "Observed Data"), size = 3) +
  # Highlight the 2011 methodology transition
  geom_vline(xintercept = 2011, linetype = "dashed", color = "firebrick", size = 0.8) +
  annotate("text", x = 2011.2, y = max(eulachon_pred$eulachon_ssb_est, na.rm = TRUE) * 0.9, 
           label = "2011: Modern Survey Era", hjust = 0, color = "firebrick", fontface = "bold") +
  scale_color_manual(values = c("Observed Data" = "black", "Spline Fit" = "steelblue")) +
  labs(
    title = "Eulachon Spawning Stock Biomass (SSB) Spline Reconstruction",
    subtitle = "Smoothing spline fit with 95% CIs pre- and post-2011 methodology shift",
    x = "Year",
    y = "Eulachon SSB Estimate",
    color = "Series"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")


#fit a sine curve to eulachon-------
library(tidyverse)

# 1. Filter observed non-NA data
eul_obs <- lre_dat_clean %>% 
  filter(!is.na(eulachon_ssb_est))

# 2. Get initial parameter guesses
mean_y <- mean(eul_obs$eulachon_ssb_est)
amp_y  <- (max(eul_obs$eulachon_ssb_est) - min(eul_obs$eulachon_ssb_est)) / 2

# 3. Fit Sine Wave Model using Non-linear Least Squares (nls)
# We test a baseline period guess of ~8-10 years (typical forage fish cycle)
sine_fit <- nls(
  eulachon_ssb_est ~ A * sin((2 * pi / P) * year + phi) + C,
  data = eul_obs,
  start = list(
    A   = amp_y,
    P   = 9,          # Starting period guess in years
    phi = 0,
    C   = mean_y
  ),
  control = nls.control(maxiter = 500, warnOnly = TRUE)
)

# Display fitted model parameters (Period, Amplitude, Phase)
summary(sine_fit)

# 4. Predict across the ENTIRE time series (extrapolating back to start)
lre_dat_eulachon_sine <- lre_dat_clean %>%
  mutate(
    # Generate predicted sine wave values
    eulachon_sine_pred = predict(sine_fit, newdata = .),
    # Enforce non-negativity constraint (SSB cannot be < 0)
    eulachon_sine_pred = pmax(0, eulachon_sine_pred)
  )

# 5. Plot Observed vs Extended Sine Wave
ggplot(lre_dat_eulachon_sine, aes(x = year)) +
  # Extended Sine Curve (Full Time Series)
  geom_line(aes(y = eulachon_sine_pred, color = "Extrapolated Sine Curve"), size = 1.2) +
  # Observed Data Points (2011 onward)
  geom_point(aes(y = eulachon_ssb_est, color = "Observed Data"), size = 3) +
  # Methodological transition line
  geom_vline(xintercept = 2011, linetype = "dashed", color = "firebrick", size = 0.8) +
  annotate("text", x = 2011.2, y = max(lre_dat_eulachon_sine$eulachon_sine_pred) * 0.9, 
           label = "2011 Survey Start", color = "firebrick", fontface = "bold") +
  scale_color_manual(values = c("Observed Data" = "black", "Extrapolated Sine Curve" = "darkgreen")) +
  labs(
    title = "Eulachon SSB Sine Wave Fit & Back-Extrapolation",
    subtitle = "Harmonic model fitted to post-2011 data and projected back to start of series",
    x = "Year",
    y = "Eulachon SSB Estimate",
    color = "Series"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

