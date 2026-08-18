#go back to full time series
output_dir <- "copilot/outputs_csl_cr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Load jake's data-----------------------------------------------------------------------------
jake_path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Jake Marshall/"
load(paste0(jake_path,"lre_dat_yearly.RData"),verbose=T)
head(lre_dat_yearly)
names(lre_dat_yearly)

head(lre_dat_eulachon_sine)
lre<-lre_dat_eulachon_sine

Eulachon_data_Gustavson_2010_status_review.xlsx

path<-"data_Lisa/Eulachon_data_Gustavson_2010_status_review.xlsx"

 eulachon_count<-read_xlsx(path=path,sheet=2,skip=1)  %>% 
   rename(eulachon_cr_pounds=`Total landings\r\n(pounds)`,
   eulachon_cr_nfish_10.8_per_pound=`Number of fish at\r\n10.8 per pound`,
    eulachon_cr_nfish_12.3_per_pound=`Number of fish at\r\n12.3 per pound`)
 head(eulachon_count)

 
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
 
 #Change Jake's processing---------
 jake_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Jake Marshall/"
 load(paste0(jake_path, "lre_dat_yearly.RData"), verbose = TRUE)

 odfw_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices/LisaDataProcessScripts 2025/mammals/ODFW atlas count Columbia River v 20250218.xlsx"

 odfw <- read_xlsx(path = odfw_path, sheet = 1) %>%
   filter(
     spp == "ZC",
     location == "COLUMBIA RIVER-EAST MOORING BASIN"
   ) %>%
   mutate(
     date_parsed = ymd(datemil),
     year        = year(date_parsed),
     week        = week(date_parsed)
   )

 # Calculate Lisa's Spring Annual CSL Mean (Weeks 10–26)
 csl_year <- odfw %>%
   filter(week %in% 10:26) %>%
   group_by(year) %>%
   summarize(
     csl_nonpup_total_emb_lisa = floor(mean(nonpup_total, na.rm = TRUE)),
     csl_nobs_emb_lisa         = n(),
     .groups                   = "drop"
   )

 # Calculate Lisa's Spring Weekly CSL Mean (Weeks 10–26)
 csl_week <- odfw %>%
 #  filter(week %in% 10:26) %>%
   group_by(year, week) %>%
   summarize(
     csl_nonpup_total_emb_lisa = floor(mean(nonpup_total, na.rm = TRUE)),
     csl_nobs_emb_lisa         = n(),
     .groups                   = "drop"
   )
 
 write.csv(
   csl_week,
   file.path(output_dir, "csl_week_all_census.csv"),
   row.names = FALSE
 )
 
 # 
 # # Merge with long-term yearly dataset
 # lre_dat_yearly_lisa <- left_join(lre_dat_yearly, csl_year, by = "year")
 # 
 # # -----------------------------------------------------------------------------
 # # Part 2: Extract & Format 50-Year Annual CSL Baseline (1976–2024)--------
 # # -----------------------------------------------------------------------------
 # csl_annual_baseline_50yr <- lre_dat_yearly_lisa %>%
 #   select(year, csl_nonpup_total_emb_lisa, eulachon_ssb_est) %>%
 #   filter(year >= 1976, year <= 2024) %>%
 #   rename(csl_annual_mean = csl_nonpup_total_emb_lisa) %>%
 #   complete(year = 1976:2024) %>%
 #   arrange(year) %>%
 #   mutate(
 #     csl_annual_mean  = zoo::na.approx(csl_annual_mean, na.rm = FALSE, rule = 2),
 #     log_annual_scale = log(pmax(1, csl_annual_mean))
 #   )
 # 
 # cat("--- Long-Term CSL Annual Baseline (1976–2024: Weeks 10–26 Spring Mean) ---\n")
 # print(csl_annual_baseline_50yr, n = 50)
 # 
 # write.csv(
 #   csl_annual_baseline_50yr,
 #   file.path(output_dir, "csl_annual_baseline_1976_2024.csv"),
 #   row.names = FALSE
 # )
 
 csl_annual_baseline_50yr <- read.csv(file.path(output_dir, "csl_annual_baseline_1976_2024.csv"),row.names = NULL)
 head(csl_annual_baseline_50yr)
 
 
 #Jake vs me---------
 
 # year                          1998 1999 2000 2001 2002 2003 2004 2005 2006  2007  2008  2009  2010  2011  2012  2013  2014  2015  2016  2017  2018  2020  2021  2022  2023  2024
 # csl_nonpup_total_emb_AIRPLANE    0   20   23    7   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
 # csl_nonpup_total_emb_VISUAL     59   NA   30   83  113  133   45   64   53    59    53    59   100    89    87   342   442   935   882   550   567    NA   413   226    41    NA
 # csl_nonpup_total_emb_UAV        NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA   576   480   396   709   307
 # csl_nobs_emb_AIRPLANE            1    3    1    2   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
 # csl_nobs_emb_VISUAL             14   NA   48   42   50   22   46   91   50    46    99    87    88    73    59    56    70    54    55    64    14    NA    12     7     1    NA
 # csl_nobs_emb_UAV                NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA     1     2     3     6     5
 
 t(csl_year)
 # year                 1998 1999 2000 2001 2002 2003 2004 2005 2006  2007  2008  2009  2010  2011  2012  2013  2014  2015  2016  2017  2018  2020  2021  2022  2023  2024
 # csl_nonpup_total_emb   55   20   30   79  113  133   45   64   53    59    53    59   100    89    87   342   442   935   882   550   567   576   423   277   614   307
 # csl_nobs_emb           15    3   49   44   50   22   46   91   50    46    99    87    88    73    59    56    70    54    55    64    14     1    14    10     7     5
 
 t(csl_dat_year_jake)
 
 # year                 1998 2000 2001 2002 2003 2004 2005 2006 2007  2008  2009  2010  2011  2012  2013  2014  2015  2016  2017  2018  2021  2022  2023
 # csl_nonpup_total_emb   46   30   83  107   89   40   67   58   67    64    65   115   102    97   371   406   727   692   570   617   543   226    41
 # csl_nobs_emb            8   32   28   44   16   34   66   39   33    70    66    64    55    47    39    51    39    38    45    12     7     7     1
 
 #intermediate products--------
 #note, eulachon_input was created in ssl_cr_predict_csl_phenology v2.r to fill in missing data:
 # UNGROUPED INTERPOLATION: Fills Eulachon continuously across newly created 2020 rows
 #eulachon_input       = zoo::na.approx(eulachon_ssb_4week_est, na.rm = FALSE, rule = 2),
 week_final_scaled<-read.csv("csl_week_filled_in_2011_2024.csv",row.names = NULL)

  #RESULTS:
 week_final_scaled<-read.csv("csl_week_filled_in_2011_2024.csv",row.names = NULL)
 #and remake plot:
 csl_weekly_fill_2011_2024_plot <-
   ggplot(week_final_scaled , aes(x = week)) +
   geom_ribbon(aes(ymin = csl_lwr_95, ymax = csl_upr_95), fill = "steelblue", alpha = 0.25) +
   geom_line(aes(y = csl_gam_pred, color = "Reconstructed Curve"), linewidth = 1) +
   geom_point(aes(y = csl_nonpup_total_emb, color = "Observed Points"), size = 2.5, na.rm = TRUE) +
   facet_wrap(~ year, ncol = 4, scales = "free_y") +
   scale_color_manual(values = c("Observed Points" = "black", "Reconstructed Curve" = "steelblue")) +
   labs(
     title = "Reconstructed Sea Lion Phenology Across COVID Gap (2018–2021)",
     subtitle = "2020 explicitly expanded via complete() and predicted from GAM",
     x = "Week Number (Apr–Jun)",
     y = "CSL Non-Pup Count",
     color = "Legend"
   ) +
   theme_minimal(base_size = 13) +
   theme(legend.position = "top")
 print(csl_weekly_fill_2011_2024_plot)
 
 