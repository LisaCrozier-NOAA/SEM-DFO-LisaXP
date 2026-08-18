#Recalculate annual metrics
#Jake only used visual counts, which produced lots of missing weeks and years, and seemingly wildly inaccurate annual totals 
#(e.g., 2023, when there as only 1 visual survey w/ 14, but 6 uav that said 709). 
#I included all census types (airplane, uav, visual), so my annuals dont match his, but they seem much more realistic.
#I also broadened "spring" to 10-26 weeks for more samples


library(lubridate)
path<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices/LisaDataProcessScripts 2025/mammals/ODFW atlas count Columbia River v 20250218.xlsx"

odfw <- read_xlsx(path = path, sheet = 1)  %>%
    filter(spp == "ZC", 
     #    type == "VISUAL", 
         location == "COLUMBIA RIVER-EAST MOORING BASIN") %>%
    mutate(
      # 1. Convert integer YYYYMMDD (e.g. 19830531) to standard Date
      date_parsed = ymd(datemil),
      
      # 2. Extract calendar year and week number
      year = year(date_parsed),
      week = week(date_parsed)   # Use isoweek(date_parsed) if aligned to ISO weeks
    )




    
  csl_dat_year_jake <- odfw  %>% 
    filter(week %in% 14:25,
           type=="VISUAL") %>%
  # Group by year, and week
  group_by(year, type) %>% 
  # Calculate metrics per location-year-week
  summarize(
    csl_nonpup_total_emb = floor(mean(nonpup_total, na.rm = TRUE)),
    csl_nobs_emb = n(),  # Number of observations
    .groups = "drop"
  ) %>%
  select(-type)
  
  t(csl_dat_year_jake)
  t(csl_dat_year_jake)
  
  # year                 1998 2000 2001 2002 2003 2004 2005 2006 2007  2008  2009  2010  2011  2012  2013  2014  2015  2016  2017  2018  2021  2022  2023
  # csl_nonpup_total_emb   46   30   83  107   89   40   67   58   67    64    65   115   102    97   371   406   727   692   570   617   543   226    41
  # csl_nobs_emb            8   32   28   44   16   34   66   39   33    70    66    64    55    47    39    51    39    38    45    12     7     7     1
  
  

# # Reshape to wide format with separate columns for each census type
csl_dat <- odfw  %>% 
  filter(week %in% 10:26,year>=1998) %>%
  # Group by year, and week
  group_by(year, type) %>% 
  # Calculate metrics per location-year-week
  summarize(
    csl_nonpup_total_emb = floor(mean(nonpup_total, na.rm = TRUE)),
    csl_nobs_emb = n(),  # Number of observations
    .groups = "drop"
  )  %>%
  pivot_wider(
    id_cols = c(year),
    names_from = type,
    values_from = c(csl_nonpup_total_emb, csl_nobs_emb),
    names_sep = "_"
  )

t(csl_dat)

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

csl_year <- odfw  %>% 
  filter(week %in% 10:26,year>=1998) %>%
  # Group by year, and week
  group_by(year) %>% 
  # Calculate metrics per location-year-week
  summarize(
    csl_nonpup_total_emb = floor(mean(nonpup_total, na.rm = TRUE)),
    csl_nobs_emb = n(),  # Number of observations
    .groups = "drop"
  ) 

csl_week <- odfw  %>% 
  filter(week %in% 10:26,year>=1998) %>%
  # Group by year, and week
  group_by(year,week) %>% 
  # Calculate metrics per location-year-week
  summarize(
    csl_nonpup_total_emb = floor(mean(nonpup_total, na.rm = TRUE)),
    csl_nobs_emb = n(),  # Number of observations
    .groups = "drop"
  ) 


write.csv(
  csl_week,
  file.path(output_dir, "csl_week_observed.csv"),
  row.names = FALSE
)
write.csv(
  csl_year,
  file.path(output_dir, "csl_year_observed.csv"),
  row.names = FALSE
)

weekly_df <- left_join(weekly_df, csl_dat_weekly, by = c("year", "week"))
yearly_df <- left_join(yearly_df, csl_dat_yearly, by = c("year"))

  