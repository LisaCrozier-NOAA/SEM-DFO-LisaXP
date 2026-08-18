########################################################
########################################################
## Processing Data from 2011-2024
## 1. Adult passage at Bonneville
##      - Shad
##      - Lamprey - taken out -- makes up > 1% of biomass of only 8% of weeks across time period and is nearly always under 5%
##      - Adult Steelhead
##      - Adult Chinook
## 2. Juvenile passage at Bonneville
##      - Juvenile Steelhead
##      - Juvenile Yearling Chinook
##      - Juvenile Subyearling Chinook 
## 3. Eulachon
##    a. Eulachon SSB Estimates from WDFW - weekly mean
## 4. Harbor Seals
##    a. counts at Desdemona Sands and South Jetty - weekly mean
##    b. abundance-based diet rank
##    c. biomass-based diet rank
## 5. CSL count at EMB - weekly mean 
## 6. SSL count at SJ and EMB - weekly mean
## 7. CSL and SSL at Bonn


## Graveyard -- old data not used in final analysis
##    Eulachon harvest
##    CPS presence offshore of river mouth - weekly mean
##    Hake Harvest - Annual index
##    River Covariates


########################################################
########################################################

setwd('C:/Users/jakem/OneDrive/Desktop/salmon_R/')

# make time series grid from 1976 through 2024 to fit all variables into
weekly_df <- tibble(
  date = seq(as.Date("1976-01-01"), as.Date("2024-12-17"), by = "weeks")) %>%
  mutate(
    year = year(date),
    month = month(date),
    doy = yday(date),
    week = ceiling(doy/7) 
  ) %>% 
  select(year, week, month, date)

# make time series grid from 1976 through 2024 to fit all variables into
yearly_df <- tibble(
  year = 1976:2024
) 


#########################################################################################
#########################################################################################
## 1. Shad Passage at Bonneville
#########################################################################################
#########################################################################################

prey_passage <- read_excel("./data/fish/Adult_passage_1976_2024.xlsx")

# convert time and select useful variables
prey_passage_dat <- prey_passage %>%
  mutate(date = as.Date(paste0(year, "-", substr(`mm-dd`, 6, 10)), format="%Y-%m-%d"),
         year = year(date),
         doy = yday(date),
         week = ceiling(doy/7), 
         # deal with "NA" counts without throwing warning
         count = as.numeric(ifelse(count == "NA", NA, count))) %>% 
  select(date, year, week, parameter, count)
  

# Calculate species-specific weekly passage totals
prey_passage_weekly <- prey_passage_dat %>% 
  group_by(year, week) %>% 
  summarize(
    # Adult_chin_bonn_pass = sum(if_else(parameter == "Chin", 
    #                                    ifelse(is.na(count), 0, count), 
    #                                    0), na.rm = TRUE), # assume 0 passage if NA
    # Adult_stlhd_bonn_pass = sum(if_else(parameter == "Stlhd", 
    #                                     ifelse(is.na(count), 0, count), 
    #                                     0), na.rm = TRUE),
    Lmpry_bonn_pass = sum(if_else(parameter == "Lmpry",
                                  ifelse(is.na(count), 0, count),
                                  0), na.rm = TRUE),
    Shad_bonn_pass = sum(if_else(parameter == "Shad", 
                                 ifelse(is.na(count), 0, count), 
                                 0), na.rm = TRUE),
    .groups = "drop"
  )
  
# Calculate species-specific yearly passage totals
prey_passage_yearly <- prey_passage_weekly %>% 
  group_by(year) %>% 
  summarize(
    Shad_bonn_pass = sum(Shad_bonn_pass),
    Lmpry_bonn_pass = sum(Lmpry_bonn_pass),
    .groups = "drop"
  )


weekly_df <- left_join(weekly_df, prey_passage_weekly, by = c("year", "week"))
yearly_df <- left_join(yearly_df, prey_passage_yearly, by = c("year"))



#########################################################################################
## Adult Chinook and Steelhead passage at Bonneville------
#########################################################################################

onc_passage <- read_csv("./data/fish/Adult_oncorhynchus_bonn_1976_2024.csv")

# Pink and Chum runs are < 0.01% of all fish
# Jack Coho is 0.7% of all fish
# Coho is 7% of all fish
# Jack Chinook is 7% of all fish
# Sockeye is 15%
# Steelhead is 23.5% of all fish
# Chinook is 46.4% of all fish

# clean and make date variables
onc_passage_dat <- onc_passage %>%
  
  ## Remove Feb 29 from non-leap years
  mutate(
    date_string = paste(`mm-dd`, year),
    is_invalid_feb29 = grepl("Feb-29", date_string) & 
      !leap_year(as.numeric(str_extract(date_string, "\\d{4}$")))
  ) %>%
  # Remove those dates and change Feb-29 to 29-Feb for correct parsing of leap year days
  filter(!is_invalid_feb29) %>%
  mutate(date_string = ifelse(`mm-dd` == "Feb-29", paste("29-Feb", year), date_string)) %>% 
  # Parse the remaining (valid) dates
  mutate(
    date = dmy(date_string),
    year = year(date),
    month = month(date),
    doy = yday(date),
    week = ceiling(doy/7),
    
    # develop season according to CRITFIC definitions -- link in time with length measurements 
    season = case_when(
      week < week_of_june_1  ~ "spring",
      between(week, week_of_june_1, week_of_july_31) ~ "summer",
      week > week_of_july_31  ~ "fall"
      )
  ) %>% 
  select(date, year, week, doy, season, parameter, value) %>% 
  # select only chinook, steelhead, sockeye, coho
  filter(parameter %in% c("Chin", "JChin", "Stlhd"))
  

# Group together species passage totals weekly
onc_passage_weekly <- onc_passage_dat %>% 
  group_by(year, week) %>% 
  summarize(
    season = first(season),
    
    # All salmon
    oncorhynchus_bonn_pass = sum(value, na.rm = T),
    
    # Aggregate all chinook
    Adult_chin_bonn_pass = sum(if_else(parameter %in% c("Chin", "JChin"),
                                  ifelse(is.na(value), 0, value), 
                                  0), na.rm = TRUE), 
    
    # Achin is non-jack adults
    Achin_bonn_pass = sum(if_else(parameter %in% c("Chin"), ## checked -- correctly adds jack and adult chinook
                                       ifelse(is.na(value), 0, value), 
                                       0), na.rm = TRUE), # assume 0 passage if NA
    Jchin_bonn_pass = sum(if_else(parameter %in% c("JChin"),
                                      ifelse(is.na(value), 0, value),
                                      0), na.rm = TRUE),
    
    # Aggregate chinook grouped by season of run -- Totals and just Jacks
    Spring_Achin_bonn_pass = sum(if_else((parameter %in% c("Chin") & season == "spring"), 
                                       ifelse(is.na(value), 0, value), 
                                       0), na.rm = TRUE), 
    Spring_Jchin_bonn_pass = sum(if_else((parameter %in% c("JChin") & season == "spring"), 
                                        ifelse(is.na(value), 0, value), 
                                        0), na.rm = TRUE), 
    Summer_Achin_bonn_pass = sum(if_else((parameter %in% c("Chin") & season == "summer"), 
                                        ifelse(is.na(value), 0, value), 
                                        0), na.rm = TRUE),
    Summer_Jchin_bonn_pass = sum(if_else((parameter %in% c("JChin") & season == "summer"), 
                                        ifelse(is.na(value), 0, value), 
                                        0), na.rm = TRUE), 
    Fall_Achin_bonn_pass = sum(if_else((parameter %in% c("Chin") & season == "fall"), 
                                        ifelse(is.na(value), 0, value), 
                                        0), na.rm = TRUE), 
    Fall_Jchin_bonn_pass = sum(if_else((parameter %in% c("JChin") & season == "fall"), 
                                      ifelse(is.na(value), 0, value), 
                                      0), na.rm = TRUE),
    # checked that sum of seasonal chinook variables is equal to total Adult_chin_bonn_pass (jack and adult)
    
    # Calculate Steelhead passage
    Adult_stlhd_bonn_pass = sum(if_else(parameter == "Stlhd", 
                                        ifelse(is.na(value), 0, value), 
                                        0), na.rm = TRUE),
    
    .groups = "drop"
  ) %>% 
  select(year, week, season, oncorhynchus_bonn_pass, 
         Adult_chin_bonn_pass,
         Achin_bonn_pass, Jchin_bonn_pass, 
         Spring_Achin_bonn_pass, Spring_Jchin_bonn_pass, 
         Summer_Achin_bonn_pass, Summer_Jchin_bonn_pass,
         Fall_Achin_bonn_pass, Fall_Jchin_bonn_pass,
         Adult_stlhd_bonn_pass)


# Calculate yearly passage totals of oncorhynchus
onc_passage_yearly <- onc_passage_dat %>% 
  group_by(year) %>% 
  summarize(
    
    season = first(season),
    
    # All salmon
    oncorhynchus_bonn_pass = sum(value, na.rm = T),
    
    # Aggregate all chinook
    Adult_chin_bonn_pass = sum(if_else(parameter %in% c("Chin", "JChin"),
                                       ifelse(is.na(value), 0, value), 
                                       0), na.rm = TRUE), 
    # Achin is non-jack adults
    Achin_bonn_pass = sum(if_else(parameter %in% c("Chin"), 
                                       ifelse(is.na(value), 0, value), 
                                       0), na.rm = TRUE),
    Jchin_bonn_pass = sum(if_else(parameter %in% c("JChin"),
                                      ifelse(is.na(value), 0, value),
                                      0), na.rm = TRUE),
    
    # Aggregate chinook grouped by season of run -- Totals and just Jacks
    Spring_Achin_bonn_pass = sum(if_else((parameter %in% c("Chin") & season == "spring"), 
                                        ifelse(is.na(value), 0, value), 
                                        0), na.rm = TRUE), 
    Spring_Jchin_bonn_pass = sum(if_else((parameter %in% c("JChin") & season == "spring"), 
                                            ifelse(is.na(value), 0, value), 
                                            0), na.rm = TRUE), 
    Summer_Achin_bonn_pass = sum(if_else((parameter %in% c("Chin") & season == "summer"), 
                                        ifelse(is.na(value), 0, value), 
                                        0), na.rm = TRUE),
    Summer_Jchin_bonn_pass = sum(if_else((parameter %in% c("JChin") & season == "summer"), 
                                            ifelse(is.na(value), 0, value), 
                                            0), na.rm = TRUE), 
    Fall_Achin_bonn_pass = sum(if_else((parameter %in% c("Chin") & season == "fall"), 
                                      ifelse(is.na(value), 0, value), 
                                      0), na.rm = TRUE), 
    Fall_Jchin_bonn_pass = sum(if_else((parameter %in% c("JChin") & season == "fall"), 
                                          ifelse(is.na(value), 0, value), 
                                          0), na.rm = TRUE),
    
    # # calculate jack proportion in each season
    # Jack_chin_proportion = Jack_chin_bonn_pass / Adult_chin_bonn_pass,
    # Spring_JChin_proportion = JackSpring_chin_bonn_pass / Spring_Achin_bonn_pass,
    # Summer_Jhin_proportion = JackSummer_chin_bonn_pass / Summer_Achin_bonn_pass,
    # Fall_Jhin_proportion = JackFall_chin_bonn_pass / Fall_Achin_bonn_pass,
    
    # Calculate Steelhead passage
    Adult_stlhd_bonn_pass = sum(if_else(parameter == "Stlhd", 
                                        ifelse(is.na(value), 0, as.numeric(value)), 
                                        0), na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  select(year, oncorhynchus_bonn_pass, 
         Adult_chin_bonn_pass,
         Achin_bonn_pass, Jchin_bonn_pass, 
         Spring_Achin_bonn_pass, Spring_Jchin_bonn_pass, 
         Summer_Achin_bonn_pass, Summer_Jchin_bonn_pass,
         Fall_Achin_bonn_pass, Fall_Jchin_bonn_pass,
         Adult_stlhd_bonn_pass)

weekly_df <- left_join(weekly_df, onc_passage_weekly, by = c("year", "week"))
yearly_df <- left_join(yearly_df, onc_passage_yearly, by = c("year"))




#########################################################################################
#########################################################################################
## 2. Juvenile passage at Bonneville-----
#########################################################################################
#########################################################################################

# add data together into one data frame
juve_passage_dat1 <- read.csv("./data/fish/Smolt_index_onerow_2010_2024.csv")
juve_passage_dat2 <- read.csv("./data/fish/Smolt_index_onerow_2004_2009.csv")
juve_passage_dat3 <- read.csv("./data/fish/Smolt_index_onerow_1985_2003.csv")
juve_passage <- rbind(juve_passage_dat1, juve_passage_dat2, juve_passage_dat3)

# clean up data
juve_passage_dat <- juve_passage %>%
  
  ## Remove Feb 29 from non-leap years
  mutate(
    date_string = paste(mm.dd, year),
    is_invalid_feb29 = grepl("Feb-29", date_string) & 
      !leap_year(as.numeric(str_extract(date_string, "\\d{4}$")))
  ) %>%
  # Remove those dates and change Feb-29 to 29-Feb for correct parsing of leap year days
  filter(!is_invalid_feb29) %>%
  mutate(date_string = ifelse(mm.dd == "Feb-29", paste("29-Feb", year), date_string)) %>% 
  # Parse the remaining (valid) dates
  mutate(
    date = dmy(date_string),
    year = year(date),
    month = month(date),
    doy = yday(date),
    week = ceiling(doy/7)
  ) %>% 
  select(year, week, parameter, value) 
  
  
# Calculate species-specific weekly juvenile passage
juve_passage_weekly <- juve_passage_dat %>% 
  group_by(year, week) %>% 
  summarize(
    Ch1_bonn_pass = sum(if_else(parameter == "Ch1", 
                                 ifelse(is.na(value), 0, value), # turn NAs into 0s
                                 0), na.rm = TRUE),
    Ch0_bonn_pass = sum(if_else(parameter == "Ch0",  
                                 ifelse(is.na(value), 0, value), 
                                 0), na.rm = TRUE),
    Juve_stlhd_bonn_pass = sum(if_else(parameter == "Sthd",  
                                        ifelse(is.na(value), 0, value), 
                                        0), na.rm = TRUE),
    .groups = "drop"
  ) 

# Calculate species-specific yearly juvenile passage
juve_passage_yearly <- juve_passage_weekly %>% 
  group_by(year) %>% 
  summarize(
    Ch1_bonn_pass = sum(Ch1_bonn_pass, na.rm = TRUE),
    Ch0_bonn_pass = sum(Ch0_bonn_pass, na.rm = TRUE),
    Juve_stlhd_bonn_pass = sum(Juve_stlhd_bonn_pass, na.rm = TRUE),
    .groups = "drop"
  ) 


weekly_df <- left_join(weekly_df, juve_passage_weekly, by = c("year", "week"))
yearly_df <- left_join(yearly_df, juve_passage_yearly, by = c("year"))



#########################################################################################
#########################################################################################
## 3. Eulachon SSB estimates---------
#########################################################################################
########################################################################################


# new SSB estimates from WDFW
eul_dat <- read_excel('./data/fish/CR_Weekly_Eulachon_SSB_Estimates_2011_2024.xlsx') %>%
  mutate(
    date = `Date -4 weeks`, # use 4-week lag as date for eulachon ssb in estuary
    year = year(date),
    doy = yday(date),
    week = ceiling(doy/7)
  )

weekly_eul_dat <- eul_dat %>%
  group_by(year, week) %>%
  summarize(
    eulachon_ssb_4week_est = mean(`Weekly SSB Estimate`, na.rm = T),
    year = first(year),
    week = first(week),
    .groups = "drop"
  )
# rename(eulachon_ssb_4week_est = `Weekly SSB Estimate`)

yearly_eul_dat <- eul_dat %>%
  group_by(year) %>%
  summarize(
    eulachon_ssb_est = sum(`Weekly SSB Estimate`),
    year = first(year)
  )


weekly_df <- left_join(weekly_df,
                       weekly_eul_dat,
                       by = c("year", "week")) %>% 
  mutate(eulachon_ssb_4week_est = replace_na(eulachon_ssb_4week_est, 0))

yearly_df <- left_join(yearly_df,
                       yearly_eul_dat,
                       by = c("year"))


#########################################################################################
#########################################################################################
# CPS Stock Biomass Estimate
#########################################################################################
#########################################################################################

# Estimated Stock Biomass
load("./data/fish/biomass_timeseries_final.Rdata")
spp_northern_names <- c("Engraulis mordax", "Sardinops sagax")
spp_all_names <- c("Clupea pallasii", "Trachurus symmetricus")

# Configure data for merging with annual time series
cps_bmass <- biomass.ts %>%
  ungroup() %>%
  # Select for stocks we're interested in
  filter((species %in% spp_northern_names & stock == "Northern") | (species %in% spp_all_names & stock == "All")) %>%
  select(year, species, stock, biomass) %>%

  # make variables for each species
  pivot_wider(names_from = c(species, stock),
            values_from = biomass,
            names_glue = "{species}_{stock}_biomass") %>%

  # Rename variables to simpler names
  rename_with(~ case_when(
    str_detect(.x, "Sardinops") ~ "sardine_northern_bmass",
    str_detect(.x, "Engraulis") ~ "anchovy_northern_bmass",
    str_detect(.x, "Clupea") ~ "herring_bmass",
    str_detect(.x, "Trachurus") ~ "jack_mackerel_bmass",
    TRUE ~ .x
  ))


# northern anchovy, herring variables have nothing before 2015, so interpolation doesn't work

# Merge with dfs
yearly_df <- left_join(yearly_df, cps_bmass, by = c("year"))




#########################################################################################
#########################################################################################
## Hake Harvest
#########################################################################################
#########################################################################################

# spatial coverage goes down to a little south of newport
hake <- read.csv('./data/fish/hake_harvest.csv')

lat.max <- 47.3
lat.min <- 45.3


# take all hake harvest data in the latitude range
hake_dat <- hake %>%
  mutate(latitude = mean(c(lat, lat2)),
         longitude = mean(c(lon, lon2)),
         date = as.Date(substr(DEPLOYMENT_DATE, 1, 10), format="%Y-%m-%d"),
         week = week(date)) %>%
  filter(between(latitude, lat.min, lat.max)) %>%

  # Calculate yearly means of aggregated CPUE across range
  group_by(year) %>%
  summarize(
    hake_cpue = 1000 * mean(cpue, na.rm = TRUE),
    # CPUE is in units of metric tons/minute
    # above line transforms CPUE to kg/minute
    .groups = "drop"
  )

yearly_df <- left_join(yearly_df, hake_dat, by = c("year"))


#########################################################################################
#########################################################################################
## Sturgeon Abundance below Bonneville Dam
#########################################################################################
#########################################################################################

#Estimated and projected abundance of 42–60-inch total length (38–54-inch fork length) white sturgeon in the lower Columbia River, 1987–2024.
sturg <- read.csv('./data/fish/Sturgeon_abund.csv') 

yearly_df <- left_join(yearly_df, sturg, by = c("year"))



#########################################################################################
#########################################################################################
## 4. HS counts
#########################################################################################
#########################################################################################

# ODFW data
setwd('C:/Users/jakem/OneDrive/Desktop/salmon_R/')
odfw_dat <- read.csv('./data/env/ODFW_Pinn_Counts.csv')

## Configure date into year, month, week variables, and select for Desdemona Sands only
odfw_dat <- odfw_dat %>% 
  mutate(
    dates = as.Date(as.character(datemil), format = "%Y%m%d"),
    year = year(dates),
    month = month(dates),
    doy = yday(dates),
    week = ceiling(doy/7),
    dates = as.Date(as.character(datemil), format = "%Y%m%d")
  )


# filter for HS counts at Desdemona Sands taken with airplane method (to match WDFW's methods)
odfw_combine <- odfw_dat %>% 
  mutate(location = str_sub(location, 16, -1), # rename locations for ease
         source = 1 # ODFW source indicator
  ) %>% 
  filter(location %in% c("DESDEMONA SANDS"), spp == "PV", type == "AIRPLANE") %>% 
  select(year, month, doy, week, spp, nonpup_total, source)



# # look at observation types
table(odfw_dat$type[which(odfw_dat$spp == "PV")], odfw_dat$location[which(odfw_dat$spp == "PV")])
# # most visual observations were taken at EMB
# # all SJ and DS observations before 2023 were taken using aerial methods
# # nearly all EMB observations were taken using visual methods



# WDFW data
wdfw_hs <- read.csv('./data/env/WDFW_Pinn_Counts.csv')
site_codes <- read_xlsx('./data/env/WDFW_site_code_info.xlsx')

# add site names to observations and change variable names
wdfw_combine <- left_join(wdfw_hs %>% select(Sitecode, Julian, Year, Day, Species, 
                                        Count.nonpup, Count.total, Count.pups), 
                     site_codes %>% select(`Site code`, Area) %>% 
                       rename(Sitecode = `Site code`) %>% 
                       distinct(Sitecode, .keep_all = TRUE),  # Keep only unique site code
                     by = c("Sitecode")) %>% 
  
  mutate(month = month(as.Date(Day, format = "%m/%d/%Y")),
         week = ceiling(Julian/7)) %>% 
  # select and rename key variables to fit ODFW dataset
  select(Area, Year, month, Julian, week, Species, Count.nonpup, Count.total, Count.pups) %>% 
  rename(location = Area, year = Year, doy = Julian, spp = Species, nonpup_total =  Count.nonpup, 
         total = Count.total, pup_total = Count.pups) %>% 
  # change location to uppercase to match odfw
  mutate(location = factor(toupper(location), levels = c("DESDEMONA SANDS",            
                                                         "TAYLOR SANDS",               
                                                         "SOUTH OF MILLER SANDS",   
                                                         "CATHLAMET BAY",              
                                                         "E PUGET ISLAND",             
                                                         "GRAYS BAY" ,                 
                                                         "NE OF WELCH ISLAND",         
                                                         "TONGUE PT. SANDS",           
                                                         "LOWER WOODY/SNAG IS.",       
                                                         "THREE TREE POINT AREA",      
                                                         "CHINOOK/BAKER BAY",          
                                                         "TIP OF SOUTH JETTY",         
                                                         "WALLACE IS./EUREKA BAR",     
                                                         "LARGE NAVIGATION BOUY",      
                                                         "COWLITZ R./CARROLL SLOUGH",  
                                                         "HAULED OUT UNKNOWN",         
                                                         "ASTORIA/E MOORING BASIN",    
                                                         "PHOCA ROCK (COLUMBIA RIVER)")), 
         source = 0
  ) %>% 
  # filter for HS counts at Desdemona Sands taken with airplane method (to match WDFW's methods)
  filter(location %in% c("DESDEMONA SANDS"), spp == "PV") %>% 
  select(year, month, doy, week, spp, nonpup_total, source)


# Combine data from ODFW and WDFW
hs_dat_combined <- rbind(odfw_combine, wdfw_combine) 

# Clean 
# -- take out very low counts during weeks with high counts
hs_dat_cleaned <- hs_dat_combined %>% 
  filter(week %in% 14:25) %>%
  group_by(year, week) %>%
  mutate(
    week_max = max(nonpup_total, na.rm = TRUE),
    # Remove near-zero values if there are substantial values that week
    keep_row = case_when(
      # Remove points if <30% of week max
      week_max > 100 & nonpup_total < (week_max * 0.3) ~ FALSE, 
      TRUE ~ TRUE
    )
  ) %>%
  filter(keep_row) %>%
  select(-keep_row) %>%
  ungroup()

# # # uncomment filter and select above to plot removed points and percentage removed
# ggplot(hs_dat_cleaned)+
#   geom_point(aes(x = week, y = nonpup_total, color = keep_row))+
#   facet_wrap(~year, nrow = 3)
# 
# 1 - (sum(hs_dat_cleaned$keep_row) / nrow(hs_dat_cleaned)) # 11% of points removed


# Form yearly data frame for mean weekly seal counts and survey source ratio
hs_dat_yearly <- hs_dat_cleaned %>% 
  # average counts over each year for each location
  group_by(year) %>% 
  summarize(
    hs_nonpup_total_ds = floor(mean(nonpup_total, na.rm = TRUE)),
    hs_odfw_ratio_ds = sum(source)/n(),  # Proportion of ODFW observations
    hs_nobs_ds = n(), 
    .groups = "drop"
  )

yearly_df <- left_join(yearly_df, hs_dat_yearly, by = c("year"))


# # Form weekly data frame for mean weekly seal counts and survey source ratio
# hs_dat_weekly <- hs_dat_cleaned %>% 
#   # average counts over each year for each location
#   group_by(year, week) %>% 
#   summarize(
#     hs_nonpup_total_ds = floor(mean(nonpup_total, na.rm = TRUE)),
#     hs_odfw_ratio_ds = sum(source)/n(),  # Proportion of ODFW observations
#     hs_nobs_ds = n(),
#     .groups = "drop"
#   ) 

# add to df
# weekly_df <- left_join(weekly_df, hs_dat_weekly, by = c("year", "week"))


## issue: don't have times for ODFW data, but do have time and tide for WDFW (none of which are higher than 2.61 ft? m?)
## -- RESOLVED: Bryan Wright said that all spring can be assumed to be +- 2 hours of low tide at Desdemona Sands
##              - see email from him




#########################################################################################
#########################################################################################
## 5. CSL at EMB
#########################################################################################
#########################################################################################

###   all visual counts of CSL at East Mooring Basin 


# filter ODFW pinniped data for CSL at EMB
csl_dat_weekly <- odfw_dat  %>% 
  filter(spp == "ZC", 
         type == "VISUAL", 
         location == "COLUMBIA RIVER-EAST MOORING BASIN") %>% 
  select(-location) %>% 
  # Group by year, and week
  group_by(year, week) %>% 
  # Calculate metrics per location-year-week
  summarize(
    csl_nonpup_total_emb = floor(mean(nonpup_total, na.rm = TRUE)),
    csl_nobs_emb = n(),  # Number of observations
    .groups = "drop"
  ) 

# quantify yearly as spring for now
csl_dat_yearly <- odfw_dat %>% 
  filter(spp == "ZC", 
         type == "VISUAL", 
         location == "COLUMBIA RIVER-EAST MOORING BASIN",
         week  %in% 14:25) %>% 
  select(-location) %>% 
  # Group by year
  group_by(year) %>% 
  # Calculate metrics per location-year-week
  summarize(
    csl_nonpup_total_emb = floor(mean(nonpup_total, na.rm = TRUE)),
    csl_nobs_emb = n(),  # Number of observations
    .groups = "drop"
  )

# # Reshape to wide format with separate columns for each location
# csl_dat <- csl_dat %>%
#   pivot_wider(
#     id_cols = c(year, week),
#     names_from = location,
#     values_from = c(csl_nonpup_total, csl_nobs),
#     names_sep = "_"
#   )

weekly_df <- left_join(weekly_df, csl_dat_weekly, by = c("year", "week"))
yearly_df <- left_join(yearly_df, csl_dat_yearly, by = c("year"))

#########################################################################################
#########################################################################################
## 6. SSL at SJ
#########################################################################################
#########################################################################################

###   all aerial counts of SSL at South Jetty ??? between 15th and 25th week of year
# average on yearly basis
ssl_dat_yearly <- odfw_dat %>%
  filter(spp == "EJ", 
         type == "AIRPLANE",
         location == "COLUMBIA RIVER-SOUTH JETTY") %>% 

  # Group by location, year, and week
  group_by(year) %>% 
  
  # Calculate metrics
  summarize(
    ssl_nonpup_total_sj = floor(mean(nonpup_total, na.rm = TRUE)),
    ssl_nobs_sj = n(),  # Number of observations
    .groups = "drop"
  ) 

# average on weekly basis
ssl_dat_weekly <- odfw_dat %>%
  filter(spp == "EJ", 
         type == "AIRPLANE",
         location == "COLUMBIA RIVER-SOUTH JETTY") %>% 
  
  # Group by location, year, and week
  group_by(year, week) %>% 
  
  # Calculate metrics
  summarize(
    ssl_nonpup_total_sj = floor(mean(nonpup_total, na.rm = TRUE)),
    ssl_nobs_sj = n(),  # Number of observations
    .groups = "drop"
  ) 

yearly_df <- left_join(yearly_df, ssl_dat_yearly, by = c("year"))
weekly_df <- left_join(weekly_df, ssl_dat_weekly, by = c("year", "week"))



# #########################################################################################
# #########################################################################################
# # Create data.frame with raw pinniped counts
# #########################################################################################
# #########################################################################################
# 
# raw_pinn_dat <- hs_dat_combined %>% mutate(location = "Desdemona Sands") %>% select(-source)# hs_dat_combined is all hs counts at Desdemona Sands from ODFW and WDFW
# 
# ssl_sj_dat <- odfw_dat %>%
#   filter(spp == "EJ", 
#          type == "AIRPLANE",
#          location == "COLUMBIA RIVER-SOUTH JETTY") %>% 
#   select(year, month, doy, week, spp, nonpup_total, location)
#   
# csl_emb_dat <- odfw_dat  %>% 
#   filter(spp == "ZC", 
#          type == "VISUAL", 
#          location == "COLUMBIA RIVER-EAST MOORING BASIN") %>% 
#   select(year, month, doy, week, spp, nonpup_total, location)
#   
# names(raw_pinn_dat)
# 
# 
# raw_pinn_dat <- rbind(raw_pinn_dat, ssl_sj_dat, csl_emb_dat)
# 
# save(raw_pinn_dat, file = "./data/env/raw_pinn_dat.RData")

#########################################################################################
#########################################################################################
## 7. CSL, SSL at Bonn
#########################################################################################
#########################################################################################

# From USACE
bonn_dat <- read_xlsx('./data/BON_pinniped_abundance_2002-2025.xlsx')

# take out interpolated data, all have yes or interpolated, except for special cases below
pattern <- paste(c("yes", "interpolated"), collapse = "|")

bonn_dat <- bonn_dat %>% 
  filter(!grepl(pattern, Interpolated, ignore.case = TRUE) & !Interpolated %in% c("no count performed","Interploated")) %>%
  mutate(date = substr(DATE, 1, 10),
         year = YEAR, 
         month = MONTH, 
         doy = yday(date),
         week = ceiling(doy/7)) 

bonn_dat_weekly <- bonn_dat %>% 
  group_by(year, week) %>% 
  summarize(
    csl_bonn = replace_na(floor(mean(`ZCA Abundance`, na.rm = T)), 0), # assume 0 count when survey isn't running
    ssl_bonn = replace_na(floor(mean(`EJU Abundance`, na.rm = T)), 0),
    total_bonn = replace_na(floor(mean(`Total Pinniped Abundance`, na.rm = T)), 0), 
    .groups = "drop"
  ) 
  
bonn_dat_yearly <- bonn_dat_weekly %>% 
  filter(week %in% 14:25) %>% 
  group_by(year) %>% 
    summarize(
      csl_bonn = floor(mean(csl_bonn)), 
      ssl_bonn = floor(mean(ssl_bonn)),
      total_bonn = floor(mean(total_bonn)), 
      .groups = "drop"
    ) 
  

weekly_df <- left_join(weekly_df, bonn_dat_weekly, by = c("year", "week"))
yearly_df <- left_join(yearly_df, bonn_dat_yearly, by = c("year"))
  

#
lethal_removals <- data.frame(
  year = 2008:2024,
  csl_removed = c(7,15,14,1,13,4,15,32,59,24,28,19,0,21,14,22,27),
  ssl_removed = c(0,0,0,0,0,0,0,0,0,0,0,0,6,37,9,25,17)
)

weekly_df <- left_join(weekly_df, lethal_removals, by = c("year"))
yearly_df <- left_join(yearly_df, lethal_removals, by = c("year"))


#########################################################################################
## save
#########################################################################################
lre_dat_weekly <- weekly_df
lre_dat_yearly <- yearly_df


save(lre_dat_weekly, file = './data/lre_dat_weekly.RData')
save(lre_dat_yearly, file = './data/lre_dat_yearly.RData')



# Couple of analyses

# take note of temporal coverage of each variable
first_nonzero_years <- lre_dat_yearly %>% 
  select(year, all_of(names(lre_dat_yearly)[c(2,3,5,14:17, 21:26)])) %>%
  pivot_longer(-year, names_to = "variable", values_to = "value") %>%
  filter(value > 0 & !is.na(value)) %>%
  group_by(variable) %>%
  summarise(first_nonzero_year = min(year), .groups = "drop")






# Migration window analysis



# CSL
# spring 2011 and 2012
lre_dat_weekly %>% filter(year %in% 2011:2012 & week  %in% 10:21) %>% pull(csl_nonpup_total_emb) %>% mean(.,na.rm = T)
# spring 2014 and 2015
lre_dat_weekly %>% filter(year %in% 2014:2015 & week  %in% 10:21) %>% pull(csl_nonpup_total_emb) %>% mean(.,na.rm = T)


### Adult Steelhead
weeks_migrating_sthd <- lre_dat_weekly %>%
  group_by(year) %>%
  filter(Adult_stlhd_bonn_pass > 50) %>%
  summarize(min_week = min(week),
            max_week = max(week))
range_sthd <- c(weeks_migrating_sthd %>% select(min_week) %>% min(), 
                weeks_migrating_sthd %>% select(max_week) %>% max())

# Winter Steelhead -- how many pass weekly before week 12
winter_sthd <- lre_dat_weekly %>% 
  filter(week<=10) %>% 
  summarize(
    mean_passage = mean(Adult_stlhd_bonn_pass), 
  )
lre_dat_weekly %>% select(year, week, Adult_stlhd_bonn_pass) %>% pivot_wider(names_from = year, values_from = Adult_stlhd_bonn_pass) %>% filter(week < 10)

summer_sthd <- lre_dat_weekly %>% 
  filter(week %in% 25:35) %>% 
  summarize(
    mean_passage = mean(Adult_stlhd_bonn_pass), 
  )

winter_chin <- lre_dat_weekly %>% 
  filter(week<=10) %>% 
  group_by(year) %>% 
  summarize(
    mean_passage = mean(Adult_chin_bonn_pass), 
  )
lre_dat_weekly %>% select(year, week, Adult_chin_bonn_pass) %>% pivot_wider(names_from = year, values_from = Adult_chin_bonn_pass) %>% filter(week > 40)


## Juveniles
# what months are juveniles usually migrating
weeks_migrating_ch1 <- lre_dat_weekly %>%
  group_by(year) %>%
  filter(Ch1_bonn_pass > 1000) %>%
  summarize(min_week = min(week),
            max_week = max(week))
range_ch1 <- c(weeks_migrating_ch1 %>% select(min_week) %>% min(), 
               weeks_migrating_ch1 %>% select(max_week) %>% max())

weeks_migrating_ch0 <- lre_dat_weekly %>%
  group_by(year) %>%
  filter(Ch0_bonn_pass > 1000) %>%
  summarize(min_week = min(week),
            max_week = max(week))
range_ch0 <- c(weeks_migrating_ch0 %>% select(min_week) %>% min(), 
               weeks_migrating_ch0 %>% select(max_week) %>% max())

weeks_migrating_sthd <- lre_dat_weekly %>%
  group_by(year) %>%
  filter(Juve_stlhd_bonn_pass > 1000) %>%
  summarize(min_week = min(week),
            max_week = max(week))
range_juve_sthd <- c(weeks_migrating_sthd %>% select(min_week) %>% min(), 
                     weeks_migrating_sthd %>% select(max_week) %>% max())

# juvenile migration happens mostly between weeks 15 and 25





## end of script














#### GRAVEYARD




# #########################################################################################
# #########################################################################################
# ## CPS presence offshore of river mouth
# #########################################################################################
# #########################################################################################
# # CPS data comes from model outputs from latitude 34 to 48 degrees north, -133.9 to -115.5 E, from 2011 through most of 2024 at a ~ 5 day time grid
# 
# # load('./data/_namerica_state_wgs_df.RData')
# # load('./data/_bath_200_utm_df.RData')
# 
# process_cps_data <- function(data){
# 
#   cps_dat <- data %>%
#     slice(-1) %>%
#     select(time, latitude, longitude, anchovy_GAM, sardine_GAM, herring_GAM) %>%
#     mutate(
#       # Organize time variables
#       date = as.Date(substr(time, 1, 10), format="%Y-%m-%d"),
#       year = year(date),
#       month = month(date),
#       doy = yday(date),
#       week = ceiling(doy/7),
# 
#       # Change vars to numeric
#       latitude = as.numeric(latitude),
#       longitude = as.numeric(longitude),
#       anchovy_GAM = as.numeric(anchovy_GAM),
#       sardine_GAM = as.numeric(sardine_GAM),
#       herring_GAM = as.numeric(herring_GAM),
#     )
# 
#   # make weekly grid
#   time_grid <- tibble(
#     date = seq(min(cps_dat$date), max(cps_dat$date), by = "weeks")) %>%
#     mutate(
#       year = year(date),
#       month = month(date),
#       doy = yday(date),
#       week = ceiling(doy/7)
#     ) %>%
#     select(year, week)
# 
# 
#   # join cps_dat to time grid
#   cps_weekly <- left_join(time_grid, cps_dat, by = c("year", "week")) %>%
# 
#     # Aggregate by taking averages over each week
#     group_by(year, week) %>%
#     summarize(
#       anchovy_w = mean(anchovy_GAM, na.rm = TRUE),
#       sardine_w = mean(sardine_GAM, na.rm = TRUE),
#       herring_w = mean(herring_GAM, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
# 
#     # linearly interpolate across weeks with missing CPS probabilities
#     mutate(
#       anchovy_presence = na.approx(anchovy_w, na.rm = FALSE),
#       sardine_presence = na.approx(sardine_w, na.rm = FALSE),
#       herring_presence = na.approx(herring_w, na.rm = FALSE)
#     ) %>%
#     select(year, week, anchovy_presence, sardine_presence, herring_presence)
# 
#   # Aggregate by taking averages over each year
#   cps_yearly <- cps_weekly %>%
#     group_by(year) %>%
#     summarize(
#       anchovy_presence = mean(anchovy_presence, na.rm = TRUE),
#       sardine_presence = mean(sardine_presence, na.rm = TRUE),
#       herring_presence = mean(herring_presence, na.rm = TRUE),
#     )
#   return(list(cps_weekly, cps_yearly))
# }
# 
# 
# # cps_dat_whole_region <- read.csv('./data/fish/CPS_whole_region.csv')
# # cps_dat_SFtoCF <- read.csv('./data/fish/CPS_SDMs_CMtoCF.csv')
# cps_dat_CMtoCF <- read.csv('./data/fish/CPS_SDMs_CMtoCF.csv')
# 
# # cps_whole_region <- process_cps_data(cps_dat_whole_region)
# # cps_SFtoCF <- process_cps_data(cps_dat_SFtoCF)
# cps_CMtoCF <- process_cps_data(cps_dat_CMtoCF)
# 
# # cps_whole_region_weekly <- process_cps_data(cps_dat_whole_region)[1]
# # cps_SFtoCF_weekly <- process_cps_data(cps_dat_SFtoCF)[1]
# # cps_CMtoCF_weekly <- process_cps_data(cps_dat_CMtoCF)[1]
# 
# # save(cps_whole_region_weekly, file = "./data/fish/CPS_cleaned_whole_region.RData")
# # save(cps_SFtoCF_weekly, file = "./data/fish/CPS_cleaned_SFtoCapeFlatt.RData")
# # save(cps_CMtoCF_weekly, file = "./data/fish/CPS_cleaned_CapeMendtoCapeFlatt.RData")
# 
# weekly_df <- left_join(weekly_df, as.data.frame(cps_CMtoCF[1]), by = c("year", "week"))
# yearly_df <- left_join(yearly_df, as.data.frame(cps_CMtoCF[2]), by = c("year"))
# 
# 



# #########################################################################################
# #########################################################################################
# ## Eulachon egg larvae counts
# #########################################################################################
# ########################################################################################


# eul <- read.csv('./data/fish/eulachon_run.csv')
# eul <- as.matrix(eul)[1:53,]
# 
# # add in weekly eulachon egg/larvae count data
# eul_dat <- melt(eul[,-1]) %>%
#   mutate(year = as.numeric(substr(Var2, 2, 5)),
#          week = Var1,
#          eulachon_egglarvae = ifelse(is.na(value), 0, value)) %>%
#   select(year, week, eulachon_egglarvae)
# 
# 
# ## Scale weekly egg/larvae data by annual ssb estimates
# # annual ssb data from https://wdfw.wa.gov/sites/default/files/2025-01/2025-sturgeon-smelt-joint-staff-report.pdf
# eulachon_ssb_data <- data.frame(
#   year = 2011:2024,
#   Weeks_sampled_for_SSB = c(19, 25, 29, 22, 33, 25, 18, 13, 16, 10, 17, 19, 20, 17),
#   annual_ssb_and_harvest = c(3300000, 3200000, 9600000, 16600000, 11400000,
#                                        5100000, 1600000, 400000, 4205000, (4205000+9000000)/2, # linearly interpolated missing 2020 data
#                                        9000000, 18300000, 17000000, 10400000))
# 
# eul_weekly <- eul_dat %>%
#   group_by(year) %>%
#   mutate(
#     # make weekly proportions of egg/larvae quantity, summing to 1 across each year
#     weekly_proportion = eulachon_egglarvae / sum(eulachon_egglarvae, na.rm = T),
#     total_yearly = sum(weekly_proportion) # check for accurate scaling -- Looks good!
#   ) %>%
# 
#   # add in eulachon ssb data and multiply into week proportions
#   left_join(eulachon_ssb_data, by = "year") %>%
#   mutate(
#     eulachon_ssb_est_pounds = weekly_proportion * annual_ssb_and_harvest,
#     eulachon_ssb_est_grams = eulachon_ssb_est_pounds * 453.592,  # pounds to grams
#     eulachon_ssb_est_abund = eulachon_ssb_est_grams / 39.1  # pounds to abundance -- using 39.1 grams per fish
#   )
# 
# # # verify accurate scaling -- looks good!
# # eul_dat %>%
# #   group_by(year) %>%
# #   summarise(
# #     reconstructed_annual_ssb = sum(eulachon_ssb_est_pounds, na.rm = TRUE),
# #     original_annual_ssb = first(annual_ssb_and_harvest),
# #     difference = reconstructed_annual_ssb - original_annual_ssb
# #   )
# 
# weekly_df <- left_join(weekly_df,
#                        eul_weekly %>% select(year, week, eulachon_ssb_est_grams, eulachon_ssb_est_pounds, eulachon_ssb_est_abund, eulachon_egglarvae),
#                        by = c("year", "week"))




# #########################################################################################
# #########################################################################################
# ## Eulachon harvest
# #########################################################################################
# ########################################################################################
# 
# ## Annual Harvest 
# ## Commercial landings of eulachon (in mt) and estimated total number of days the fishery was open in the Columbia River from 1935 to 2009.
# eul <- read.csv('./data/fish/eulachon_catch_with_days_open.csv')
# 
# eul_dat <- eul %>% 
#   rename(
#     eulachon_catch_CR = catch_CR,
#     eulachon_days_open_CR = days_CR, 
#     eulachon_catch_Cow = catch_Cow, 
#     eulachon_days_open_Cow = days_Cow
#   )
# 
# yearly_df <- left_join(yearly_df, eul_dat, by = c("year"))


##  Weekly harvest
# eul <- read.csv('./data/fish/eulachon_harvest.csv')
# eul_df <- eul %>% 
#   select(-c("CPUE", "Pounds.2"))
# 
# eul_harv_dat <- eul_df %>%
#   pivot_longer(
#     cols = starts_with("week"),   
#     names_to = "week",            
#     values_to = "eulachon_harvest" 
#   ) %>%
#   # Extract just the numeric part from week
#   mutate(
#     week = as.numeric(str_replace(week, "week", "")),
#     Year = as.numeric(Year)
#   ) %>%
#   # Rename Year to lowercase for consistency with your weekly_df
#   rename(year = Year) %>%
#   # Optionally filter out NA values
#   filter(!is.na(eulachon_harvest))
# 
# weekly_df <- left_join(weekly_df, eul_harv_dat, by = c("year", "week"))




# 
# #########################################################################################
# #########################################################################################
# ## River covariates
# #########################################################################################
# ########################################################################################
# 
# ## Load Data
# river_dat <- read.csv('./data/env/Bonn_river_vars_1976_2024.csv')
# 
# 
# river_dat <- river_dat %>%
#   ## Remove Feb 29 from non-leap years
#   mutate(
#     date_string = paste(mm.dd, year),
#     is_invalid_feb29 = grepl("Feb-29", date_string) & 
#       !leap_year(as.numeric(str_extract(date_string, "\\d{4}$")))
#   ) %>%
#   # Remove those dates and change Feb-29 to 29-Feb for correct parsing of leap year days
#   filter(!is_invalid_feb29) %>%
#   mutate(date_string = ifelse(mm.dd == "Feb-29", paste("29-Feb", year), date_string)) %>% 
#   # Parse the remaining (valid) dates
#   mutate(
#     date = dmy(date_string),
#     year = year(date),
#     month = month(date),
#     doy = yday(date),
#     week = ceiling(doy/7)
#   ) %>% 
#   select(year, week, parameter, value)
# 
# 
# # Calculate river covariates weekly
# river_weekly <- river_dat %>% 
#   group_by(year, week) %>% 
#   summarize(
#     temp = mean(if_else(parameter == "tempc", value, NA_real_), na.rm = TRUE),
#     flow = mean(if_else(parameter == "outflow", value, NA_real_), na.rm = TRUE),
#     spill = mean(if_else(parameter == "spill", value, NA_real_), na.rm = TRUE),
#     .groups = "drop"
#   ) 
# 
# # Calculate river covariates yearly
# river_yearly <- river_dat %>% 
#   group_by(year) %>% 
#   summarize(
#     temp = mean(if_else(parameter == "tempc", value, NA_real_), na.rm = TRUE),
#     flow = mean(if_else(parameter == "outflow", value, NA_real_), na.rm = TRUE),
#     spill = mean(if_else(parameter == "spill", value, NA_real_), na.rm = TRUE),
#     .groups = "drop"
#   ) 
# 
# weekly_df <- left_join(weekly_df, river_weekly, by = c("year", "week"))
# yearly_df <- left_join(yearly_df, river_yearly, by = c("year"))
# 
# 



