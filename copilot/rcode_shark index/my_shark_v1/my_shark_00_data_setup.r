suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
  library(glue)
})

# ---------------------------
# 0) Paths
# ---------------------------
path_sem   <- "data_Lisa/sem_master_data.csv"
path_shark <- "data_Lisa/AKshark.table19.3_2022assess.csv"

path_egoa  <- "Ferris-DFAIndicators-goa/data/EGOA_EcoState_Data_Jan2023.csv"
path_wgoa  <- "Ferris-DFAIndicators-goa/data/WGOA_EcoState_Data_Jan2023.csv"

out_dir    <- "copilot/outputs_2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# data_Lisa/goa_prey_clim_raw_trends_avg.csv for AK data already combined see Lisa.coalate.DFAtrends.r for source file
# data_Lisa/ak_yr for all raw AK data

        # ---------------------------
        # 1) Helper functions for loading data
        # ---------------------------
        guess_year_col <- function(df) {
          nms <- names(df)
          hit <- nms[str_detect(nms, regex("^(year|yr)$|year", ignore_case = TRUE))]
          if (length(hit) == 0) return(NA_character_)
          hit[1]
        }
        
        to_annual_numeric <- function(df, df_name = "data") {
          yc <- guess_year_col(df)
          if (is.na(yc)) stop(df_name, ": could not infer a year column.")
          
          df %>%
            rename(year = all_of(yc)) %>%
            mutate(year = suppressWarnings(as.integer(year))) %>%
            filter(!is.na(year)) %>%
            group_by(year) %>%
            summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
            arrange(year)
        }
        
        safe_read_csv_clean <- function(path, name) {
          if (!file.exists(path)) stop("Missing file: ", path, " (", name, ")")
          read_csv(path, show_col_types = FALSE) %>% clean_names()
        }

# ---------------------------
# 2) Read raw data
# ---------------------------
sem_raw   <- safe_read_csv_clean(path_sem, "SEM master")
shark_raw <- safe_read_csv_clean(path_shark, "Shark table")
egoa_raw  <- safe_read_csv_clean(path_egoa, "EGOA indicators")

#calculate running means before you lose the full time series
library(tidyverse)

shark_processed <- shark_raw %>%
  select(fmp,year,spiny_dogfish,pacific_sleeper_shark,salmon_shark) %>%
  # 1. Ensure data is ordered chronologically within each spatial region
  arrange(fmp, year) %>%
  
  # 2. Group by fmp so averages never cross spatial boundaries
  group_by(fmp) %>%
  
  # 3. Calculate forward (t, t+1) and backward (t, t-1) 2-year averages
  mutate(
    # Forward 2-year average: (Year t + Year t+1) / 2
    across(
      .cols  = spiny_dogfish:salmon_shark, 
      .fns   = ~ (.x + lead(.x, 1)) / 2, 
      .names = "{.col}_2y_avg_forward"
    ),
    # Backward 2-year average: (Year t + Year t-1) / 2
    across(
      .cols  = spiny_dogfish:salmon_shark, 
      .fns   = ~ (.x + lag(.x, 1)) / 2, 
      .names = "{.col}_2y_avg_backward"
    )
  ) %>%
 
  # 4. Clean up the grouping metadata
  ungroup()

shark_wide <- shark_processed %>%
  # 1. Pivot the data columns to wide format using fmp in the column names
  pivot_wider(
    id_cols = year,
    names_from = fmp,
    values_from = c(spiny_dogfish, pacific_sleeper_shark, salmon_shark, 
                    pacific_sleeper_shark_2y_avg_forward, salmon_shark_2y_avg_forward,
                    pacific_sleeper_shark_2y_avg_backward, salmon_shark_2y_avg_backward
    ),
    names_glue = "{fmp}_{.value}"
  ) %>%
  # 2. Keep the table sorted chronologically
  arrange(year)

# View the result
head(shark_wide)
tail(shark_wide[,8:ncol(shark_wide)])
names(shark_wide)

write.csv(shark_wide,"data_Lisa/shark_wide.csv",row.names = FALSE)


# -----------------------------------------------------------------------------
# 1. Extract Top Shark Index & Clean Names
# -----------------------------------------------------------------------------
        out_dir <- "copilot/outputs_4"
        
        #also in data_Lisa/
        ak_dat<-read.csv("copilot/outputs_2/data_all_tested_columns_annual.csv");names(ak_dat)
        ak2<-ak_dat %>% select(year,sst_wgoa_coastwatch_junjulaug,pdo_djf,enso_dj,mid_il_capelin,stka_herr_matbiom)
        
        
        sem_master_data<-read.csv("outputs_4/sem_master_data.csv",row.names=1);head(sem_master_data)
        guild.dfas1<-read.csv("outputs_4/guild.dfa.NCC.AK.csv",row.names=1) %>%
          mutate( sarSR=sem_master_data$sarSR,
                  sarUC=sem_master_data$sarUC) %>%
          inner_join(trend_df_1998_topbio_allwgoa,by="year") %>%
          inner_join(trend_df_1998_herring_egoa,by="year") %>%
          inner_join(ak_dat %>% select(year,sst_wgoa_coastwatch_junjulaug,pdo_djf,enso_dj,mid_il_capelin,stka_herr_matbiom),by="year")
        names(guild.dfas1)
        
        guild.dfasAK<-guild.dfas1 %>% rename(
          X21_sst_wgoa_coastwatch_junjulaug=sst_wgoa_coastwatch_junjulaug,
          X21_pdo_djf=pdo_djf,
          X21_enso_dj=enso_dj,
          
          X13_wgoa_cap.pcod=wgoa_cap.pcod,
          X13_egoa_herring = egoa_herring,
          X13_mid_il_capelin = mid_il_capelin,
          X13_stka_herr_matbiom=stka_herr_matbiom)
        
        write.csv(guild.dfasAK,file.path("data_Lisa/guild.dfasAK.csv"))


# Created in my_shark_03_plot_bestmodels.r, 'all_sims' contains the scenario output from the top model (enso_dj | z_roll2)
      out_dir <- "copilot/outputs_4"
      
      all_sims<-read.csv(file.path(out_dir, "shark_predictions.csv"))
      
      top_shark_index <- all_sims %>%
        filter(Scenario == "enso_dj | z_roll2") %>%
        select(year, x15_shark_enso_roll2 = I_Shark)


# Load, clean, scale guild dataset
      guild.dfas1 <- read.csv("data_Lisa/guild.dfasAK.csv", row.names = NULL) %>% 
        clean_names() %>% 
        select(-x10_harbor_seal_cr_2yr_lead) %>%
        left_join(top_shark_index, by = "year") %>%
        mutate(across(-year, ~ as.vector(scale(.))))

