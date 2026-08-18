#Best models for sharks
#currently, M_t and O_t are fitted to the same temperature index and they are both almost 1, so they aren't reallly doing anything.
#it would be better if M_t was a local temperature index, and O_t was a large scale spatial overlap index, like HCI.
#right now, it is really just a very subtle weighting of the shark abundance by temperature


suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
})


#Shark index function:
#source("functions/shark_index_fxn.r") -- now named v1

#This is the original function created for my_pred_03_shark_index_fxn.r on Aug 6, 2026
#the next version should take the scaling out of this fxn and 
#take scaled values as input for more stability

shark_index_fxn <- function(data_base, 
                            col_shark, 
                            shark_temp_col, 
                            Q10 = 2, 
                            overlap_form, 
                            overlap_slope) {
  
  # 1. Validate Input Columns
  if (!col_shark %in% names(data_base) || !shark_temp_col %in% names(data_base)) {
    return(NULL)
  }
  
  shark_vec   <- data_base[[col_shark]]
  T_shark_raw <- data_base[[shark_temp_col]]
  
  # 2. Temperature Driver & Metabolic Response
  T_shark    <- safe_scale(T_shark_raw)
  Tref_shark <- mean(T_shark_raw, na.rm = TRUE)
  M_t        <- Q10^((T_shark_raw - Tref_shark) / 10)
  
  # 3. Spatial Overlap Form
  O_t <- if (overlap_form == "constant") {
    rep(1, nrow(data_base))
  } else if (overlap_form == "linear") {
    pmax(0.1, 1 + overlap_slope * T_shark)
  } else if (overlap_form == "logistic") {
    plogis(overlap_slope * T_shark) * 2
  } else {
    return(NULL)
  }
  
  # 4. Construct Integrated Dataframe
  shark_index_dat <- data_base %>%
    mutate(
      # Transformations
      shark_z       = safe_scale(shark_vec),
      shark_z_roll2 = rowMeans(cbind(shark_z, dplyr::lag(shark_z, 1)), na.rm = TRUE),
      
      # Temperature and Ecological Drivers
      T_shark       = T_shark,
      M_t           = M_t,
      O_t           = O_t,
      
      # Integrated Indices
      I_shark_z       = shark_z * M_t * O_t,
      I_shark_z_roll2 = shark_z_roll2 * M_t * O_t
    )
  
  return(shark_index_dat)
}
# -----------------------------------------------------------------------------
# 1. Load & Base Data Preparation 
# -----------------------------------------------------------------------------
# Safe scaling utility function
safe_scale <- function(x) {
  if (all(is.na(x))) return(x)
  as.vector(scale(x))
}




        salmon_dat<-read.csv(file.path("data_Lisa/sem_master_data.csv"))%>% 
          clean_names()  %>% 
          select(year,contains("x07"),contains("x16_sar"),contains("sar_"),contains("x09_dfa_hake"))
        names(salmon_dat)
        
        shark_dat<-read.csv(file.path("data_Lisa/shark_wide.csv"))%>% 
          clean_names() %>%
          select(year,contains("goa"))
        names(shark_dat)
        
        sst_dat<-read.csv(file.path("data_Lisa/goa_prey_clim_raw_trends_avg.csv"),row.names = 1) %>% 
          clean_names() %>%
          select(c(1,19:32))
        names(sst_dat)


data_base <- salmon_dat %>%
  left_join(shark_dat %>% select(year,goa_pacific_sleeper_shark), by = "year") %>%
  left_join(sst_dat %>% select(year,enso_dj,pdo_djf,sst_wgoa_coastwatch_junjulaug, sst_egoa_coastwatch_junjulaug), by = "year") %>%
  dplyr::rename(
    Hake      = x09_dfa_hake_age5plus,
    CPUE_Juv  = x07_dfa_cpue_int_spr_jun_hw,
    SAR_Adult = x16_sar
  )


#shark transformation function
    newdat<-shark_index_fxn(data_base=data_base, 
                          col_shark="goa_pacific_sleeper_shark", 
                          shark_temp_col="enso_dj", 
                          Q10=2, 
                          overlap_form="logistic", 
                          overlap_slope=1)

head(newdat)

#explore plots
  dat=newdat[ , c("year","shark_z", "shark_z_roll2","I_shark_z","I_shark_z_roll2") ]
      n=ncol(dat)-1
      matplot(dat[,1],scale(dat[,-1]),type='b',col=1:n,lty=1:n)
      legend("topright",legend=paste(1:n,names(dat[,-1])),col=1:n,bty='n',cex=0.7)

  dat=newdat[ , c("year","enso_dj", "pdo_djf","M_t","O_t") ]
      n=ncol(dat)-1
      matplot(dat[,1],(dat[,-1]),type='b',col=1:n,lty=1:n)
      legend("top",legend=paste(1:n,names(dat[,-1])),col=1:n,bty='n',cex=0.7)
      
