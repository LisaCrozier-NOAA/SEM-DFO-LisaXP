
#This function was created for my_pred_03_shark_index_fxn.r on Aug 6, 2026
#shark_index_fxn.v2(shark_scaled=-2:2,temp_raw_Mt=2,temp_ref_Mt=0,temp_Ot=2)

#v2 returns just the shark index, doesn't need a database input
shark_index_fxn.v2 <- function(shark_scaled, 
                            temp_raw_Mt,
                            temp_ref_Mt=0, #or mean(temp_raw_Mt,na.rm=T)
                            temp_Ot,
                            Q10 = 2, 
                            overlap_form="logistic", 
                            overlap_slope=0.5) {
  
  # 1. Shark Input
  shark_vec   <- shark_scaled
      shark_z       = shark_vec
      shark_z_roll2 = rowMeans(cbind(shark_z, dplyr::lag(shark_z, 1)), na.rm = TRUE)
  
  # 2. Temperature Driver & Metabolic Response -- must be raw temperatures for metabolic fxn to work
  T_shark_raw <- temp_raw_Mt
  Tref_shark <- temp_ref_Mt
#  Tref_shark <- mean(T_shark_raw, na.rm = TRUE)
  M_t        <- Q10^((T_shark_raw - Tref_shark) / 10)

  # 3. Spatial Overlap Form
  T_shark    <- temp_Ot
  O_t <- if (overlap_form == "constant") {
          rep(1, length(temp_Ot))
        } else if (overlap_form == "linear") {
          pmax(0.1, 1 + overlap_slope * T_shark)
        } else if (overlap_form == "logistic") {
          plogis(overlap_slope * T_shark) * 2
        } else {
          return(NULL)
        }
  
  # 4. Integrated Indices
      I_shark_z       = shark_z * M_t * O_t
      I_shark_z_roll2 = shark_z_roll2 * M_t * O_t
    
  
  return(I_shark_z)
}



#original, returns a full database-----------
shark_index_fxn_v1 <- function(data_base, 
                            col_shark, 
                            shark_temp_col, 
                            Q10 = 2, 
                            overlap_form="logistic", 
                            overlap_slope=1) {
  
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