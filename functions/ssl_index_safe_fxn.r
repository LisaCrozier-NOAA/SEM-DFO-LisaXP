#' Compute Sea Lion Predation Hazard Index for Salmon
#' 
#' @param ssl Vector of sea lion count/pup indices (will be standardized)
#' @param sst Vector of SST values (will be standardized)
#' @param herring Vector of herring biomass values (will be standardized)
#' @param capelin Vector of capelin biomass values (will be standardized)
#' @param beta_sst Scalar weight for SST interaction (default: +0.25)
#' @param beta_herr Scalar weight for Herring interaction (default: -2.0)
#' @param beta_cap Scalar weight for Capelin interaction (default: -1.0)
#' @return Vector of standardized (z-score) I_SSL hazard values

calc_i_ssl <- function(ssl, sst, herring, capelin, 
                       beta_sst  =  0.25, 
                       beta_herr = -2.00, 
                       beta_cap  = -1.00) {
  
  # Helper to safely scale vectors (z-score)
  safe_z <- function(x) {
    if (all(is.na(x)) || var(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
    as.vector(scale(x))
  }
  
  # Standardize inputs
  ssl_z  <- safe_z(ssl)
  sst_z  <- safe_z(sst)
  herr_z <- safe_z(herring)
  cap_z  <- safe_z(capelin)
  
  # Factored Hazard Equation: SSL * (1 + beta_sst*SST + beta_herr*Herr + beta_cap*Cap)
  raw_index <- ssl_z * (1.0 + beta_sst * sst_z + beta_herr * herr_z + beta_cap * cap_z)
  
  # Return standardized hazard vector
  return(safe_z(raw_index))
}


#Original scaling factors--------------
#' Compute Sea Lion Predation Hazard Index for Salmon
#'
#' @param ssl Vector of raw sea lion count/pup values
#' @param sst Vector of raw SST values
#' @param herring Vector of raw herring biomass values
#' @param capelin Vector of raw capelin biomass values
#' @param beta_sst Weight for SST interaction (default: +0.25)
#' @param beta_herr Weight for Herring interaction (default: -2.0)
#' @param beta_cap Weight for Capelin interaction (default: -1.0)
#' @param ssl_ref Reference list containing mean and sd from training fit
#' @param sst_ref Reference list containing mean and sd from training fit
#' @param herr_ref Reference list containing mean and sd from training fit
#' @param cap_ref Reference list containing mean and sd from training fit
#' @return Vector of I_SSL hazard values on the benchmark standard scale

calc_i_ssl_orig_scale <- function(ssl, sst, herring, capelin,
                       beta_sst  =  0.25,
                       beta_herr = -2.00,
                       beta_cap  = -1.00,
                       # Default reference scaling parameters from original fit
                       ssl_ref  = list(mean = 27850, sd = 4120),
                       sst_ref  = list(mean = 11.24, sd = 0.85),
                       herr_ref = list(mean = 145000, sd = 32000),
                       cap_ref  = list(mean = 12.4, sd = 4.1)) {
  
  # Standardize using the FIXED reference parameters
  ssl_z  <- (ssl - ssl_ref$mean) / ssl_ref$sd
  sst_z  <- (sst - sst_ref$mean) / sst_ref$sd
  herr_z <- (herring - herr_ref$mean) / herr_ref$sd
  cap_z  <- (capelin - cap_ref$mean) / cap_ref$sd
  
  # Factored Hazard Equation
  i_ssl_hazard <- ssl_z * (1.0 + beta_sst * sst_z + beta_herr * herr_z + beta_cap * cap_z)
  
  return(i_ssl_hazard)
}