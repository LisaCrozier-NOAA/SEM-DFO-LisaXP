#' Compute Shark Predation Hazard Index for Salmon
#' 
#' @param shark Vector of shark abundance counts (will be standardized)
#' @param temp_raw_Mt Vector of RAW water temperature in °C (NOT scaled)
#' @param temp_Ot Vector of spatial overlap climate index (e.g., PDO, ENSO)
#' @param alpha_ot Scalar weight for spatial overlap term (default: +0.37)
#' @param Q10 Thermal sensitivity coefficient (default: 2.0)
#' @param overlap_form Type of overlap function: "logistic", "linear", or "constant"
#' @return Vector of standardized (z-score) I_Shark hazard values

calc_i_shark <- function(shark, temp_raw_Mt, temp_Ot,
                         alpha_ot     = 0.372,
                         Q10          = 2.0,
                         overlap_form = "logistic") {
  
  # Helper to safely scale vectors (z-score)
  safe_z <- function(x) {
    if (all(is.na(x)) || var(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
    as.vector(scale(x))
  }
  
  # 1. Standardize Shark Abundance
  shark_z <- safe_z(shark)
  
  # 2. Metabolic Multiplier (Mt) using raw Celsius temperatures
  temp_ref <- mean(temp_raw_Mt, na.rm = TRUE)
  M_t      <- Q10^((temp_raw_Mt - temp_ref) / 10)
  
  # 3. Spatial Overlap Multiplier (Ot)
  if (overlap_form == "logistic") {
    # Bounded between 0 and 2, centered at 1.0
    O_t <- plogis(alpha_ot * temp_Ot) * 2
  } else if (overlap_form == "linear") {
    # Unbounded linear slope
    O_t <- pmax(0.1, 1.0 + alpha_ot * temp_Ot)
  } else {
    # Constant overlap
    O_t <- rep(1.0, length(temp_Ot))
  }
  
  # 4. Compound Hazard Index
  raw_index <- shark_z * M_t * O_t
  
  # Return standardized hazard vector
  return(safe_z(raw_index))
}


#original scaling factors-------------
#' Compute Shark Predation Hazard Index for Salmon
#'
#' @param shark Vector of raw shark abundance counts
#' @param temp_raw_Mt Vector of RAW water temperature in °C (NOT scaled)
#' @param temp_Ot Vector of spatial overlap climate index (e.g., PDO, ENSO)
#' @param alpha_ot Scalar weight for spatial overlap term (default: +0.372)
#' @param Q10 Thermal sensitivity coefficient (default: 2.0)
#' @param overlap_form Type of overlap function: "logistic", "linear", or "constant"
#' @param shark_ref Reference list containing mean and sd from training fit
#' @param temp_ref_Mt Baseline reference temperature °C for Mt (default: 11.24)
#' @return Vector of I_Shark hazard values on the benchmark standard scale

calc_i_shark_orig_scale <- function(shark, temp_raw_Mt, temp_Ot,
                         alpha_ot     = 0.372,
                         Q10          = 2.0,
                         overlap_form = "logistic",
                         # Default reference parameters from original fit
                         shark_ref    = list(mean = 142.5, sd = 38.2),
                         temp_ref_Mt  = 11.24) {
  
  # 1. Standardize Shark Abundance using FIXED reference parameters
  shark_z <- (shark - shark_ref$mean) / shark_ref$sd
  
  # 2. Metabolic Multiplier (Mt) anchored to the FIXED reference baseline temp
  M_t <- Q10^((temp_raw_Mt - temp_ref_Mt) / 10)
  
  # 3. Spatial Overlap Multiplier (Ot)
  if (overlap_form == "logistic") {
    O_t <- plogis(alpha_ot * temp_Ot) * 2
  } else if (overlap_form == "linear") {
    O_t <- pmax(0.1, 1.0 + alpha_ot * temp_Ot)
  } else {
    O_t <- rep(1.0, length(temp_Ot))
  }
  
  # 4. Compound Hazard Index
  i_shark_hazard <- shark_z * M_t * O_t
  
  return(i_shark_hazard)
}