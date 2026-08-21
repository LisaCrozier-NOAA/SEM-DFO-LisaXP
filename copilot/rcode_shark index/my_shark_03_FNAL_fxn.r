#' Calculate Climate-Integrated Shark Predation Index
#'
#' @param N_shark Numeric vector of raw shark abundance/counts.
#' @param T_mt Numeric vector of raw temperatures (e.g., SST in °C) for the metabolic multiplier (Mt).
#' @param T_ot Numeric vector of raw temperatures (e.g., SST in °C) for the spatial overlap multiplier (Ot).
#' @param Q10 Numeric value for metabolic scaling (default = 2).
#' @param overlap_slope Numeric value for the logistic overlap steepness (default = 0.8).
#' @param transform_shark Character string indicating shark transformation: "log1p" (default) or "raw".
#' @param scale_final Logical, whether to return the Z-scored index (TRUE, default) or raw index (FALSE).
#'
#' @return A numeric vector representing the integrated shark index across years.

ishark_fxn <- function(N_shark, 
                                             T_mt,
                                             T_mt_ref=mean(T_mt, na.rm = TRUE), 
                                             T_ot, 
                                             Q10 = 2, 
                                             overlap_slope = 1, 
                                             transform_shark = "log1p", 
                                             scale_final = TRUE) {
  
  # 1. Input checks
  if (length(N_shark) != length(T_mt) || length(N_shark) != length(T_ot)) {
    stop("Input vectors (N_shark, T_mt, T_ot) must be of equal length.")
  }
  
  # 2. Transform Abundance (Keep non-negative)
  if (transform_shark == "log1p") {
    minv <- suppressWarnings(min(N_shark, na.rm = TRUE))
    N_pos <- log1p(N_shark - minv)
  } else if (transform_shark == "raw") {
    N_pos <- N_shark
  } else {
    stop("transform_shark must be either 'log1p' or 'raw'.")
  }
  
  # 3. Calculate Metabolic Multiplier (Mt) using raw T_mt differential
  M_t     <- Q10^((T_mt - Tref_mt) / 10)
  
  # 4. Calculate Spatial Overlap Multiplier (Ot) using Z-scaled T_ot for Logistic curve
  # Safe Z-scaling function
  scale_vector <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }
  
  T_ot_z <- scale_vector(T_ot)
  O_t    <- plogis(overlap_slope * T_ot_z) * 2
  
  # 5. Compute Integrated Predation Index in Positive Space
  I_Shark_raw <- N_pos * M_t * O_t
  
  # 6. Standardize or Return Raw Index
  if (scale_final) {
    I_Shark_out <- scale_vector(I_Shark_raw)
  } else {
    I_Shark_out <- I_Shark_raw
  }
  
  return(I_Shark_out)
}