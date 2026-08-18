# Calculate Fixed Integer Weight Hazard Index
# Hazard = -1 * [ 1.0*SSL - 0.3*(SSL*SST) + 2.0*(SSL*Herr) + 1.0*(SSL*Cap) ]
# The output might need to be rescaled for SEM input, but it is currently not returning that


ssl_index_fxn <- function(ssl_scaled, sst_scaled, herring_scaled, capelin_scaled,
                                  scalar_ssl     =  1.0,  # Positive: Base hazard multiplier
                                  scalar_sst     =  0.3,  # Positive: Warm water increases predation
                                  scalar_herring = -2.0,  # Negative: Herring buffers salmon predation
                                  scalar_capelin = -1.0) { # Negative: Capelin buffers salmon predation
  
  # Explicit Salmon Predation Hazard with ssl_scaled factored out
  I_SSL_salmon_hazard = ssl_scaled * ( scalar_ssl 
                                       + scalar_sst * sst_scaled 
                                       + scalar_herring * herring_scaled 
                                       + scalar_capelin * capelin_scaled )
  
  return(I_SSL_salmon_hazard)
}
#signs reversed, ssl_scaled duplicated----------
# ssl_index_fxn<-function (ssl_scaled, sst_scaled, herring_scaled, capelin_scaled,
#                  scalar_ssl = 1,
#                  scalar_sst = 0.3,
#                  scalar_herring = -2,
#                  scalar_capelin = -1) {
#   
#     # Uncalibrated Fixed-Weight Hazard Sum
#     # I_SSL_raw = -1 * ( 1.0 * ssl_scaled 
#     #                    - 0.3 * (ssl_scaled * sst_scaled) 
#     #                    + 2.0 * (ssl_scaled * herring_scaled) 
#     #                    + 1.0 * (ssl_scaled * capelin_scaled) )
# 
#     I_SSL_raw =       scalar_ssl * ssl_scaled 
#                        + scalar_sst * (ssl_scaled * sst_scaled) 
#                        + scalar_herring * (ssl_scaled * herring_scaled) 
#                        + scalar_capelin * (ssl_scaled * capelin_scaled) 
#     
#     # Standardize final index for SEM stability
#     I_SSL_simple = as.vector(scale(I_SSL_raw))
#   
#     return(I_SSL_raw)
# }
# 
# 
