plot_layout_builder_fxn <-function(model_name, node_map) {
#build_dynamic_layout <- function(model_name, node_map) {
  
  # Safe extraction helper: returns the long name, or "" if pruned
  get_v <- function(role_name) {
    if (role_name %in% names(node_map)) node_map[[role_name]] else ""
  }
  
  # Extract all structural components with your exact text formats
  preyNCC   <- get_v("PreyNCC")
  predNCC   <- get_v("PredNCC")
  preyAK    <- get_v("PreyAK")
  predAK    <- get_v("PredAK")
  growth    <- get_v("Growth")
  abundance <- get_v("Abundance")
  sar       <- get_v("SAR")
  
  # --- Build the 5x5 layout matrix based on the model type ---
  if (grepl("1A", model_name)) {
    mat <- matrix(c(
      preyNCC, ""    , ""   , ""       , predNCC,  
      ""     , growth, ""   , abundance, ""     ,
      ""     , ""    , ""   , ""       , ""     ,
      preyAK , ""    , ""   , ""       , predAK , 
      ""     , ""    , sar  , ""       , ""     
    ), nrow = 5, ncol = 5, byrow = TRUE)
    
  } else if (grepl("1B", model_name)) {
    mat <- matrix(c(
      preyNCC, ""    , ""   , ""       , predNCC,  
      ""     , growth, ""   , abundance, ""     ,
      ""     , ""    , ""   , ""       , ""     ,
      preyAK , ""    , ""   , ""       , predAK , 
      ""     , ""    , sar  , ""       , ""     
    ), nrow = 5, ncol = 5, byrow = TRUE)
    
  } else if (grepl("1C", model_name)) {
    mat <- matrix(c(
      preyNCC, ""    , ""    , ""       , predNCC,  
      ""     , ""    , ""    , abundance, ""     ,
      ""     , ""    , growth, ""       , ""     ,
      preyAK , ""    , ""    , ""       , predAK , 
      ""     , ""    , sar   , ""       , ""     
    ), nrow = 5, ncol = 5, byrow = TRUE)
    
  } else {
    stop("Model type must contain 1A, 1B, or 1C!")
  }
  
  # Return the pure matrix directly—tidySEM loves this format!
  return(mat)
}

test_layout <- build_dynamic_layout("DAG1A_short_topmodel", test_map)
print(test_layout)