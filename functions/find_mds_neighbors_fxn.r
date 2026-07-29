# 2. Function to find the nearest neighbors in the 2D MDS space
find_mds_neighbors <- function(target_var, mds_data, n_neighbors = 5) {

  # Get target coordinates
  target_coords <- mds_data %>%
    filter(Lisaname == target_var) %>%
    select(Dim1, Dim2)

  if(nrow(target_coords) == 0) return(NULL)

  # Calculate Euclidean distance to all other variables
  mds_data %>%
    filter(Lisaname != target_var) %>%
    mutate(
      Distance = round(sqrt((Dim1 - target_coords$Dim1)^2 + (Dim2 - target_coords$Dim2)^2),3),
      # NEW: Pull the raw correlation value between the target and this neighbor
      Pairwise_Cor_With_Target = round(cor_matrix_full_years[target_var, Lisaname],2),
      
      # NEW: Classify the direction of that relationship
      Relationship_To_Target = if_else(Pairwise_Cor_With_Target >= 0, "Positive (+)", "Negative (-)")
    ) %>%
    arrange(Distance)  %>%
    mutate(Target_Variable = target_var) %>%
    select(Target_Variable,Lisaname, SEMnode,Relationship_To_Target,Pairwise_Cor_With_Target, Distance, contains("SAR_Correlation")) %>%
    head(n_neighbors)
}
