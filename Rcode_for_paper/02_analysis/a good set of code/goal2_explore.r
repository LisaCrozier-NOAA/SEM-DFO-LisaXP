#Explore results for goal 2

head(indicator_diagnostic_summary)
indicator_diagnostic_summary %>% 
  filter(H_code=="H1a", model_id=="DAG1A_short")
#61 % top models chose forage fish, that was always negative
#39% chose jsoes, and that was always positive -- switch result for forage fish