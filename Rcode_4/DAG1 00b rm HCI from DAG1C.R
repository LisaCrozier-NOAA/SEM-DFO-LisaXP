#Remove HCI from variable list, esp DAB1C

#2. DAG1 Merge aic and importance files from all models ----
path0<- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_06_29_SEM_AKPred/"
path<- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/"

# Define the models and folders
model_names <- c("DAG1C_long", "DAG1C_short")

# --- Function to read and label files ---
load_sem_data <- function(name, file_type) {
  file_path <- paste0(path,name, "/", ifelse(file_type == "results", "SEMresultsByClus.csv", "varImportance.csv"))
  
  read_csv(file_path) %>%
    mutate(
      model_id = name,
      dag_type = str_extract(name, "DAG1[ABC]"), # Extract DAG1A or DAG1B
      length   = str_extract(name, "long|short") # Extract long or short
    )
}



# --- Load everything into two main tables ----
all_performance <- map_df(model_names, ~load_sem_data(.x, "results"))
all_importance  <- map_df(model_names, ~load_sem_data(.x, "importance"))

importance_summary <- all_importance %>%
  # 1. Reshape and sort the initial data
  select(model_id, SEMnode, var, importance) %>%
  pivot_wider(names_from = model_id, values_from = importance) %>%
  arrange(SEMnode, var) %>%
  
  # 2. Add the custom SAR row, setting the 1s directly inside the new row
  # (Using your 'model_names' vector ensures every model column gets a 1 safely)
  add_row(
    SEMnode = "SAR", 
    var = "16.SAR", 
    !!!setNames(rep(1, length(model_names)), model_names)
  ) %>%
  
  # 3. Join the metadata lookup
  left_join(var_lookup_NCC_AK %>% select(var, Lisaname), by = "var") %>%
  
  # 4. Clean, round, and calculate max importances
  mutate(across(DAG1C_long:DAG1C_short, as.numeric)) %>% 
  mutate(across(DAG1C_long:DAG1C_short, ~ round(., 3))) %>%
  rowwise() %>%
  mutate(maxImport = round(max(c_across(DAG1C_long:DAG1C_short), na.rm = TRUE), 2)) %>%
  ungroup() %>%
  
  # 5. Final positioning and sorting
  relocate(SEMnode, Lisaname, var, maxImport, DAG1C_long:DAG1C_short) %>%
  arrange(desc(maxImport))
importance_summary



importance_crop<-importance_summary %>% filter(maxImport>0.1)

importance_crop %>% arrange(SEMnode,desc(DAG1C_long))
importance_crop %>% arrange(desc(DAG1C_long))

importance_crop %>% arrange(desc(DAG1C_long))

print(importance_summary)#,n=Inf)

print(importance_summary %>% filter(maxImport>0.2))


my_performance<-all_performance %>%  select(modNum, model_id,length, dag_type, AIC, CFI, AGFI, p_Chi2,
                                            PreyNCCindNames,       PredNCCindNames, GrowthIndNames, AbundanceIndNames, PreyAKindNames, PredAKindNames) 

topmodel1<-my_performance %>%
  group_by(model_id) %>%
  slice(1) %>% 
  ungroup() %>%
  arrange(length)
topmodel1

print(topmodel1)

topmodels<-my_performance %>%
  group_by(model_id) %>%
  slice(1:100) %>% 
  ungroup() %>%
  arrange(length)

print(topmodels)

write.csv(topmodels,"LisaXP/outputs_4/Doug.topmodels100.csv")


write.csv(importance_summary,paste("LisaXP/outputs_4/Doug","06_29_SEM_AKPred","importance_summary.csv",sep="_"))

