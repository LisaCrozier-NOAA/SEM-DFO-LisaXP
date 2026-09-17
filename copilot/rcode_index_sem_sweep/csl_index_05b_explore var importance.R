


out_dir <- "copilot/outputs_7"
all_models<- read.csv(file.path(out_dir, "sem_global_2var_models_full.csv"), row.names = NULL)


#remove csl/chinook ratio index

mycol=c("model_type","terms","aic","sar_r2")
mycol_b=c("model_type","terms","aic","sar_r2","b1","p1","b2","p2","b3","p3")
mycol_b_only=c("b1","p1","b2","p2","b3","p3")
all_models_aic[1:5,mycol]
round(all_models_aic[1:5,mycol_b_only],3)
# ---------------------------
# 7) AIC weights + variable importance
# ---------------------------
mytitle<-"2 & 3 variables, all AK drivers"
plotname<-"alldrivers"
      all_models_sub <- all_models #%>%
        # filter(!grepl("x15_icsl_cslchin_ratio",terms)) %>%
        # filter(!grepl("x21",terms)) %>%
       # filter(model_type=="2var_global") 
      #filter(model_type=="3var_global") 


mytitle<-"2 variables, all AK drivers"
plotname<-"alldrivers_2var"

all_models_sub <- all_models %>%
 # filter(!grepl("x15_icsl_cslchin_ratio",terms)) %>%
 # filter(!grepl("x21",terms)) %>%
   filter(model_type=="2var_global") 
  #filter(model_type=="3var_global") 

mytitle<-"2 variables, no climate"
plotname<-"noclimate_2var"

all_models_sub <- all_models %>%
  # filter(!grepl("x15_icsl_cslchin_ratio",terms)) %>%
   filter(!grepl("x21",terms)) %>%
  filter(model_type=="2var_global") 
#filter(model_type=="3var_global") 


mytitle<-"2 variables, no climate, no ratio"
plotname<-"noclimate_noratio_2var"

all_models_sub <- all_models %>%
   filter(!grepl("x15_icsl_cslchin_ratio",terms)) %>%
  filter(!grepl("x21",terms)) %>%
  filter(model_type=="2var_global") 
#filter(model_type=="3var_global") 


all_models_sub[1:5,mycol]


assign(paste0("topmodels_",plotname),all_models_sub[1:5,],envir=.GlobalEnv)
print(paste0("topmodels_",plotname))

#run the rest----------
  all_models_sub <- all_models_sub %>%
  mutate(
    delta_aic_global = aic - min(aic, na.rm = TRUE),
    aic_weight_raw = exp(-0.5 * delta_aic_global),
    aic_weight = aic_weight_raw / sum(aic_weight_raw, na.rm = TRUE)
  ) %>%
  select(-aic_weight_raw)

# long model-term table for importance -- remove index
model_terms_long <- all_models_sub %>%
  select(model_label, model_type, n_terms, terms, aic, aic_weight, term1, term2, term3,
         has_ishark, has_issl, has_icsl, has_any_index, index_count) %>%
  pivot_longer(cols = c(term1, term2, term3), names_to = "slot", values_to = "variable") %>%
  filter(!is.na(variable))

var_importance <- model_terms_long %>%
  group_by(variable) %>%
  summarise(
    importance = sum(aic_weight, na.rm = TRUE),
    n_models = n(),
    in_any_index_model_weight = sum(aic_weight[has_any_index], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(importance))

# optional variable tags
var_importance <- var_importance %>%
  mutate(
    var_group = case_when(
      str_detect(variable, "ishark|issl|icsl") ~ "Index",
      str_detect(variable, "^x10_|^x15_") ~ "Predator",
      str_detect(variable, "^x21_") ~ "Climate",
      str_detect(variable, "^x12_|^x13_") ~ "Prey",
      str_detect(variable, "^x14_") ~ "Competitor",
      TRUE ~ "Other"
    )
  )


# ---------------------------
# 9) Variable importance plot
# ---------------------------
importance_plot_df <- var_importance %>%
  mutate(variable = fct_reorder(variable, importance))

p_imp1 <- ggplot(importance_plot_df, aes(x = variable, y = importance, fill = var_group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~var_group, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(
    title = mytitle,
    subtitle = "Importance = summed AIC weights across models containing variable",
    x = "Variable",
    y = "AIC-weighted importance",
    fill = NULL
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_imp1


assign(paste0("plot_importance_",plotname),p_imp,envir=.GlobalEnv)
print(paste0("plot_importance_",plotname))


assign(paste0("table_importance_",plotname),importance_plot_df,envir=.GlobalEnv)
print(paste0("table_importance_",plotname))

importance_plot_df
importance_plot_df %>%
filter(var_group=="Index")



ggsave(
  file.path(out_dir, paste0("sem_var_importance_",plotname,".png")),
  p_imp1, width = 10, height = 12, dpi = 170
)
#print results---------
print(paste0("sem_var_importance_",plotname,".png"))
print(paste0("topmodels_",plotname))
print(all_models_sub[1:5,mycol])
print(round(all_models_sub[1:5,mycol_b_only],3))
print(paste0("table_importance_",plotname))
print(importance_plot_df %>%  filter(importance>0.1))

