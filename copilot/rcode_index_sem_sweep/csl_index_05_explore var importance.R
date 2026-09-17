


out_dir <- "copilot/outputs_6"
all_models<- read.csv(file.path(out_dir, "sem_global_2var_3var_models_sorted_aic.csv"), row.names = NULL)
write.csv(all_models_r2,   file.path(out_dir, "sem_global_2var_3var_models_sorted_r2.csv"), row.names = FALSE)
write.csv(var_importance,  file.path(out_dir, "sem_global_variable_importance_aicweights.csv"), row.names = FALSE)


#remove csl/chinook ratio index

all_models_aic[1,mycol]
names(all_models_aic)
mycol=c("model_type","terms","aic","sar_r2")

# ---------------------------
# 7) AIC weights + variable importance
# ---------------------------

all_models_sub <- all_models_aic %>%
 # filter(!grepl("x15_icsl_cslchin_ratio",terms)) %>%
 # filter(!grepl("x21",terms)) %>%
  # filter(model_type=="2var_global") 
  filter(model_type=="3var_global") 

all_models_sub[1:5,]  

all_models_sub[1:5,mycol]
all_models_sub_3var<-all_models_sub
all_models_sub_2var<-all_models_sub
all_models_sub_2var_include_cslchin_ratio<-all_models_sub

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
    title = "Variable importance from global 2- and 3-term SEM sweep",
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

importance_plot_df
importance_plot_df %>%
filter(var_group=="Index")

p_imp1_noratio_noclimate_3var


p_imp1_noratio_noclimate_2var
# variable                        importance n_models in_any_index_model_weight var_group
# 1 x15_issl_z                         0.487         34                   0.487   Index    
# 2 x15_ishark_df_asleeper_z           0.141         34                   0.141   Index    
# 3 x15_ishark_z                       0.0211        34                   0.0211  Index    
# 4 x15_icsl_risk_eulachon75_shad25    0.00874       34                   0.00874 Index    

p_imp1_noratio_noclimate_3var<-p_imp1
#p_imp1_noratio_noclimate_2var<-p_imp1
#p_imp1_noratio<-p_imp1


ggsave(
  file.path(out_dir, "sem_global_variable_importance_aicweights.png"),
  p_imp, width = 10, height = 12, dpi = 170
)

ggsave(
  file.path(out_dir, "sem_2var_noratio_noclimate_importance.png"),
  p_imp1_noratio_noclimate_2var, width = 10, height = 12, dpi = 170
)
ggsave(
  file.path(out_dir, "sem_3var_noratio_noclimate_importance.png"),
  p_imp1_noratio_noclimate_3var, width = 10, height = 12, dpi = 170
)
