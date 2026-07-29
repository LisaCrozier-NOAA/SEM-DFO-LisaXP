#Dump, this was my example that I fed into gemini

#take home here is that you have to have hakeage5Plus in the model
#everything added to that seems kind of spurious to me, but maybe we could make a latent variable that improves on hake
#or maybe we can figure out what was driving hake, and use that as a different route
#should check whether these are pos or neg coef -- direct or indirect mechanisms?
  


# Define your roles
response_var <- "X16_SAR"

candidates <- var_lookup_NCC_AK %>% 
  filter(!node_id %in% c("X06","X07","X16")) %>% 
  pull(Lisaname)

#candidates <- setdiff(var$Lisaname, c(response_var, fixed_predictor))

# Build the Global Formula
full_formula <- as.formula(
  paste(response_var, "~", paste(c( candidates), collapse = " + "))
)


# # Calculate absolute correlation matrix for all predictors
cor_matrix <- abs(cor(guild.dfas1[, candidates], use = "complete.obs"))

# Mark pairs with correlation > 0.5
exclude_pairs <- cor_matrix > 0.5
diag(exclude_pairs) <- FALSE  # Don't exclude a variable for being correlated with itself



# Fit the global model
# Reminder: na.fail is required for dredge to work!
global_model <- lm(full_formula, data = guild.dfas1, na.action = "na.fail")

# Execute dredge
dredge_results <- dredge(
  global_model,
  rank = "AICc",
  fixed = c(fixed_predictor),
  m.lim = c(1, 4)#,            # Total variables in model (including fixed)
#  subset = !exclude_pairs     # Logical matrix to prevent colinearity
)

# View top models---------
head(dredge_results)
#write.csv(as.data.frame(dredge_results),file="LisaXP/NCC.dredge_results.csv")
dredge_results<-read.csv("LisaXP/outputs_4/NCC_AK_dredge_results.csv",row.names = 1)


# 1. Calculate importance (Sum of Weights)-----------
# This returns a named vector of the relative importance of each variable
imp_vector <- sw(dredge_results) 

# 2. Convert to a clean, "pretty" data frame
importance_table <- data.frame(
  Variable = names(imp_vector),
  Importance = as.numeric(imp_vector),
  N_Models = sw(dredge_results) %>% attr("n.models") # Optional: how many models it appeared in
) %>%
  mutate(Importance = round(Importance, 2)) %>% # Rounding step
  arrange(desc(Importance))

# 3. View the table
print(importance_table)

# 4. Save it
write.csv(importance_table, "LisaXP/outputs_4/LisaDredge_NCC_AK_variable_importance.csv", row.names = FALSE)

