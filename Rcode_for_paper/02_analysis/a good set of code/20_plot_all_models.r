library(tidyverse)
library(tidySEM)

# Paths to the 4 Step 1 sweep results
sweep_paths <- list(
  "DAG1A_long"  = "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1A_long/SEMresultsByClus.csv",
  "DAG1B_long"  = "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1B_long/SEMresultsByClus.csv",
  "DAG1A_short" = "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1A_short/SEMresultsByClus.csv",
  "DAG1B_short" = "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1B_short/SEMresultsByClus.csv"
)

# Extract row 1 (best AIC) from each file
step1_winners <- imap_dfr(sweep_paths, function(path, name) {
  if (!file.exists(path)) {
    warning("File not found: ", path)
    return(tibble())
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  df %>%
    slice(1) %>%
    mutate(model_case = name) %>%
    select(model_case, modNum, AIC, p_Chi2, CFI, SRMR, RMSEA,
           PreyNCCindNames, PredNCCindNames, GrowthIndNames, 
           AbundanceIndNames, PreyAKindNames, PredAKindNames)
})

# Display the winners transposed so every indicator column is easy to read
print(t(step1_winners))


library(tidyverse)

# 1. Load the crosswalk file
crosswalk_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/crosswalk_names.csv"
cw <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

# 2. Winning indicator list from Step 1
winning_inds <- c(
  "04.CompFishNCC_smoothed",
  "02.CompJellyNCC_DFA1",
  "09.PredFishNCC_b_DFA1",
  "06.Cond1NCC_DFA1",
  "07.Cond2NCC_smoothed",
  "13.FishPreyAK_DFA1",
  "12.ZooPreyAK_0_smoothed",
  "12.ZooPreyAK_DFA1",
  "10.PredMammalNCC_0_smoothed",
  "15.PredFishAK_2_smoothed",
  "15.PredFishAK_b_DFA1",
  "15.PredFishAK_c_0_smoothed"
)

# 3. Check crosswalk entries for these winning indicators
cw_sub <- cw %>% 
  filter(DFAname %in% winning_inds | LisaName %in% winning_inds | dfa_cols %in% winning_inds)

print(cw_sub)

# 4. Check if these (or their mapped names) exist in your raw data CSV
data_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/sem_altprey_data_1998_2021_10212026.csv"
if (file.exists(data_path)) {
  df_raw <- read.csv(data_path, nrows = 5)
  cat("\nMatched directly in raw CSV names:\n")
  print(intersect(winning_inds, names(df_raw)))
  cat("\nMatched via crosswalk dfa_cols in raw CSV names:\n")
  print(intersect(cw_sub$DFAname, names(df_raw)))
} else {
  cat("\nData file not found at:", data_path, "\n")
}


library(tidyverse)

# 1. Path to Doug's DFA cluster data
clus_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DFA/clusDataDFA.rds"

if (!file.exists(clus_path)) {
  stop("File not found at: ", clus_path)
}

clusData <- readRDS(clus_path)

# 2. Inspect structure and available columns
cat("--- clusData Structure ---\n")
str(clusData, max.level = 1)
head(clusData, 3)

# 3. Check what latent variables exist in clusData
cat("\n--- Unique SEMlatent values ---\n")
print(unique(clusData$SEMlatent))

# 4. Check for the winning indicator names from Doug's sweep
winning_inds <- c(
  "04.CompFishNCC_smoothed",
  "02.CompJellyNCC_DFA1",
  "09.PredFishNCC_b_DFA1",
  "06.Cond1NCC_DFA1",
  "07.Cond2NCC_smoothed",
  "13.FishPreyAK_DFA1",
  "12.ZooPreyAK_0_smoothed",
  "12.ZooPreyAK_DFA1",
  "10.PredMammalNCC_0_smoothed",
  "15.PredFishAK_2_smoothed",
  "15.PredFishAK_b_DFA1",
  "15.PredFishAK_c_0_smoothed"
)

cat("\n--- Check matches for winning indicators in clusData$shortName ---\n")
print(intersect(winning_inds, unique(clusData$shortName)))

# Also check how many total unique series are in clusData
cat("\nTotal unique shortNames in clusData:", length(unique(clusData$shortName)), "\n")



#run Doug's model------------

library(tidyverse)
library(lubridate)

clus_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DFA/clusDataDFA.rds"
clusData  <- readRDS(clus_path)

# 1. Base filter
clusData_sub <- clusData %>% 
  filter(minYear == 1998) %>%
  mutate(year = year(ymd(date))) %>%
  filter(year >= 1998, year <= 2021)

# 2. Extract winning indicators for DAG1A_long
dag1a_long_inds <- c(
  "PreyNCC"   = "04.CompFishNCC_smoothed",
  "PredNCC"   = "09.PredFishNCC_b_DFA1",
  "Growth"    = "06.Cond1NCC_DFA1",
  "Abundance" = "07.Cond2NCC_smoothed",
  "PreyAK"    = "13.FishPreyAK_DFA1",
  "PredAK"    = "10.PredMammalNCC_0_smoothed"
)

# Helper function to extract, standardize, and rename a single node's series
pull_std_node <- function(cdat, var_node, ind_name) {
  df_node <- cdat %>%
    filter(SEMlatent == var_node, shortName == ind_name) %>%
    select(date, finalVal) %>%
    arrange(date)
  
  # Standardize
  df_node[[var_node]] <- as.vector(scale(df_node$finalVal))
  df_node %>% select(date, all_of(var_node))
}

# 3. Pull each predictor node
df_model_1a_long <- pull_std_node(clusData_sub, "PreyNCC",   dag1a_long_inds["PreyNCC"]) %>%
  full_join(pull_std_node(clusData_sub, "PredNCC",   dag1a_long_inds["PredNCC"]),   by = "date") %>%
  full_join(pull_std_node(clusData_sub, "Growth",    dag1a_long_inds["Growth"]),    by = "date") %>%
  full_join(pull_std_node(clusData_sub, "Abundance", dag1a_long_inds["Abundance"]), by = "date") %>%
  full_join(pull_std_node(clusData_sub, "PreyAK",    dag1a_long_inds["PreyAK"]),    by = "date") %>%
  full_join(pull_std_node(clusData_sub, "PredAK",    dag1a_long_inds["PredAK"]),    by = "date")

# 4. Pull logSAR and standardize
df_sar <- clusData_sub %>%
  filter(SEMlatent == "SAR") %>%
  select(date, finalVal) %>%
  arrange(date) %>%
  mutate(logSAR = as.vector(scale(finalVal))) %>%
  select(date, logSAR)

# 5. Merge and filter complete cases
df_dag1a_long <- df_model_1a_long %>%
  full_join(df_sar, by = "date") %>%
  filter(!is.na(logSAR)) %>%
  arrange(date)

# 6. Verify table dimensions and completeness
cat("Number of observations (rows):", nrow(df_dag1a_long), "\n")
cat("Any missing values remaining?:\n")
print(colSums(is.na(df_dag1a_long)))
head(df_dag1a_long, 3)


library(lavaan)

# 1. Structural model syntax for DAG 1A
mod_dag1a_syntax <- '
  Growth    ~ PreyNCC
  Abundance ~ Growth + PredNCC
  logSAR    ~ Abundance + PreyAK + PredAK
'

# 2. Fit model with lavaan
fit_dag1a_long <- sem(
  model = mod_dag1a_syntax,
  data  = df_dag1a_long
)

# 3. Compare fit measures against Doug sweep result (Target AIC: 115.7936)
fit_summary <- fitMeasures(fit_dag1a_long, c("chisq", "df", "pvalue", "aic", "cfi", "srmr", "rmsea"))
print(round(fit_summary, 4))

# 4. Inspect parameter estimates
pe_1a <- parameterEstimates(fit_dag1a_long, standardized = TRUE)
pe_1a %>% 
  filter(op == "~") %>% 
  select(lhs, op, rhs, est, se, z, pvalue, std.all) %>%
  print()


#Plot Doug's model--------------
library(ggplot2)
library(tidySEM) # or whichever package provides prepare_graph()

# If my_sem_summary is not already in your workspace:
if (!exists("my_sem_summary")) {
  my_sem_summary <- function(fit, bootstrap = FALSE) {
    list(estimates = parameterEstimates(fit, standardized = TRUE))
  }
}

# 1. Custom layout matrix matching DAG 1A topology
# Node order: PreyNCC, PredNCC, Growth, Abundance, PreyAK, PredAK, logSAR
dag1a_layout <- get_layout(
  "PreyNCC", "",       "",          "",       "PreyAK",
  "",        "Growth", "Abundance", "logSAR", "",
  "PredNCC", "",       "",          "",       "PredAK",
  rows = 3
)

# 2. Map labels to clean names
node_labels_1a <- c(
  "PreyNCC"   = "Prey (NCC)\n04.CompFish",
  "PredNCC"   = "Predators (NCC)\n09.PredFish_b (Hake)",
  "Growth"    = "Salmon Growth\n06.Cond1",
  "Abundance" = "Salmon Abund.\n07.Cond2",
  "PreyAK"    = "Prey (AK)\n13.FishPrey",
  "PredAK"    = "Predators (AK)\n10.PredMammal",
  "logSAR"    = "Salmon SAR"
)

# 3. Test graph render
p_test_1a <- Lisa_sem_graph.v2(
  fit               = fit_dag1a_long,
  title             = "DAG 1A Long (24yr) - Top Model",
  layout            = dag1a_layout,
  node_labels       = node_labels_1a,
  show_coefficients = TRUE,
  standardized      = TRUE,
  show_r2           = TRUE,
  display           = TRUE
)

#Replace Doug names with LisaNames--------
library(tidyverse)

# 1. Load the crosswalk file
crosswalk_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/crosswalk_names.csv"
cw <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

# 2. Winning indicator shortNames for DAG1A_long
dag1a_inds <- c(
  "PreyNCC"   = "04.CompFishNCC_smoothed",
  "PredNCC"   = "09.PredFishNCC_b_DFA1",
  "Growth"    = "06.Cond1NCC_DFA1",
  "Abundance" = "07.Cond2NCC_smoothed",
  "PreyAK"    = "13.FishPreyAK_DFA1",
  "PredAK"    = "10.PredMammalNCC_0_smoothed"
)

# 3. Helper to look up the LisaName from the crosswalk
get_lisa_label <- function(node_name, ind_short_name, cw_df) {
  row_match <- cw_df %>% 
    filter(DFAname == ind_short_name | dfa_cols == ind_short_name | LisaName == ind_short_name)
  
  if (nrow(row_match) > 0 && !is.na(row_match$LisaName[1]) && row_match$LisaName[1] != "") {
    clean_ind <- row_match$LisaName[1]
  } else {
    clean_ind <- ind_short_name
  }
  
  # Format with structural node name on top and LisaName below
  paste0(node_name, "\n(", clean_ind, ")")
}

# 4. Generate the mapped node_labels vector automatically
node_labels_lisa_1a <- c(
  "PreyNCC"   = get_lisa_label("Prey (NCC)", dag1a_inds["PreyNCC"], cw),
  "PredNCC"   = get_lisa_label("Predators (NCC)", dag1a_inds["PredNCC"], cw),
  "Growth"    = get_lisa_label("Salmon Growth", dag1a_inds["Growth"], cw),
  "Abundance" = get_lisa_label("Salmon Abund.", dag1a_inds["Abundance"], cw),
  "PreyAK"    = get_lisa_label("Prey (AK)", dag1a_inds["PreyAK"], cw),
  "PredAK"    = get_lisa_label("Predators (AK)", dag1a_inds["PredAK"], cw),
  "logSAR"    = "Salmon SAR"
)

cat("--- Generated Node Labels using LisaNames ---\n")
print(node_labels_lisa_1a)

# 5. Render DAG 1A graph with LisaNames
p_lisa_1a <- Lisa_sem_graph.v2(
  fit               = fit_dag1a_long,
  title             = "DAG 1A Long (24yr) - Top Model (LisaNames)",
  layout            = dag1a_layout,
  node_labels       = node_labels_lisa_1a,
  show_coefficients = TRUE,
  standardized      = TRUE,
  show_r2           = TRUE,
  display           = TRUE
)
#fix box sizes, add pretty names ------
# 1. Concise, high-readability labels tailored to DAG 1A Long
labels_1a_clean <- c(
  "PreyNCC"   = "Comp Fish (NCC)",
  "PredNCC"   = "Hake/Mackerel DFA",
  "Growth"    = "Salmon Growth",
  "Abundance" = "Salmon CPUE",
  "PreyAK"    = "Fish Prey (AK)",
  "PredAK"    = "Harbor Seal (NCC)",
  "logSAR"    = "Salmon SAR"
)

# 2. Render with your exact sizing arguments from your previous script
p_dag1a_sized <- Lisa_sem_graph.v2(
  fit               = fit_dag1a_long,
  title             = "DAG 1A Long (24yr) - Top Model",
  layout            = dag1a_layout,
  node_labels       = labels_1a_clean,
  r2_nodes          = c("Growth", "Abundance", "logSAR"),
  show_coefficients = TRUE,
  standardized      = TRUE,
  show_r2           = TRUE,
  
  # Your proven geometry & text parameters:
  width             = 10,
  height            = 6,
  node_width        = 10,
  node_height       = 1.5,
  rect_width        = 3,
  rect_height       = 0.8,
  ellipses_width    = 2,
  ellipses_height   = 1.2,
  text_size         = 4.5,
  variance_diameter = 0.7,
  curvature         = 10,
  
  display           = TRUE,
  save              = FALSE
)


#resize-------
library(tidySEM)

# 1. Wide layout matrix with blank separator columns ("") between nodes
dag1a_layout_wide <- get_layout(
  "PreyNCC", "", "",       "", "",          "", "",       "", "PreyAK",
  "",        "", "Growth", "", "Abundance", "", "logSAR", "", "",
  "PredNCC", "", "",       "", "",          "", "",       "", "PredAK",
  rows = 3
)

# 2. Readable labels
labels_1a_clean <- c(
  "PreyNCC"   = "Comp Fish\n(NCC)",
  "PredNCC"   = "Hake/Mack.\nDFA",
  "Growth"    = "Salmon\nGrowth",
  "Abundance" = "Salmon\nCPUE",
  "PreyAK"    = "Fish Prey\n(AK)",
  "PredAK"    = "Harbor Seal\n(NCC)",
  "logSAR"    = "Salmon\nSAR"
)

# 3. Render with widened grid spacing
p_dag1a_sized <- Lisa_sem_graph.v2(
  fit               = fit_dag1a_long,
  title             = "DAG 1A Long (24yr) - Top Model",
  layout            = dag1a_layout_wide,
  node_labels       = labels_1a_clean,
  r2_nodes          = c("Growth", "Abundance", "logSAR"),
  show_coefficients = TRUE,
  standardized      = TRUE,
  show_r2           = TRUE,
  
  # Adjusted dimensions: smaller box width gives edges room to breathe
  rect_width        = 1.3,
  rect_height       = 0.6,
  ellipses_width    = 1.3,
  ellipses_height   = 0.6,
  text_size         = 3.5,
  variance_diameter = 0.5,
  curvature         = 10,
  
  width             = 11,
  height            = 5,
  display           = TRUE,
  save              = FALSE
)
 
#Step 5e: Vertical Stacking Layout for DAG 1A----

library(tidySEM)

# Staggered 4-row x 7-column grid layout
dag1a_layout_stacked <- get_layout(
  "PreyNCC", "", "",          "", "",       "", "PreyAK",
  "",        "", "Growth",    "", "",       "", "",
  "",        "", "Abundance", "", "logSAR", "", "",
  "PredNCC", "", "",          "", "",       "", "PredAK",
  rows = 4
)

# Render DAG 1A with stacked intermediate nodes
p_dag1a_stacked <- Lisa_sem_graph.v2(
  fit               = fit_dag1a_long,
  title             = "DAG 1A Long (24yr) - Top Model",
  layout            = dag1a_layout_stacked,
  node_labels       = labels_1a_clean,
  r2_nodes          = c("Growth", "Abundance", "logSAR"),
  show_coefficients = TRUE,
  standardized      = TRUE,
  show_r2           = TRUE,
  
  rect_width        = 1.3,
  rect_height       = 0.6,
  ellipses_width    = 1.3,
  ellipses_height   = 0.6,
  text_size         = 3.5,
  variance_diameter = 0.5,
  curvature         = 10,
  
  width             = 10,
  height            = 6,
  display           = TRUE,
  save              = FALSE
)

#back to lisaNames----
library(tidyverse)
library(tidySEM)

# 1. Load the crosswalk file
crosswalk_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/metadata/crosswalk_names.csv"
cw <- read.csv(crosswalk_path, stringsAsFactors = FALSE)

# 2. Winning indicator shortNames for DAG1A_long
dag1a_inds <- c(
  "PreyNCC"   = "04.CompFishNCC_smoothed",
  "PredNCC"   = "09.PredFishNCC_b_DFA1",
  "Growth"    = "06.Cond1NCC_DFA1",
  "Abundance" = "07.Cond2NCC_smoothed",
  "PreyAK"    = "13.FishPreyAK_DFA1",
  "PredAK"    = "10.PredMammalNCC_0_smoothed"
)

# 3. Print the exact crosswalk entries so you can inspect them directly
cat("--- Exact Crosswalk Entries for Winning Indicators ---\n")
cw_inspect <- cw %>%
  filter(DFAname %in% dag1a_inds | dfa_cols %in% dag1a_inds) %>%
  select(DFAname, LisaName, indicator1, rankedIndicators)
print(cw_inspect)

# 4. Helper to build explicit informational labels:
#    Format: Latent Concept \n [Exact LisaName]
build_explicit_label <- function(node_title, ind_key, cw_df) {
  row_match <- cw_df %>% 
    filter(DFAname == ind_key | dfa_cols == ind_key)
  
  lisa_val <- if (nrow(row_match) > 0 && !is.na(row_match$LisaName[1]) && row_match$LisaName[1] != "") {
    row_match$LisaName[1]
  } else {
    ind_key
  }
  
  paste0(node_title, "\n[", lisa_val, "]")
}

labels_1a_explicit <- c(
  "PreyNCC"   = build_explicit_label("Prey (NCC)", dag1a_inds["PreyNCC"], cw),
  "PredNCC"   = build_explicit_label("Predators (NCC)", dag1a_inds["PredNCC"], cw),
  "Growth"    = build_explicit_label("Salmon Growth", dag1a_inds["Growth"], cw),
  "Abundance" = build_explicit_label("Salmon Abund.", dag1a_inds["Abundance"], cw),
  "PreyAK"    = build_explicit_label("Prey (AK)", dag1a_inds["PreyAK"], cw),
  "PredAK"    = build_explicit_label("Predators (AK)", dag1a_inds["PredAK"], cw),
  "logSAR"    = "Salmon SAR\n[X16_SAR]"
)

cat("\n--- Generated Node Labels ---\n")
print(labels_1a_explicit)

# 5. Expanded Layout Grid:
#    Adding more blank columns between nodes gives room for substantially larger boxes
dag1a_layout_large <- get_layout(
  "PreyNCC", "", "", "", "", "", "", "PreyAK",
  "",        "", "Growth",    "", "", "", "", "",
  "",        "", "Abundance", "", "", "logSAR", "", "",
  "PredNCC", "", "", "", "", "", "", "PredAK",
  rows = 4
)

# 6. Render with substantially larger boxes and clear text
p_dag1a_large <- Lisa_sem_graph.v2(
  fit               = fit_dag1a_long,
  title             = "DAG 1A Long (24yr) - Top Model (Explicit LisaNames)",
  layout            = dag1a_layout_large,
  node_labels       = labels_1a_explicit,
  r2_nodes          = c("Growth", "Abundance", "logSAR"),
  show_coefficients = TRUE,
  standardized      = TRUE,
  show_r2           = TRUE,
  
  # Substantially larger box dimensions
  rect_width        = 3,
  rect_height       = 0.9,
  ellipses_width    = 2.4,
  ellipses_height   = 0.9,
  text_size         = 3.2,
  variance_diameter = 0.6,
  curvature         = 10,
  
  width             = 14,
  height            = 8,
  display           = TRUE,
  save              = FALSE
)

#make insign lines dotted--------
p_dag1a_v3 <- Lisa_sem_graph.v3(
  fit               = fit_dag1a_long,
  title             = "DAG 1A Long (24yr) - Top Model (Dashed = Non-Sig)",
  layout            = dag1a_layout_large,
  node_labels       = labels_1a_explicit,
  r2_nodes          = c("Growth", "Abundance", "logSAR"),
  show_coefficients = TRUE,
  p_threshold       = 0.05,
  nonsig_linetype   = "dashed", # Change to "dotted" if preferred
  standardized      = TRUE,
  show_r2           = TRUE,
  font_family       = "serif",
  rect_width        = 3,
  rect_height       = 0.9,
  ellipses_width    = 2.4,
  ellipses_height   = 0.9,
  text_size         = 3.2,
  variance_diameter = 0.6,
  curvature         = 10,
  
  width             = 14,
  height            = 8,
  display           = TRUE,
  save              = FALSE
)
