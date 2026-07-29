library(tidyverse)

#To prove that the "stragglers" (like Common Murres and Harbor Seals) dominate the short models with mathematically strong but ecologically backward (positive) effects, we need to look under the hood of the model parameter estimates.
#We can write a script that pulls the coefficient estimates, sign direction, and significance for every environmental indicator across all 6 models. 
#By separating the results of the Long Models from the Short Models , we will generate a clean, publication-ready table that clearly contrasts the true ecological signals against the spurious stragglers.

#Script: Quantifying Predictor Performance & Sign Across Models
#This script joins your estimates_all parameter table with your top_100_keys and the var_lookup_NCC_AK metadata. 
#It calculates:
#Selection Rate: How often the variable was picked.
#Sign Direction: What percent of the time the coefficient was Positive ($+$) vs. Negative ($-$).
#Significance Rate: How often it had a robust effect ($p < 0.05$) within those selections.


# Process DAG1A-B coefficients-------- 
          #--- 1. Re-ensure rankings and parameter estimates are loaded ---
          # (Using the objects and paths you defined in your session)
          path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/"
          model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short", "DAG1C_long", "DAG1C_short")
          
#Get nearest neighbor info
          #source("LisaXP/DAG1 00c PCA and MDS plots.R")
          plot_data_mds <- plot_data_mds_guild.dfas1_X16SARcor<- read.csv("LisaXP/outputs_4/mds_coordinates__guild.dfas1_data.csv",row.names = NULL)
          # 1. Prepare the full 24-year dataset (dropping only 'year' and the final column)
          full_data <- guild.dfas1 %>%
            select(-ncol(.)) %>% # Drop the final column
            select(-year)        # Drop the timeline tracker
          
          # 2. Compute a pairwise complete correlation matrix across all 24 years
          # This handles the NAs by pairing whatever years are available for each combo
          cor_matrix_full_years <- cor(full_data, use = "pairwise.complete.obs")
          
#Find top 100 models
          cat("Loading rankings and estimates...\n")
          rankings_all <- map_dfr(model_names, function(name) {
            file_path <- paste0(path, name, "/SEMresultsByClus.csv")
            if (!file.exists(file_path)) return(tibble())
            read_csv(file_path, show_col_types = FALSE) %>% mutate(model_id = name)
          })
          
          top_100_keys <- rankings_all %>%
            group_by(model_id) %>%
            slice_min(order_by = AIC, n = 100, with_ties = FALSE) %>%
            select(model_id, modNum) %>%
            ungroup()
          
          cat("Loading parameter estimates...\n")
          estimates_all <- map_dfr(model_names, function(name) {
            file_path <- paste0(path, name, "/parameterEstimates.csv")
            if (!file.exists(file_path)) return(tibble())
            read_csv(file_path, show_col_types = FALSE) %>% mutate(model_id = name)
          })
          
          
          # --- 2. Map Generic Nodes to Actual Indicators in the Rankings ---
          # We pivot the indicator columns to create a clean "lookup" for each modNum
          rankings_mapped <- rankings_all %>%
            select(
              model_id, modNum, AIC,
              PreyNCCindNames, PredNCCindNames, GrowthIndNames, 
              AbundanceIndNames, PreyAKindNames, PredAKindNames
            ) %>%
            pivot_longer(
              cols = ends_with("indNames"),
              names_to = "generic_node",
              values_to = "actual_indicator"
            ) %>%
            # Clean up the node names (e.g., "PreyNCCindNames" becomes "PreyNCC")
            mutate(generic_node = str_remove(generic_node, "(?i)indnames"))
          
          # --- 3. Join Estimates with Rankings Mapping ---
          # We join on the model_id, modNum, and match the predictor (rhs) to the generic_node
          mapped_estimates <- estimates_all %>%
            filter(op == "~") %>% # Only look at regression paths
            inner_join(rankings_mapped, by = c("model_id", "modNum", "rhs" = "generic_node")) %>%
            # Join with your var_lookup table to bring in Lisanames (the clean biology names)
            left_join(var_lookup_NCC_AK %>% select(var, Lisaname), by = c("actual_indicator" = "var")) %>%
            # Calculate metadata tags for each run
            mutate(
              is_significant = if_else(pvalue < 0.05, 1, 0),
              is_positive = if_else(est > 0, 1, 0),
              is_top_100 = if_else(paste(model_id, modNum) %in% paste(top_100_keys$model_id, top_100_keys$modNum), TRUE, FALSE)
            )
          
          # --- 4. Summarize Indicator Selection, Sign, and Significance for Top 100 ---
          top100_indicator_behavior <- mapped_estimates %>%
            filter(is_top_100 == TRUE) %>%
            group_by(model_id, rhs, Lisaname) %>% # 'rhs' holds the SEM node (e.g. PreyNCC), 'Lisaname' is the variable
            summarise(
              Times_Selected_In_Top100 = n(),
              
              # Sign Behavior
              Times_Positive = sum(is_positive, na.rm = TRUE),
              Pct_Positive = round((Times_Positive / Times_Selected_In_Top100) * 100, 1),
              Pct_Negative = 100 - Pct_Positive,
              
              # Statistical Significance
              Times_Significant = sum(is_significant, na.rm = TRUE),
              Pct_Significant = round((Times_Significant / Times_Selected_In_Top100) * 100, 1),
              
              Mean_Estimate = round(mean(est, na.rm = TRUE), 3),
              .groups = "drop"
            ) %>%
            # Sort to put the most dominant indicators per model/node at the top
            arrange(model_id, rhs, desc(Times_Selected_In_Top100))
          
          # Save top100_indicator_behavior file-----
          #--- 5. Save results ---
          print(top100_indicator_behavior, n = Inf)
          write_csv(top100_indicator_behavior, "LisaXP/outputs_4/Top_100_Indicator_Sign_and_Significance_ABC.csv")

#Add explanations===========
          #To turn your detailed ecological diagnostics into a clean, publication-ready summary, we can write an R script that automates the filters, applies your trophic hypotheses, determines whether a variable is a DFA or a straggler, programmatically pulls the key nearest neighbors, and applies your custom explanations.
          #The Summary Construction Pipeline
          #This script runs the entire analysis in one clean pipe. 
          #It dynamically uses your find_mds_neighbors function behind the scenes to extract the qualifying neighbors ($r > 0.5$ or $r < -0.2$) for every target variable and collapses them into a single scannable cell.
          
          
          library(tidyverse)
          
          # --- 1. Helper Function to Extract Neighbors by Correlation Direction ---
          get_neighbor_summary <- function(target_var, mds_data, direction = "pos") {
            neighbors <- tryCatch({
              find_mds_neighbors(target_var, mds_data, n_neighbors = 10)
            }, error = function(e) NULL)
            
            if (is.null(neighbors) || nrow(neighbors) == 0) return("")
            
            # Filter based on the requested direction and thresholds
            if (direction == "pos") {
              filtered <- neighbors %>% filter(Pairwise_Cor_With_Target > 0.5)
            } else {
              filtered <- neighbors %>% filter(Pairwise_Cor_With_Target < -0.2)
            }
            
            if (nrow(filtered) == 0) return("")
            
            # Format as "Name (Value)"
            filtered %>%
              mutate(Formatted = sprintf("%s (%s)", Lisaname, round(Pairwise_Cor_With_Target, 2))) %>%
              pull(Formatted) %>%
              paste(collapse = ", ")
          }
          
          # --- 2. Build the Ecological Synthesis Table ---
          ecological_synthesis <- top100_indicator_behavior %>%
            # Filter 1: Keep only indicators selected more than 10 times in the top 100 models
            filter(Times_Selected_In_Top100 > 10) %>%
            rename(Times_In_Top100 = Times_Selected_In_Top100) %>%
            
            # Filter 2: Exclude DAG1C for mediators & NCC nodes
            filter(!(rhs %in% c("Abundance", "Growth", "PreyNCC", "PredNCC") & str_detect(model_id, "DAG1C"))) %>%
            
            # Add columns for Hypothesis, Type (DFA vs. Straggler), and Confirmation
            mutate(
              Node2 = case_when(
                grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", Lisaname) ~ "Competitor",
                TRUE ~ as.character(rhs) # Force character type consistency
              ),
              # Define expected ecological relationship
              Hyp = case_when(
                grepl("X04|X05|X14|X13_DFA_WGOA_DFA_midTrophic", Lisaname) ~ "Negative",
                rhs %in% c("Abundance", "Growth", "PreyNCC", "PreyAK") ~ "Positive",
                rhs %in% c("PredNCC", "PredAK") ~ "Negative",
                TRUE                            ~ "Neutral"
              ),
              
              Type = if_else(str_detect(Lisaname, "DFA"), "DFA", "Straggler"),
              
              OK = case_when(
                Hyp == "Positive" & Pct_Positive > 50 ~ 1,
                Hyp == "Negative" & Pct_Negative > 50 ~ 1,
                TRUE                                  ~ 0
              )
            ) %>%
            
            # Fetch neighbors programmatically for both directions
            rowwise() %>%
            mutate(
              Pos_Neigh = get_neighbor_summary(Lisaname, plot_data_mds, direction = "pos"),
              Neg_Neigh = get_neighbor_summary(Lisaname, plot_data_mds, direction = "neg")
            ) %>%
            ungroup() %>%
            
            # Add Your Custom Explanations
            mutate(
              Explanation = case_when(
                Lisaname == "X06_DFA_IGF_mu" ~ "Positive & significant (1B, 1C) or negative but NS (1A)",
                # Core Salmon Diagnostics (Successes)
                rhs == "Abundance" & OK == 1 ~ "Better NCC survival carries through AK stage",
                rhs == "Growth"    & OK == 1 ~ "Bigger is better",
                
                # Core Salmon Diagnostics (Discrepancies)
                rhs == "Abundance" & OK == 0 ~ "Discrepancy: Survival signal decoupled",
                rhs == "Growth"    & OK == 0 ~ "Discrepancy: Growth signal decoupled",

                
                # PreyNCC explanations
                #Confirmations
                Lisaname == "X05_DFA_abundSardine" ~ "Competitor (tracking herring / fewer sardines = more food for salmon)",
                Lisaname == "X01_DFA_sumPreyOfPrey_planktonJun" ~ "Direct trophic support (more food = more growth)",
                
                # PredNCC explanations
                #Discrepancies
                Lisaname == "X08_commonMurre_JSOES" ~ "Ecosystem proxy (Murre ~ Habitat Compression Index)",
                Lisaname == "X08_Loons_8_WS" ~ "Coincidence (Loon population co-trending with CPUE)",
                
                #Confirmations
                Lisaname == "X09_DFA_HakeAge5Plus" ~ "Direct predator",
                
                # PredAK explanations
                #Discrepancies
                Lisaname == "X10_Harbor_seal_CR_2yrLead" ~ "Spurious (Noisy / high missing data; trend-rider on short series)",
                Lisaname == "X15_PacificCodBiomass_predAK" ~ "Ecosystem proxy (Pcod had negative response to heat wave, favored other predators)",
                Lisaname == "X15_halibutBiomassAge8plus_2yrLead_predAK" ~ "Ecosystem proxy (Halibut declines in warm water favor competitors)",
                
                #Confirmations
                Lisaname == "X10_Harbour_s_2yrLead_WS" ~ "Direct predator",
                Lisaname == "X10_DFA_ssl.est.wholerange_2yrLead" ~ "Direct predator (deepwater top-down consumer)",
                Lisaname == "X15_DFA_sleeperSharkBSAI_predAK" ~ "Direct predator (deepwater top-down consumer)",
                Lisaname == "X15_ArrowtoothFlounderBiomass_predAK" ~ "Direct predator (deepwater top-down consumer)",
                
                # PreyAK explanations
                #Discrepancies
                 Lisaname == "X12_copepodBiomass_WGoA" ~ "many small copepods, sum of small and large copepods S of Seward AK",
                 Lisaname == "X12_copepodCom_WGoA" ~ "loads neg on lower trophic DFA, so pos on mid-trophic DFA, ratio of large copepods to (small+large)",
                 Lisaname == "X13_hexagram_EAI" ~ "Tufted puffin chick diet, similar to seabird DFA",
                 Lisaname == "X12_DFA_biomassEuphShelfSum" ~ "Ecosystem proxy (krill ~ pos PDO, Evans et al 2023)",
                 Lisaname == "X13_ammod_WAI" ~ "Puffin chick diet, negligible since 2010",
                 Lisaname == "X13_DFA_WGOA_DFA_seabirds" ~ "Birds have improved since ~2008"
                
                #Confirmations
                Lisaname == "X12_copepodCom_EGoA" ~ "Direct trophic support (cold-water lipid-rich copepod signal)",
                Lisaname == "X14_pinkSalmon" ~ "Competitor (more pink salmon = less food for Chinook salmon)",
                Lisaname == "X13_DFA_WGOA_DFA_midTrophic" ~ "MidTrophic = shrimp/jellyfish, potential competitors; neg loading of capelin & eulachon; neg of lower trophic larval fish",
                
                # Fallback default
                TRUE ~ "Investigate neighborhood relationships further"
              )
            ) %>%
            
            # Streamline, make Node an ordered factor, and add physical sorting column
            select(
              Model   = model_id, 
              Node    = rhs,
              Node2,
              Indicator = Lisaname, 
              Type, 
              Count   = Times_In_Top100, 
              Pct_Pos = Pct_Positive, 
              Pct_Neg = Pct_Negative, 
              Hyp, 
              OK, 
              Pos_Neigh, 
              Neg_Neigh, 
              Explanation
            ) %>%
            # 1. Set up the numeric sorting column and set factor level order
            mutate(
              order = case_when(
                Node2 == "Abundance" ~ 0,
                Node2 == "Growth"    ~ 0.5,
                Node2 == "Competitor"    ~ 0.75,
                Node2 == "PreyNCC"   ~ 1,
                Node2 == "PredNCC"   ~ 2,
                Node2 == "PreyAK"    ~ 3,
                Node2 == "PredAK"    ~ 4,
                TRUE                ~ 99
              ),
              Node = factor(Node, levels = c("Abundance", "Growth", "PreyNCC", "PredNCC", "PreyAK", "PredAK")),
              Node2 = factor(Node2, levels = c("Abundance", "Growth","Competitor", "PreyNCC", "PredNCC", "PreyAK", "PredAK"))
            ) %>%
            # 2. Arrange exactly as requested: Trophic Order -> Success (Confirmation) -> Alphabetical Indicator
            arrange(order, desc(OK), Indicator)

          #explore results----   
          print(ecological_synthesis, n = Inf)

 #summarize across models---------         
          ecological_summary <- ecological_synthesis %>% 
            group_by(Indicator) %>%
            summarise(
              # Categorical columns (takes the value from the first row of each group)
              Node2         = first(Node2),
              Hyp         = first(Hyp),
              Type          = first(Type),
              Explanation         = first(Explanation),
              order         = first(order),
              
              # Numeric summaries
              Support       = mean(OK, na.rm = TRUE),
              mean_Positive = mean(Pct_Pos, na.rm = TRUE),
              Count = sum(Count, na.rm = TRUE)
            ) %>%
            relocate(Indicator, Node2,Hyp, Type, Support, mean_Positive,Count,Explanation,order)%>%
            arrange(order,Node2,desc(Support),Type)
          
          print(ecological_summary,n=Inf)
          
          
          # --- 3. Save & View the Summary Table ---
          print(ecological_summary, n = Inf)
          write_csv(ecological_synthesis, "LisaXP/outputs_4/Indicator_Support_Table.csv")
          write_csv(ecological_summary, "LisaXP/outputs_4/Indicator_Summary_Table.csv")
          

          ecological_summary %>% filter(Support==0,Type=="DFA")
          