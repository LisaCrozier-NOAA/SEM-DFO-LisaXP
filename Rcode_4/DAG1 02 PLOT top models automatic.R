#This script uses previously built models and 2 plot layout functions and plotting names created elsewhere
#to plot all of the topmodels and the reduced forms of those models.

#load if starting from scratch----------
source("LisaXP/functions/sem_plot_fxn.r")
source("LisaXP/functions/plot_layout_builder_fxn.r")
source("LisaXP/functions/extract_nodes_for_plot_layout_fxn.r")
load(file="LisaXP/outputs_4/DAG1.topmodels_noHCI.rdata", verbose=T)
#source("pretty_names.r")

my_pretty_names <- c(
  "X01_habCompInd"  = "Habitat Compress",
  "X01.ZooPreyNCC_JSOES_DFA_sumPreyOfPrey_planktonJun"  = "JSOES Plankton",
  "X01.ZooPreyNCC_JSOES_Cassin_s_aukl_WS"  = "Cassin's auklet",
  "X03.FishPreyNCC_b_HakeAge1"           = "FishPrey_HakeAge1",
  "X04_marketsquid_GAM"  = "Market squid_GAM",
  "X05_DFA_abundSardine"  = "Sardine DFA",
  "X08_commonMurre_JSOES" = "Murre JSOES",
  "X08.PredBirdNCC_hump_Common_murre_WS" = "Common Murre WS",
  "X09_DFA_HakeAge5Plus"  = "Hake/Mackerel DFA",
  "X10_Harbour_s_2yrLead_WS"  = "Harbour_seal_2yrLead_WS",
  "X10_Harbor_seal_CR_2yrLead"  = "Harbour_seal_2yrLead_CR",
  "X10_Northern_f_s_2yrLead_WS"  = "No_fur_seal_WS",
  "X12_DFA_biomassEuphShelfSum"  = "WCVI Euph/Amphipods",
  "X12_copepodBiomass_WGoA"  = "WGOA copepod Biomass",
  "X12_copepodCom_EGoA"  = "EGOA copepod Biomass",
  "X14_pinkSalmon"  = "Pink Salmon",
  
  "X05_herring_GAM"  = "Herring_GAM",
  "X13_pollockBiomassGoAage3plus_predAK"  = "Pollock GOA",
  "X13_ammod_WAI" = "Sandlance WAI",
  "X15_DFA_sleeperSharkBSAI_predAK"  = "Sleeper Shark BSAI",
  "X15_salmonSharkGoA_predAK"  = "Salmon Shark GOA",
  "X15_PacificCodBiomass_predAK"  = "Pacific Cod GOA",
  
  "X06_Lmu_IntSprJunH"         = "Lmu_JunH",
  "X06_StomFull_May"           = "StomFull_May",
  "X06_DFA_IGF_mu"             = "IGF",
  "X07_DFA_cpue_IntSprJunHW"   = "Salmon CPUE",
  "X16_SAR"                    = "Salmon SAR"
)



#Run models-------------
# Define the 6 base structural configurations
#model_bases <- c("DAG1A_short", "DAG1A_long", "DAG1B_long", "DAG1B_short", "DAG1C_long", "DAG1C_short_noHCI")
model_bases <- c( "DAG1C_short_noHCI")

for (m_base in model_bases) {
  
  # PART A: Standard Topmodels (6 Plots)-----------
  model_name <- paste0(m_base, "_topmodel")
  
  if (exists(model_name)) {
    message(paste("Processing Standard Model:", model_name))
    fit_obj <- get(model_name)
    
    # 1. Dynamically discover surviving nodes (Step 1)
    node_map <- extract_active_layout_nodes(fit_obj, var_lookup_NCC_AK)
    
    # 2. Build the matrix layout layout for this model type (Step 2)
    path_layout <- build_dynamic_layout(model_name, node_map)
    assign(paste0("path_layout_", m_base), path_layout, envir = .GlobalEnv)
    
    # 3. Render the structural graph
    p <- Lisa_sem_graph.v2(
      fit_obj, 
      layout = path_layout, 
      title = model_name,
      node_labels = my_pretty_names, 
      r2_nodes = c("X06_DFA_IGF_mu", "X07_DFA_cpue_IntSprJunHW", "X16_SAR"),
      width = 10, height = 6, angle = 180, rect_width = 3, rect_height = 0.8,
      ellipses_width = 2, ellipses_height = 1.2, variance_diameter = .7, 
      text_size = 4.5, curvature = 10, node_width = 10, node_height = 1.5,
      save = TRUE,
      savename = paste0("LisaXP/outputs_4/", model_name, ".png")
    )
    
    # 4. Polish with ggplot elements and assign to global env
    p3 <- plot(p) + 
      labs(title = model_name) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    assign(paste0(m_base, "_Dougtopmodel_plot"), p3, envir = .GlobalEnv)
  }
  
  # PART B: Reduced Topmodels (6 Plots)-----------
  reduced_name <- paste0(m_base, "_topmodel_reduced")
  
  if (exists(reduced_name)) {
    message(paste("Processing Reduced Model:", reduced_name))
    fit_obj_red <- get(reduced_name)
    
    # 1. Discover surviving nodes for the reduced model variant
    node_map_red <- extract_active_layout_nodes(fit_obj_red, var_lookup_NCC_AK)
    
    # 2. Build the matrix layout tailored to the remaining elements
    path_layout_red <- build_dynamic_layout(reduced_name, node_map_red)
    
    # 3. Render the reduced structural graph
    p.reduced <- Lisa_sem_graph.v2(
      fit_obj_red, 
      layout = path_layout_red, 
      title = reduced_name,
      node_labels = my_pretty_names, 
      r2_nodes = c("X06_DFA_IGF_mu", "X07_DFA_cpue_IntSprJunHW", "X16_SAR"),
      width = 10, height = 6, angle = 180, rect_width = 3, rect_height = 0.8,
      ellipses_width = 2, ellipses_height = 1.2, variance_diameter = .7, 
      text_size = 4.5, curvature = 10, node_width = 10, node_height = 1.5,
      save = TRUE,
      savename = paste0("LisaXP/outputs_4/", reduced_name, ".png")
    )
    
    # 4. Polish with ggplot elements and assign to global env
    p3.reduced <- plot(p.reduced) + 
      labs(title = reduced_name) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    assign(paste0(m_base, "_Dougtopmodel_reduced_plot"), p3.reduced, envir = .GlobalEnv)
  }
}



