library(lavaan)
# Yes, there is an incredibly easy way to identify exactly which active variable names are missing from your my_pretty_names translation list.
# 
# We can fetch all the observed variable names across all 12 model objects, pool them together to find the unique ones currently in use, and then check which ones aren't keys in your my_pretty_names vector.

my_pretty_names <- c(
  "X01_habCompInd"  = "Habitat Compress",
  "X01.ZooPreyNCC_JSOES_DFA_sumPreyOfPrey_planktonJun"  = "JSOES Plankton",
  "X01.ZooPreyNCC_JSOES_Cassin_s_aukl_WS"  = "Cassin's auklet",
  "X03.FishPreyNCC_b_HakeAge1"           = "FishPrey_HakeAge1",
  "X04_marketsquid_GAM"  = "Market squid_GAM",
  "X05_DFA_abundSardine"  = "Sardine DFA",
  "X08_commonMurre_JSOES" = "Murre JSOES",
  "X08.PredBirdNCC_hump_Common_murre_WS" = "Common Murre WS",
  "X09_DFA_HakeAge5Plus"  = "Hake/Mackerel",
  "X10_Harbour_s_2yrLead_WS"  = "Harbour_seal_2yrLead_WS",
  "X10_Harbor_seal_CR_2yrLead"  = "Harbour_seal_2yrLead_CR",
  "X10_Northern_f_s_2yrLead_WS"  = "No_fur_seal_WS",
  "X12_DFA_biomassEuphShelfSum"  = "WCVI Euph/Amphipods",
  "X12_copepodBiomass_WGoA"  = "WGOA copepod Biomass",
  "X12_copepodCom_EGoA"  = "EGOA copepod Biomass",
  "X14_pinkSalmon"  = "Pink Salmon",
  
#These were added this round:  
  "X05_herring_GAM"  = "X05_herring_GAM",
  "X13_pollockBiomassGoAage3plus_predAK"  = "Pollock",
  "X15_DFA_sleeperSharkBSAI_predAK"  = "Sleeper Shark BSAI",
  "X15_salmonSharkGoA_predAK"  = "Salmon Shark GOA",
  "X15_PacificCodBiomass_predAK"  = "Pacific Cod",
    

  "X06_Lmu_IntSprJunH"         = "Lmu_JunH",
  "X06_StomFull_May"           = "StomFull_May",
  "X06_DFA_IGF_mu"             = "Salmon Growth",
  "X07_DFA_cpue_IntSprJunHW"   = "Salmon CPUE",
  "X16_SAR"                    = "Salmon SAR"
)


# 1. Collect all model objects that actually exist right now
all_model_names <- c(
  paste0(model_bases, "_topmodel"),
  paste0(model_bases, "_topmodel_reduced")
)

# 2. Extract every single variable name currently being used across those models
active_variables <- unique(unlist(lapply(all_model_names, function(m_name) {
  if (exists(m_name)) {
    lavNames(get(m_name), type = "ov")
  }
})))

# 3. Find which of those active variables are MISSING from your pretty names list
missing_names <- setdiff(active_variables, names(my_pretty_names))

# 4. Print out the culprits
if (length(missing_names) == 0) {
  message("Success! Your my_pretty_names list is complete. No names are missing.")
} else {
  message("The following variables are in your models but missing from your pretty names list:")
  print(missing_names)
}