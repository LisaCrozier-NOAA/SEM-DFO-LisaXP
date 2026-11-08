

library(tidyverse)

meta_df <- read_csv("Rcode_for_paper/Outputs/data/ak_dfa_constituents_metadata_summary.csv")

# Filter to robust series (<= 2 missing years) and inspect by subguild
meta_df %>%
  filter(missing_years <= 2) %>%
  count(region, subguild) %>%
  arrange(subguild, region)


# region subguild          n
# 1 EGOA   Lower Trophic     8
# 2 WGOA   Lower Trophic     7
# 3 EGOA   Mid Trophic       9
# 4 WGOA   Mid Trophic       5
# 5 EGOA   Seabird           4
# 6 WGOA   Seabird           1

# Check which Alaska prey variables met the importance threshold
val_ranks_wide <- read_csv("Rcode_for_paper/Outputs/data/obj2_data_variable_importance_wide.csv") 
  print(val_ranks_wide %>% filter(grepl("X14|X13|X12",LisaName)) %>% arrange(desc(DAG1B_long)),n=Inf)
  
#this one includes SEMnode (generic_node)  
val_ranks <- read_csv("Rcode_for_paper/Outputs/data/obj2_data_valid_rankings_long.csv") # or your valid_rankings_long object

val_ranks %>%
  filter(generic_node == "PreyAK") %>%
  select(var_name, LisaName, contains("DAG1B"))


library(tidyverse)

dfa_loadings_path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Ferris-DFAIndicators-goa/MARSS results/EGOA/egoa_ltl_loadings.csv"

loadings_df <- read_csv(dfa_loadings_path, show_col_types = FALSE)

# Inspect loadings ranked by absolute value
loadings_ranked <- loadings_df %>%
  mutate(abs_loading = abs(loading)) %>%
  arrange(desc(abs_loading))

print(loadings_ranked)


# Map readable names back to data column codes
name_to_code <- c(
  "Chla biomass"               = "EGOA_CHLA.BIOM",
  "Chla peak bloom"            = "EGOA_CHLA.PEAK",
  "Copepod community size"     = "EGOA_SECM.COP.COMM.SIZE",
  "Euphausiid abundance"       = "EGOA_SECM.EUPH.DENS",
  "Gastropod abundance"        = "EGOA_SECM.GAST.DENS",
  "Amphipod abundance"         = "EGOA_SECM.HYP.AMPH",
  "Large copepod abundance"    = "EGOA_SECM.LG.CAL.COP.DENS",
  "Mesozooplankton abundance"  = "EGOA_SECM.MESOZOOP.DENS",
  "Small copepod abundance"    = "EGOA_SECM.SM.CAL.COP.DENS"
)

# Pull the top 4 based on DFA loadings
top_4_names <- head(loadings_ranked$names, 4)
starting_indicators <- unname(name_to_code[top_4_names])

cat("Starting SEM Latent_PreyAK with top 4 DFA drivers:\n")
print(starting_indicators)
