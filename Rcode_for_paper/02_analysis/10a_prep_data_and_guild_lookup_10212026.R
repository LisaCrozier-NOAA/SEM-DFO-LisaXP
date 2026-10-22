

# ==============================================================================
# Script: 10a_prep_data_and_guild_lookup_10212026.R
# Directory: Rcode_for_paper/02_analysis/
# Purpose: Data prep, 19yr/24yr filtering, Min-Max [0, 1] scaling, CSV export,
#          and dynamic Guild Lookup Table creation with Hake self-exclusion 
#          and Euphausiids in the alternative prey pool.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------------------------
# STEP 0: CONFIGURATION SWITCH & PATHS
# ------------------------------------------------------------------------------
years_included <- "19yr" # Switch: "19yr" or "24yr"

proj_dir     <- file.path(getwd(), "Rcode_for_paper")
master_out   <- file.path(proj_dir, "Routput_for_paper")
data_out_dir <- file.path(master_out, "data")

orig_raw_df <- read.csv(file.path(proj_dir, "metadata/sem_altprey_data_1998_2021.csv"), stringsAsFactors = FALSE,row.names = 1);names(orig_raw_df)
names(orig_raw_df)[which(grepl("hake|sardine",names(orig_raw_df),ignore.case = TRUE))]


raw_csv_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
raw_df       <- read.csv(raw_csv_path, stringsAsFactors = FALSE,row.names = 1)
names(raw_df)

colnames(raw_df) <- ifelse(
  colnames(raw_df) != "year" & !grepl("_smoltyr$|_adultyr$", colnames(raw_df)),
  paste0(colnames(raw_df), "_smoltyr"),
  colnames(raw_df)
)

names(raw_df)


# Clean empty or unnamed columns
raw_df <- raw_df[, names(raw_df) != "" & !is.na(names(raw_df))]


# ------------------------------------------------------------------------------
# STEP 1: APPLY 19-YEAR VS 24-YEAR COLUMN FILTERING
# ------------------------------------------------------------------------------
if (years_included == "24yr") {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) == 24]
} else {
  valid_cols <- names(raw_df)[colSums(!is.na(raw_df)) >= 19]
}

filtered_df <- raw_df %>% 
  select(all_of(valid_cols)) %>% 
  filter(complete.cases(.))

cat(sprintf("[%s Selection] Retained %d columns across %d complete rows.\n", 
            years_included, ncol(filtered_df), nrow(filtered_df)))

# ------------------------------------------------------------------------------
# STEP 2: MIN-MAX SCALING [0, 1] & CSV EXPORT
# ------------------------------------------------------------------------------
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

minmax_complete_data <- filtered_df %>%
  mutate(across(where(is.numeric) & !matches("Year|year"), min_max_scale))

write_csv(minmax_complete_data, file.path(data_out_dir, paste0("minmax_complete_data_10212026_", years_included, ".csv")))

cor(minmax_complete_data$X09_DFA_HakeAge5Plus_smoltyr,minmax_complete_data$X21_DFA_SST_NCC_aprjun_DFObccm_smoltyr) #0.9307828
cor(minmax_complete_data$X09_DFA_HakeAge5Plus_smoltyr,minmax_complete_data$X21_sst_egoa_junjulaug_smoltyr) #0.8835028


all_cols <- names(filtered_df)
predator_cols <- sort(all_cols[grepl("^(X08_|X09_|X10_|X11_|X15_)", all_cols)])
 predator_cols
 
 all_cols <- names(minmax_complete_data)
 predator_cols <- sort(all_cols[grepl("^(X08_|X09_|X10_|X11_|X15_)", all_cols)])
 predator_cols
 
 # ==============================================================================
 # STEP 3: DEFINE FUNCTIONAL PREY POOLS & SPATIAL OVERLAP BY GUILD
 # ==============================================================================
 all_cols <- names(minmax_complete_data)
 
 get_prey_pool <- function(guild_num, stage_suffix = "smoltyr", pred_col = "") {
   # Clean stage suffix formatting
   sfx <- if (nzchar(stage_suffix)) paste0("_", stage_suffix) else "_smoltyr"
   
   m_squid     <- paste0("X04_marketsquid_GAM", sfx)
   anchovy     <- paste0("X05_anchovy_GAM", sfx)
   sardine     <- paste0("X05_DFA_abundSardine", sfx)
   herring_ncc <- paste0("X05_herring_NCC", sfx)
   eulachon    <- paste0("X05_eulachon_during_chinook", sfx)
   hake        <- paste0("X09_DFA_HakeAge5Plus", sfx)
   
   euphausiids <- paste0("biomassEuphShelfSum", sfx)
   
   pollock     <- paste0("X13_pollockGOA_age1and2", sfx)
   herring_ak  <- paste0("sitkaHerring_EGoA", sfx)
   capelin     <- paste0("capelin_WGoA", sfx)
   krill_ak    <- paste0("X12_egoa_krill", sfx)
   pink_ak     <- paste0("X14_pinkSalmonNorthAmerica", sfx)
   
   prey_list <- switch(
     as.character(guild_num),
     "1" = c(m_squid, anchovy, sardine, herring_ncc, eulachon),
     "2" = c(euphausiids, m_squid, anchovy, sardine, herring_ncc, eulachon, hake),
     "3" = c(pollock, herring_ak, capelin, krill_ak),
     "4" = c(pollock, herring_ak, capelin),
     "5" = c(pink_ak, pollock, herring_ak, capelin),
     c()
   )
   
   # Exclude self/predator column if present in pool
   prey_list <- setdiff(prey_list, pred_col)
   
   # Keep only columns available in current minmax dataset
   intersect(prey_list, all_cols)
 }
 
 get_spatial_overlap_col <- function(pred_col, region, stage_suffix = "smoltyr") {
   sfx <- if (nzchar(stage_suffix)) paste0("_", stage_suffix) else "_smoltyr"
   
   if (region == "NCC") {
     col <- paste0("X21_sst_egoa_junjulaug", sfx)
   } else {
     col <- paste0("X21_sst_egoa_junjulaug", sfx)
   }
   
   if (col %in% all_cols) return(col) else return(NA_character_)
 }
 
 
 # ==============================================================================
 # STEP 4: DYNAMICALLY BUILD GUILD LOOKUP TABLE
 # ==============================================================================
 predator_cols <- all_cols[grepl("^(X08_|X09_|X10_|X11_|X15_)", all_cols)]
 
 blacklist <- c("X09_DFA_ChinAbundSnakeFall")
 predator_cols <- predator_cols[!sapply(predator_cols, function(x) any(sapply(blacklist, grepl, x)))]
 
 lookup_rows <- list()
 
 for (pred in predator_cols) {
   # 1. Determine stage_suffix from predator column name
   stage_suffix <- if (grepl("_adultyr$", pred)) "adultyr" else "smoltyr"
   
   # 2. Determine region
   region <- if (grepl("^(X08_|X09_|X10_)|Harbour_p", pred, ignore.case = TRUE)) {
     "NCC"
   } else {
     "AK"
   }
   
   # 3. Assign guild numbers based on region and predator type
   guild_num <- case_when(
     region == "NCC" & grepl("^(X08_|Harbour_p)", pred, ignore.case = TRUE) ~ 1,
     region == "NCC" ~ 2,
     region == "AK" & grepl("(ssl|sablefish|spinyDogfish)", pred, ignore.case = TRUE) ~ 3,
     region == "AK" & grepl("(halibut|PacificCod|Arrowtooth)", pred, ignore.case = TRUE) ~ 4,
     region == "AK" & grepl("shark", pred, ignore.case = TRUE) ~ 5,
     .default = 3
   )
   
   # 4. Retrieve candidate prey pool and spatial overlap column with explicit argument names
   prey_candidates <- get_prey_pool(
     guild_num    = guild_num, 
     stage_suffix = stage_suffix, 
     pred_col     = pred
   )
   
   spatial_col <- get_spatial_overlap_col(
     pred_col     = pred, 
     region       = region, 
     stage_suffix = stage_suffix
   )
   
   # Combine non-empty candidates into string
   all_cands <- c(spatial_col, prey_candidates)
   all_cands <- all_cands[!is.na(all_cands) & all_cands != ""]
   
   lookup_rows[[length(lookup_rows) + 1]] <- tibble(
     region              = region,
     guild_num           = guild_num,
     pred_data_col       = pred,
     stage_suffix        = stage_suffix,
     spatial_overlap_col = spatial_col,
     all_candidates_str  = paste(all_cands, collapse = ", ")
   )
 }
 
 guild_pred_lookup <- bind_rows(lookup_rows)
 
 print(guild_pred_lookup, n = Inf)
 
  write_csv(guild_pred_lookup, file.path(data_out_dir, paste0("guild_predator_lookup_table_", years_included, ".csv")))
  