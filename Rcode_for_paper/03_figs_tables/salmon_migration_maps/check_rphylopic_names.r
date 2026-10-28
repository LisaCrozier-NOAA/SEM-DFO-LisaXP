# ==============================================================================
# Script: test_phylopic_uuids_ncc.R
# Purpose: Check and verify PhyloPic UUID availability for Northern California
#          Current / Oregon / Washington coast ecosystem icons.
# ==============================================================================

library(rphylopic)
library(purrr)
library(tibble)

# Candidate queries: testing common plural names first, with scientific fallbacks
candidate_queries <- list(
  # --- Zooplankton (Newport Line) ---
  copepod       = c("copepods", "Copepoda", "Calanoida"),
  crab_megalopa = c("crabs", "Brachyura", "Decapoda", "crab"),
  larval_fish   = c("larval fishes", "Actinopterygii", "Teleostei", "ichthyoplankton"),
  krill         = c("krill", "Euphausiacea", "Euphausia pacifica"),
  
  # --- Gelatinous & Invertebrate ---
  sea_nettle    = c("jellyfishes", "Scyphozoa", "Chrysaora", "sea nettle"),
  market_squid  = c("squids", "Doryteuthis opalescens", "Loliginidae", "squid"),
  
  # --- Seabirds ---
  common_murre  = c("murres", "Uria aalge", "Alcidae", "auks"),
  shearwater    = c("shearwaters", "Ardenna grisea", "Puffinus", "Procellariidae"),
  pelican       = c("pelicans", "Pelecanus", "Pelecanidae"),
  
  # --- Marine Fishes ---
  hake          = c("hakes", "Merluccius productus", "Merlucciidae"),
  rockfish      = c("rockfishes", "Sebastes", "Scorpaenidae"),
  forage_fish   = c("anchovies", "Engraulis mordax", "smelts", "Osmeridae"),
  
  # --- Columbia River Pinnipeds ---
  harbor_seal   = c("hair seals", "Phoca vitulina", "Phocidae"),
  sea_lion      = c("sea lions", "Zalophus", "Otariidae")
)

check_uuid <- function(group_name, queries) {
  for (q in queries) {
    res <- tryCatch(get_uuid(name = q), error = function(e) NULL)
    if (!is.null(res) && length(res) > 0) {
      return(tibble(
        Taxon_Group = group_name,
        Working_Query = q,
        UUID = res,
        Status = "FOUND"
      ))
    }
  }
  return(tibble(
    Taxon_Group = group_name,
    Working_Query = paste(queries, collapse = " | "),
    UUID = NA_character_,
    Status = "NOT FOUND"
  ))
}

cat("Querying PhyloPic for NCC / Oregon / Columbia River taxa...\n\n")

uuid_check_results <- imap_dfr(candidate_queries, ~ check_uuid(.y, .x))

print(uuid_check_results, n = nrow(uuid_check_results))
