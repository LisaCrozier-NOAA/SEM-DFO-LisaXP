



# ==============================================================================
# Script: 02_analysis/06b_build_sem_latent_modeling_dataset.R
# Purpose: Extract all variables retained from the cluster analysis screening
#          (max > 0.05, missing <= 2 years across 1998-2021), preserve explicit
#          'Peff_' prefixes for all predators, include core DAG anchors,
#          and export raw and standardized modeling dataframes.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------------------------
# 0. SETUP & INGEST RAW SOURCES
# ------------------------------------------------------------------------------
proj_dir <- "Rcode_for_paper"
data_dir <- file.path(proj_dir, "Outputs/data")

alt_data_path <- file.path(proj_dir, "metadata/sem_altprey_data_1998_2021_10212026.csv")
peff_ts_path  <- file.path(data_dir, "obj3_data_peff_timeseries_24yr.csv")
audit_path    <- file.path(data_dir, "obj4_semnode_data_filtering_audit.csv")

if (!file.exists(alt_data_path)) stop("Missing: ", alt_data_path)
if (!file.exists(peff_ts_path))  stop("Missing: ", peff_ts_path)
if (!file.exists(audit_path))    stop("Missing audit table. Run cluster filtering first: ", audit_path)

df_raw   <- read.csv(alt_data_path, row.names = 1, stringsAsFactors = FALSE)
peff_df  <- read_csv(peff_ts_path, show_col_types = FALSE)
audit_df <- read_csv(audit_path, show_col_types = FALSE)

# Align to 1998-2021 (n = 24 years)
yr_col <- grep("^year$", names(df_raw), value = TRUE, ignore.case = TRUE)
df_raw$year <- if (length(yr_col) > 0) as.numeric(df_raw[[yr_col[1]]]) else as.numeric(rownames(df_raw))
df_24 <- df_raw %>% filter(year >= 1998 & year <= 2021)

yr_peff_col <- grep("^year$", names(peff_df), value = TRUE, ignore.case = TRUE)
if (length(yr_peff_col) > 0) {
  peff_df <- peff_df %>% filter(.data[[yr_peff_col[1]]] >= 1998 & .data[[yr_peff_col[1]]] <= 2021)
}

find_col <- function(pattern, df) {
  matches <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
  if (length(matches) == 0) stop("Could not find column matching: ", pattern)
  matches[1]
}

# ------------------------------------------------------------------------------
# 1. INITIALIZE MASTER DATASET WITH CORE ANCHORS
# ------------------------------------------------------------------------------
df_master <- tibble(
  year        = df_24$year,
  SAR         = df_24[[find_col("sar_UCsprCh", df_24)]],
  CPUE        = df_24[[find_col("cpue_IntSprJunHW", df_24)]],
  Growth      = df_24[[find_col("IGF_mu", df_24)]],
  Pink_Salmon = df_24[[find_col("pinkSalmonNorthAmerica", df_24)]],
  Hake_DFA = df_24[[find_col("X09_DFA_HakeAge5Plus_smoltyr", df_24)]]
)

# ------------------------------------------------------------------------------
# 2. EXTRACT RETAINED VARIABLES KEEPING EXPLICIT PEFF PREFIXES
# ------------------------------------------------------------------------------
retained_vars <- audit_df %>%
  filter(Status == "RETAINED") %>%
  pull(LisaName)

cat(sprintf("Extracting %d cluster-screened variables...\n", length(retained_vars)))

for (v in retained_vars) {
  clean_name <- gsub("^X[0-9]+_", "", v)
  
  # Determine if this variable is a predator node
  node_type <- audit_df %>% filter(LisaName == v) %>% pull(Node)
  is_predator <- grepl("Pred", node_type)
  
  ts_vec <- NULL
  final_col_name <- clean_name
  
  if (is_predator) {
    # Match in peff_df, ensuring 'Peff_' prefix is preserved
    peff_match <- grep(paste0("^Peff_.*", clean_name), names(peff_df), value = TRUE, ignore.case = TRUE)
    if (length(peff_match) == 0) {
      peff_match <- grep(clean_name, names(peff_df), value = TRUE, ignore.case = TRUE)
      peff_match <- grep("^Peff_", peff_match, value = TRUE)
    }
    
    if (length(peff_match) > 0) {
      ts_vec <- peff_df[[peff_match[1]]]
      final_col_name <- peff_match[1] # Keeps Peff_ intact
    } else {
      # Fallback if present in df_raw
      raw_match <- grep(clean_name, names(df_24), value = TRUE, ignore.case = TRUE)
      if (length(raw_match) > 0) {
        ts_vec <- df_24[[raw_match[1]]]
        final_col_name <- paste0("Peff_", clean_name)
      }
    }
  } else {
    # Non-predator (prey/competitor) -> extract from df_raw
    raw_match <- grep(paste0("^", clean_name, "$"), names(df_24), value = TRUE, ignore.case = TRUE)
    if (length(raw_match) == 0) {
      raw_match <- grep(clean_name, names(df_24), value = TRUE, ignore.case = TRUE)
    }
    
    if (length(raw_match) > 0) {
      ts_vec <- df_24[[raw_match[1]]]
      final_col_name <- raw_match[1]
    }
  }
  
  if (!is.null(ts_vec)) {
    df_master[[final_col_name]] <- ts_vec
  } else {
    warning("Could not locate timeseries vector for: ", v)
  }
}

# ------------------------------------------------------------------------------
# 3. GENERATE STANDARDIZED (Z-SCORE) VERSION
# ------------------------------------------------------------------------------
df_master_scaled <- df_master %>%
  mutate(across(-year, ~ as.numeric(scale(.x))))

# ------------------------------------------------------------------------------
# 4. EXPORT MASTER MODELING DATASETS
# ------------------------------------------------------------------------------
raw_out_path    <- file.path(data_dir, "sem_modeling_data_screened_raw_1998_2021.csv")
scaled_out_path <- file.path(data_dir, "sem_modeling_data_screened_zscored_1998_2021.csv")

write_csv(df_master, raw_out_path)
write_csv(df_master_scaled, scaled_out_path)

cat("\n====================================================================\n")
cat("FINAL SEM MODELING DATASET GENERATED\n")
cat("====================================================================\n")
cat(sprintf("Observations (Years): %d (1998-2021)\n", nrow(df_master)))
cat(sprintf("Total Columns:        %d\n", ncol(df_master)))
cat("\nColumn Names in Dataset:\n")
print(names(df_master))

cat(sprintf("\nSaved unstandardized dataset to:\n  %s\n", raw_out_path))
cat(sprintf("Saved standardized dataset to:\n  %s\n", scaled_out_path))

