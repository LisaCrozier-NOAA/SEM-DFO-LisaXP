# Script: 01_data_prep/04_harvest_bird_pinniped_data.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lavaan)
})

# ------------------------------------------------------------------------------
# STEP 0: PATHS
# ------------------------------------------------------------------------------
proj_dir     <- file.path(getwd(), "Rcode_for_paper")
meta_dir     <- file.path(proj_dir, "metadata")
master_out   <- file.path(proj_dir, "Outputs")
data_out_dir <- file.path(master_out, "data")

dir.create(data_out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. PINNIPED PREDATION (Bonneville Tailrace)
# ------------------------------------------------------------------------------
pinn_path <- file.path(meta_dir, "Braun et al 2025 ACOE Pinn Pred Bonn Tailrace Table5.xlsx")

pinn_raw <- read_xlsx(pinn_path, sheet = 1, skip = 3) %>%
  rename(
    pinn_pred_raw = `Percent of\r\nRun Consumed`,
    nChin_eaten   = `Chinook Consumption Estimate`,
    nChin_bonn    = `Bonneville Dam Spring Chinook\r\nPassage`,
    Year          = Year
  ) %>%
  mutate(
    Year = as.numeric(Year),
    smoltyr = Year - 2,
    pinn_pred = if_else(pinn_pred_raw > 1, pinn_pred_raw / 100, pinn_pred_raw)
  ) %>%
  filter(!is.na(smoltyr)) %>%
  select(smoltyr, Year, nChin_bonn, nChin_eaten, pinn_pred)

# Prepend 1998 & 1999 external baselines (preserves numeric typing)
pinn_ext <- tibble(
  smoltyr     = c(1998, 1999),
  Year        = c(2000, 2001),
  nChin_bonn  = c(NA_real_, NA_real_),
  nChin_eaten = c(NA_real_, NA_real_),
  pinn_pred   = c(0.003, 0.003)
) %>%
  bind_rows(pinn_raw) %>%
  arrange(smoltyr) %>%
  filter(smoltyr >= 1998, smoltyr <= 2021) %>%
  mutate(
    f_pinn = -log(pmax(1e-4, 1 - pinn_pred))
  )

# ------------------------------------------------------------------------------
# 2. HARVEST (Zones 1-5, Lagged 2 Years)
# ------------------------------------------------------------------------------
harvest_path <- file.path(meta_dir, "JointStaffReport 2025 SPRINGCHINOOK.table5.csv")
harvest <- read.csv(harvest_path, stringsAsFactors = FALSE) %>%
  rename(
    harvest_z1_to_z5 = Harvest.impact.downstream.of.Bonneville.Dam..Zones.1.5..Grand.Total,
    nmouth           = Upriver.Run,
    nbonn            = Bonneville.Dam.Count
  ) %>%
  filter(!id %in% c("80-84", "85-89", "90-94")) %>%
  mutate(
    return_yr   = as.numeric(id),
    smoltyr     = return_yr - 2,
    u_harv_z1_5 = harvest_z1_to_z5 / nmouth,
    f_harv_z1_5 = -log(pmax(1e-4, 1 - u_harv_z1_5))
  ) %>%
  filter(smoltyr >= 1998, smoltyr <= 2021) %>%
  select(smoltyr, return_yr, nmouth, nbonn, harvest_z1_to_z5, u_harv_z1_5, f_harv_z1_5)

# ------------------------------------------------------------------------------
# 3. ESTUARY AVIAN PREDATION (Smolt Outmigration)
# ------------------------------------------------------------------------------
bird_path <- file.path(meta_dir, "aggregate_annual_estuary_avian_predation_on_chinook_fromBON_1998_2018.csv")
bird <- read.csv(bird_path, stringsAsFactors = FALSE) %>%
  group_by(year) %>%
  summarize(birdmort = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  rename(smoltyr = year) %>%
  mutate(
    f_bird = -log(pmax(1e-4, 1 - birdmort))
  )

# Append 2019, 2020, 2021 baseline projections
bird_ext <- bind_rows(
  bird,
  tibble(
    smoltyr  = c(2019, 2020, 2021),
    birdmort = c(0.0192, 0.0192, 0.0192),
    f_bird   = -log(pmax(1e-4, 1 - c(0.0192, 0.0192, 0.0192)))
  )
) %>%
  arrange(smoltyr) %>%
  filter(smoltyr >= 1998, smoltyr <= 2021)

# ------------------------------------------------------------------------------
# 4. MASTER MERGE & SINGLE IN-RIVER MORTALITY CALCULATION
# ------------------------------------------------------------------------------
data_path <- file.path(meta_dir, "sem_altprey_data_1998_2021_10212026.csv")
model_dat <- read.csv(data_path, stringsAsFactors = FALSE, row.names = 1)

combined_dat <- model_dat %>%
  rename(smoltyr = year) %>%
  left_join(bird_ext %>% select(smoltyr, birdmort, f_bird), by = "smoltyr") %>%
  left_join(harvest %>% select(smoltyr, return_yr, u_harv_z1_5, f_harv_z1_5), by = "smoltyr") %>%
  left_join(pinn_ext %>% select(smoltyr, pinn_pred, f_pinn), by = "smoltyr") %>%
  mutate(
    # --- Single Combined In-River Instantaneous Mortality ---
    # F_total = F_bird (smolt) + F_pinn (adult tailrace) + F_harv (adult Z1-5)
    f_inriver_total = f_bird + f_pinn + f_harv_z1_5,
    
    # Standardize to Z-score for the arena model
    z_f_inriver_total = as.numeric(scale(f_inriver_total)),
    
    # Retain individual standardized terms if needed for diagnostics
    z_f_bird = as.numeric(scale(f_bird)),
    z_f_pinn = as.numeric(scale(f_pinn)),
    z_f_harv = as.numeric(scale(f_harv_z1_5))
  )

# Diagnostic inspection
print(head(combined_dat %>% 
             select(smoltyr, X16_DFA_sar_UCsprCh_smoltyr, f_inriver_total, z_f_inriver_total)))

summary(combined_dat)
# Save output
out_file <- file.path(meta_dir, "sem_data_with_harvest_bird_pinn_predation.csv")
write.csv(combined_dat, out_file, row.names = FALSE)
cat(sprintf("Merged dataset successfully written to: %s\n", out_file))
