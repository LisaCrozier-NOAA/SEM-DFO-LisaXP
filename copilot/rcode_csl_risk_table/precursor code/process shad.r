# ==============================================================================
# Dedicated Script: Clean & Interpolate Daily Bonneville Shad Passage
# ==============================================================================

library(tidyverse)
library(lubridate)
library(zoo) # For linear interpolation (na.approx)

# -----------------------------------------------------------------------------
# 1. Read Raw Wide File & Drop Feb 29
# -----------------------------------------------------------------------------

shad_raw <- read.csv("data_Lisa/SHAD_bonn_19382026.csv", row.names = NULL, stringsAsFactors = FALSE)

# Filter out Feb 29 upfront
shad_clean_wide <- shad_raw %>%
  filter(!grepl("29-Feb|29.Feb|Feb.29", mm.dd, ignore.case = TRUE))

# -----------------------------------------------------------------------------
# 2. Pivot to Long Format & Standardize Column Names
# -----------------------------------------------------------------------------

shad_long <- shad_clean_wide %>%
  pivot_longer(
    cols      = -mm.dd,
    names_to  = "year_col",
    values_to = "raw_count"
  ) %>%
  mutate(
    # Parse 4-digit year from column names (e.g. "X2024.BON..." -> 2024)
    year = parse_number(year_col),
    
    # Strip spaces and commas from count values
    count_clean = gsub("[^0-9.]", "", trimws(raw_count)),
    
    # Convert empty strings or non-numerics to NA for interpolation
    count_num   = if_else(count_clean == "", NA_real_, as.numeric(count_clean))
  ) %>%
  filter(year >= 1998, year <= 2024)

# -----------------------------------------------------------------------------
# 3. Construct Safe Dates & Linearly Interpolate Missing Counts
# -----------------------------------------------------------------------------

shad_interpolated <- shad_long %>%
  # Force standard year-day string format
  mutate(
    date_str   = paste(mm.dd, year, sep = "-"),
    date_clean = dmy(date_str)
  ) %>%
  # Drop any residual unparseable dates
  filter(!is.na(date_clean)) %>%
  arrange(year, date_clean) %>%
  group_by(year) %>%
  mutate(
    # Linear interpolation for interior NAs; rule=2 carries boundary values
    count_interp = na.approx(count_num, x = date_clean, na.rm = FALSE, rule = 2),
    # If a full year is NA, default to 0
    count_interp = replace_na(count_interp, 0),
    
    week_num = week(date_clean)
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 4. Aggregate Weekly Sums (Weeks 10–26)
# -----------------------------------------------------------------------------

shad_weekly_final <- shad_interpolated %>%
  filter(week_num >= 10, week_num <= 26) %>%
  group_by(year, week = week_num) %>%
  summarise(
    Shad_weekly_count = sum(count_interp, na.rm = TRUE),
    Shad_weekly_mean  = mean(count_interp, na.rm = TRUE),
    .groups           = "drop"
  )

# Ensure full grid for weeks 10–26 across all years
shad_weekly_grid <- expand_grid(
  year = 1998:2024,
  week = 10:26
) %>%
  left_join(shad_weekly_final, by = c("year", "week")) %>%
  mutate(
    Shad_weekly_count = replace_na(Shad_weekly_count, 0),
    Shad_weekly_mean  = replace_na(Shad_weekly_mean, 0)
  )

cat("--- Processed Shad Dataset Summary (1998–2024) ---\n")
print(head(shad_weekly_grid))
print(summary(shad_weekly_grid$Shad_weekly_count))


write.csv(shad_weekly_grid, file.path(output_dir, "bonn_shad_weekly_interpolated_1998_2024.csv"), row.names = FALSE)
