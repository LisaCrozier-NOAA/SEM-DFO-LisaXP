# Read current dataset
sem_data_all <- read.csv("Rcode_for_paper/metadata/sem_data_all_1998_2021.csv", stringsAsFactors = FALSE)

# 1. Identify all X15 columns
all_x15_cols <- grep("^X15_", names(sem_data_all), value = TRUE, ignore.case = TRUE)

# 2. Filter for X15 columns that DO NOT end with / contain "smoltYear"
x15_adult_timing_cols <- all_x15_cols[!grepl("smoltyear", all_x15_cols, ignore.case = TRUE)]

# 3. Identify all X10 columns
all_x10_cols <- grep("^X10_", names(sem_data_all), value = TRUE, ignore.case = TRUE)

# 2. Filter for X15 columns that DO NOT end with / contain "smoltYear"
x10_adult_timing_cols <- all_x10_cols[!grepl("smoltyear", all_x10_cols, ignore.case = TRUE)]

# Print results
cat("Found", length(x15_adult_timing_cols), "X15 adult-timing columns without smoltYear:\n\n")
print(x15_adult_timing_cols)
print(x10_adult_timing_cols)


#Data review
#1. make sure clusDataDFA_wide_1998_2021_LisaNames.csv matches Doug's completeness file
