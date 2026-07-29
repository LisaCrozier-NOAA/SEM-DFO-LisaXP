# Check completeness of data input into SEM
#
# Doug Jackson
# doug@QEDAconsulting.com

###########################################################################
# Constants
###########################################################################
workingDir <- "/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices"

###########################################################################
# Run
###########################################################################
setwd(workingDir)

outputDir <- file.path(workingDir, "output")

clusDataDFA <- readRDS(file.path(outputDir, "DFA", "clusDataDFA.rds"))

clusDataDFA <- clusDataDFA %>% select(shortName, date, finalVal) %>% 
    pivot_wider(id_cols=date, values_from=finalVal, names_from=shortName)
clusDataDFA$complete <- complete.cases(clusDataDFA)

View(clusDataDFA)

write.csv(clusDataDFA, file.path(outputDir, "completeness.csv"), row.names=F)

