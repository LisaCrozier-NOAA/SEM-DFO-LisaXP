# Specifications for screening criteria
screenStartDate <- "01JAN1998"
screenEndDate <- "31DEC2021"

# Datasets to use
subDirs <- c("Western_Aleutian_Islands", "Central_Aleutian_Islands", "Eastern_Aleutian_Islands",
             "Western_Gulf_of_Alaska", "Eastern_Gulf_of_Alaska", "CCIEA", "WCVI", "SEM_data_2024",
             "pinkSalmon", "BonnPinn", "ColumbiaRiverPinn", "SEM_data_2025", "pinkSalmon_2025",
             "data_2025", "NHL_2025", "planktonJuneNCC_2025", "plankton_2025b",
             "plankton_2026", "NHL_2026", "JSOES_2026", "Westport_2026", "WGOA_DFA_2026", "lingcodSA_2026",
             "predAK_2026", "orca_2026")

getInd <- function(indicators, thisIndicator, thisDataset) {
    thisInd <- indicators %>% filter(indicator==thisIndicator, dataset==thisDataset)
    return(thisInd)
}

# Fill in missing values, e.g., using interpolation
impute <- function(thisInd, thisData) {
    if(nrow(thisInd)>0 & thisInd$impute!="none") {
        thisData$imputed <- thisData$value
        if(thisInd$impute=="linear") {
            thisData$imputed[thisData$imputed<=0] <- NA
        } else if(thisInd$impute=="linearKeepZeros") {
            thisData$imputed[thisData$imputed<0] <- NA   
        } else if(thisInd$impute=="linearKeepNegative") {
            # Do nothing
        }
        isFiniteIndices <- which(is.finite(thisData$imputed))
        thisData$imputed <- na_interpolation(thisData$imputed)
        thisData$imputed[1:max(1, isFiniteIndices[1]-1)] <- NA
        thisData$imputed[min(max(isFiniteIndices)+1, nrow(thisData)):nrow(thisData)] <- NA
        thisData$imputed[1] <- thisData$value[1]
        thisData$imputed[nrow(thisData)] <- thisData$value[nrow(thisData)]
    } else {
        thisData$imputed <- thisData$value
    }
    return(thisData)
}

# Natural log transform
logTransform <- function(thisInd, thisData) {
    if(nrow(thisInd)>0 & thisInd$logTransform=="simple") {
        thisData$logTransformed <- log(thisData$imputed)
        transformed <- TRUE
    } else {
        thisData$logTransformed <- NA
        transformed <- FALSE
    }
    return(list(thisData=thisData, transformed=transformed))
}

