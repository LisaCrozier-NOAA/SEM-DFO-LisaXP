# Preprocess indicator data from various sources
#
# Doug Jackson
# doug@QEDAconsulting.com
library(tidyverse)
library(R.utils)
library(readxl)
###########################################################################
# Constants
###########################################################################
workingDir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices"

AKdataDir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/data/Alaskan_indices"
CCIEAdataDir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/queryERDDAP/CCIEA_download_09dec24"
rawDataDir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/raw_data"
WCVIdataDir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/data/WCVIraw"

# New data from Lisa with more recent years
LisaDir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/SEM-DFO-NMFS/NCCWorkflow/data_2024"
LisaFiles <- c(file.path(LisaDir, "NMFSCondition1.csv"),
               file.path(LisaDir, "NMFSCondition2.csv"),
               file.path(LisaDir, "NMFSPredators.csv"),
               file.path(LisaDir, "NMFSPrey.csv"))

pinkSalmonFile1 <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/data/Ruggerone_et_al/mcf210023-sup-0001-tables1-s24.xlsx"
pinkSalmon1AsiaCols <- c("Korea", "Japan", "M&I", "WKam", "EKam")
pinkSalmon1NorthAmericaCols <- c("WAK", "SPen", "Kod", "CI", "PWS", "SEAK", "NBC", "SBC", "WA")

pinkSalmonFile2 <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/data/Ruggerone_et_al/NPAFC_tech_report.csv"

pinkSalmon2025file <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/data/Ruggerone_et_al/Ruggerone_pink_salmon_2025.csv"

# Bonneville pinniped data file
BonnPinnFile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/data/Bonneville_pinniped/250408_BON_pinniped abundance_2002-2025.xlsx"

# Columbia River pinniped data files
bonFile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/data/Columbia_River_pinniped/Tidwell_250408_BON_pinniped abundance_2002-2025.xlsx"
lcrFile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/data/Columbia_River_pinniped/ODFW atlas count Columbia River v 20250218.xlsx"
wdfwFile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/data/Columbia_River_pinniped/WDFW Harbor Seal Aerial Survey Data - Columbia River Only.csv"

SARfile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/SEM-DFO-NMFS/EvalIndices_DFA_Oct2024data/data/data_ready/SAR.csv"
newSARfile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/SEM-DFO-NMFS/NCCWorkflow/data_2024/sar.19982021.csv"

# New pinniped and salmon data from Lisa
Lisa2025dir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025"

# New versions of SEM data
SEM2025dir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/fromLisa/SEM_data_28may25"

# New data, 2025
# Black Rockfish data file
blackRockfishFile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/fromLisa/BlackRockfish_10jul25/BlackRockfish from Status Report 2023 WA State.xlsx"
FRAMfile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_FRAM/annualCPUE.csv"
JSOESjelliesFile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_JSOES_jellies/NCC2025.CompFishJellies.csv"
bongoFile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_Bongo/Bongo.aggregateannual.4possiblegroupings.csv"
CPSfile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/DougDataProcessScripts_2025/output_CPS/byYear_WA.csv"
hakeASHOPfile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_hakeFromASHOP/hake.springfishery.annualcpue.csv"

# NHL data
krillFile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_krill/krill_NHL.csv"
planktonNHLfile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/planktonNHL.csv"

# PlanktonJuneNCC data
planktonJuneNCCfile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_planktonJuneNCC/planktonJuneNCC.csv"

# plankton_2026 data
plankton2026file <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_plankton/plankton.csv"

# JSOES data
birdFile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_bird/birdNCC.csv"

# Westport survey data
WS2026file <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_Westport/WS_2026.csv"
WS2026groupedFile <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_Westport/WS_2026_grouped.csv"

# WGOA DFA trends
WGOA_DFA_2026file <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_Ferris_AK_DFA/WGO_DFA_2026.csv"

lingcodSA_2026file <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_FRAM/lingcodSA.csv"

predAK_2026file <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/output_processPredAK/predAK_2026.csv"

orca_2026file <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices/LisaDataProcessScripts 2025/raw_data/NRKW_pop_abundance.csv"

# Range of years to include when adding missing years and dates
preprocessMinDate <- "01JAN1998"
preprocessMaxDate <- "31DEC2026"

###########################################################################
# Functions
###########################################################################
addMissingYears <- function(dat) {
    # Datetimes defining preprocessed data range
    preprocessMinDatetime <- dmy(preprocessMinDate)
    preprocessMaxDatetime <- dmy(preprocessMaxDate)
    yearDF <- data.frame(Year=year(seq(preprocessMinDatetime, preprocessMaxDatetime, by="1 year")), ID=1)
    
    dat <- full_join(dat, yearDF, by="Year")
    dat <- dat %>% filter(Year>=year(preprocessMinDatetime), Year<=year(preprocessMaxDatetime)) %>% select(Year, Value) %>% arrange(Year)
    return(dat)
}

addMissingDates <- function(dat) {
    # Datetimes defining preprocessed data range
    preprocessMinDatetime <- dmy(preprocessMinDate)
    preprocessMaxDatetime <- dmy(preprocessMaxDate)
    dateDF <- data.frame(date=seq(preprocessMinDatetime, preprocessMaxDatetime, by="1 year"), ID=1)
    
    dat$date <- as.character(dat$date)
    dateDF$date <- as.character(dateDF$date)
    
    dat <- full_join(dat, dateDF, by="date")
    dat <- dat %>% filter(date>=preprocessMinDatetime, date<=preprocessMaxDatetime) %>% select(date, value) %>% 
        arrange(date)
    return(dat)
}

###########################################################################
# Run
###########################################################################
setwd(workingDir)

source("functions.R")

dataDir <- file.path(workingDir, "data")
success <- unlink(dataDir, recursive=T)
cat("Deleting existing dataDir: success = ", success, "\n")
dir.create(dataDir, showWarnings=F, recursive=T)
CCIEAoutputDir <- file.path(dataDir, "CCIEA")
dir.create(CCIEAoutputDir, showWarnings=F, recursive=T)
WCVIoutputDir <- file.path(dataDir, "WCVI")
dir.create(WCVIoutputDir, showWarnings=F, recursive=T)
pinkSalmonOutputDir <- file.path(dataDir, "pinkSalmon")
dir.create(pinkSalmonOutputDir, showWarnings=F, recursive=T)
LisaSEMoutputDir <- file.path(dataDir, "SEM_data_2024")
dir.create(LisaSEMoutputDir, showWarnings=F, recursive=T)
BonnPinnOutputDir <- file.path(dataDir, "BonnPinn")
dir.create(BonnPinnOutputDir, showWarnings=F, recursive=T)
ColumbiaRiverPinnOutputDir <- file.path(dataDir, "ColumbiaRiverPinn")
dir.create(ColumbiaRiverPinnOutputDir, showWarnings=F, recursive=T)
SEM2025outputDir <- file.path(dataDir, "SEM_data_2025")
dir.create(SEM2025outputDir, showWarnings=F, recursive=T)
pinkSalmon2025outputDir <- file.path(dataDir, "pinkSalmon_2025")
dir.create(pinkSalmon2025outputDir, showWarnings=F, recursive=T)
data2025outputDir <- file.path(dataDir, "data_2025")
dir.create(data2025outputDir, showWarnings=F, recursive=T)
NHL2025outputDir <- file.path(dataDir, "NHL_2025")
dir.create(NHL2025outputDir, showWarnings=F, recursive=T)
planktonJuneNCCoutputDir <- file.path(dataDir, "planktonJuneNCC_2025")
dir.create(planktonJuneNCCoutputDir, showWarnings=F, recursive=T)
plankton2026outputDir <- file.path(dataDir, "plankton_2026")
dir.create(plankton2026outputDir, showWarnings=F, recursive=T)
NHL2026outputDir <- file.path(dataDir, "NHL_2026")
dir.create(NHL2026outputDir, showWarnings=F, recursive=T)
JSOES2026outputDir <- file.path(dataDir, "JSOES_2026")
dir.create(JSOES2026outputDir, showWarnings=F, recursive=T)
WS2026outputDir <- file.path(dataDir, "Westport_2026")
dir.create(WS2026outputDir, showWarnings=F, recursive=T)
WGOA_DFA_2026outputDir <- file.path(dataDir, "WGOA_DFA_2026")
dir.create(WGOA_DFA_2026outputDir, showWarnings=F, recursive=T)
lingcodSA_2026outputDir <- file.path(dataDir, "lingcodSA_2026")
dir.create(lingcodSA_2026outputDir, showWarnings=F, recursive=T)
predAK_2026outputDir <- file.path(dataDir, "predAK_2026")
dir.create(predAK_2026outputDir, showWarnings=F, recursive=T)
orca_2026outputDir <- file.path(dataDir, "orca_2026")
dir.create(orca_2026outputDir, showWarnings=F, recursive=T)

indicators <- read.csv(file.path(workingDir, "indicators.csv"))
indicators <- indicators[indicators$category!="", ]

# Copy AK data into dataDir
subDirs <- list.dirs(AKdataDir, full.names=F)
for(subDir in subDirs) {
    if(subDir!="") {
        copyDirectory(file.path(AKdataDir, subDir), file.path(dataDir, subDir))
    }
}

###########################################################################
# Extract CCIEA data
CCIEAdataLoc <- read.csv(file.path(workingDir, "CCIEAdataLoc.csv"))

shortNameToUrlList <- list()
for(i in 1:nrow(CCIEAdataLoc)) {
    thisInd <- CCIEAdataLoc[i, "indicator"]
    thisDir <- CCIEAdataLoc[i, "dir"]
    thisFile <- CCIEAdataLoc[i, "file"]
    thisLabel <- CCIEAdataLoc[i, "label"]
    thisTimeCol <- CCIEAdataLoc[i, "timeCol"]
    thisLabelCol <- CCIEAdataLoc[i, "labelCol"]
    thisValueCol <- CCIEAdataLoc[i, "valueCol"]
    
    thisShortName <- indicators %>% filter(indicator==thisInd, dataset=="CCIEA") %>% select(shortName) %>% pull(shortName)
    thisSource <- indicators %>% filter(indicator==thisInd, dataset=="CCIEA") %>% select(source) %>% pull(source)
    
    # Skip this data set if its not in indicators.csv
    if(length(thisShortName)==0) {
        next
    }
    
    if(thisDir=="CCIEAdataDir") {
        thisData <- read.csv(file.path(CCIEAdataDir, paste0(thisFile, ".csv")))
    } else if(thisDir=="rawDataDir") {
        thisData <- read_excel(file.path(rawDataDir, paste0(thisFile, ".xlsx")))
        thisData[, thisTimeCol] <- ymd(thisData[, thisTimeCol][[1]], truncated=2L)
        thisData$url <- "NA"
    }
    
    thisData <- thisData[thisData[, thisLabelCol]==thisLabel, ]
    
    thisUrl <- unique(thisData$url)[1]
    shortNameToUrlList[[length(shortNameToUrlList)+1]] <- data.frame(shortName=thisShortName, url=thisUrl)
    
    thisData <- thisData[, c(thisTimeCol, thisValueCol)]
    names(thisData) <- c("date", "value")
    
    thisData <- addMissingDates(thisData)
    
    thisOutputFile <- file.path(CCIEAoutputDir, paste0(thisShortName, ".csv"))
    cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
    cat("Region:,\n", file=thisOutputFile, append=T)
    cat(paste0("Source:,", thisSource, "\n"), 
        file=thisOutputFile, append=T)
    cat("Contact:,\n", file=thisOutputFile, append=T)
    suppressWarnings(write.table(thisData, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))
}

shortNameToUrl <- bind_rows(shortNameToUrlList)
indicatorsWithUrl <- left_join(indicators, shortNameToUrl, by="shortName")
write.csv(indicatorsWithUrl, file=file.path(workingDir, "CCIEAindicatorsWithUrl.csv"), row.names=F)

###########################################################################
# Extract WCVI (DFO) data
WCVIraw <- read.csv(file.path(WCVIdataDir, "dfo_data_out.csv"))

# Eliminate space in metrics
WCVIraw$metric <- sapply(WCVIraw$metric, function(x) gsub(" ", "_", x))

# Create non-lagged versions of anchovy, herring, and sardine
for(ind in c("anchovy", "pacific_herring", "pacific_sardine")) {
    thisRaw <- WCVIraw %>% filter(metric==paste0("lag_log_abund_", ind)) %>% 
        mutate(metric=paste0("log_abund_", ind), year=year-1)
    WCVIraw <- bind_rows(WCVIraw, thisRaw)
}

WCVIind <- unique(WCVIraw[ ,c("metric", "node")])

# Only include predator and prey nodes for now
WCVIind <- WCVIind %>% filter(node %in% c("predator", "prey"))
for(i in 1:nrow(WCVIind)) {
    thisInd <- WCVIind[i, "metric"]
    thisNode <- WCVIind[i, "node"]
    
    thisData <- WCVIraw %>% filter(metric==thisInd, node==thisNode) %>% select(year, value) %>% rename(Year=year, Value=value)
    
    # Ensure that the data contain all years
    thisData <- addMissingYears(thisData)
    
    thisShortName <- indicators %>% filter(indicator==gsub(" ", "_", thisInd), dataset=="WCVI") %>% select(shortName) %>% pull(shortName)
    
    thisOutputFile <- file.path(WCVIoutputDir, paste0(thisShortName, ".csv"))
    cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
    cat("Region:,\n", file=thisOutputFile, append=T)
    cat("Source:,DFO\n", file=thisOutputFile, append=T)
    cat("Contact:,\n", file=thisOutputFile, append=T)
    suppressWarnings(write.table(thisData, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))
}

###########################################################################
# Process pink salmon data
pinkSalmon1 <- read_excel(pinkSalmonFile1, sheet="ST 9-12 Total ret (nos) 52-15", range="A6:P70")
pinkSalmon1Asia <- pinkSalmon1[, c("Year", pinkSalmon1AsiaCols)] %>% mutate(Asia=rowSums(select(., -Year))) %>% select(Year, Asia)
pinkSalmon1NorthAmerica <- pinkSalmon1[, c("Year", pinkSalmon1NorthAmericaCols)] %>% mutate(NorthAmerica=rowSums(select(., -Year))) %>%
    select(Year, NorthAmerica)
pinkSalmon1total <- pinkSalmon1 %>% select(Year, Total)

pinkSalmon1check <- left_join(pinkSalmon1Asia, pinkSalmon1NorthAmerica, by="Year")
pinkSalmon1check <- left_join(pinkSalmon1check, pinkSalmon1total, by="Year")
pinkSalmon1check <- pinkSalmon1check %>% mutate(checkSum=rowSums(select(., Asia, NorthAmerica))-Total)
cat("pinkSalmon1check max difference: ", max(abs(pinkSalmon1check$checkSum)), "\n")

names(pinkSalmon1Asia) <- c("Year", "Value")
names(pinkSalmon1NorthAmerica) <- c("Year", "Value")
names(pinkSalmon1total) <- c("Year", "Value")

pinkSalmon2 <- read_csv(pinkSalmonFile2, show_col_types=F)
pinkSalmon2Asia <- pinkSalmon2 %>% select(Year, "Asia - pink") %>% rename(Value="Asia - pink")
pinkSalmon2NorthAmerica <- pinkSalmon2 %>% select(Year, "North America - pink") %>% rename(Value="North America - pink")
pinkSalmon2total <- pinkSalmon2 %>% select(Year, "Total - pink") %>% rename(Value="Total - pink")

pinkSalmonAsia <- bind_rows(pinkSalmon1Asia, pinkSalmon2Asia)
pinkSalmonNorthAmerica <- bind_rows(pinkSalmon1NorthAmerica, pinkSalmon2NorthAmerica)
pinkSalmonTotal <- bind_rows(pinkSalmon1total, pinkSalmon2total)

thisOutputFile <- file.path(pinkSalmonOutputDir, paste0("pinkSalmonAsia.csv"))
cat(paste0("Dataset:,pinkSalmonAsia\n"), file=thisOutputFile)
cat("Region:,Asia,\n", file=thisOutputFile, append=T)
cat("Source:,Ruggerone et al.\n", file=thisOutputFile, append=T)
cat("Contact:,\n", file=thisOutputFile, append=T)
suppressWarnings(write.table(pinkSalmonAsia, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))

thisOutputFile <- file.path(pinkSalmonOutputDir, paste0("pinkSalmonNorthAmerica.csv"))
cat(paste0("Dataset:,pinkSalmonNorthAmerica\n"), file=thisOutputFile)
cat("Region:,North America,\n", file=thisOutputFile, append=T)
cat("Source:,Ruggerone et al.\n", file=thisOutputFile, append=T)
cat("Contact:,\n", file=thisOutputFile, append=T)
suppressWarnings(write.table(pinkSalmonNorthAmerica, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))

thisOutputFile <- file.path(pinkSalmonOutputDir, paste0("pinkSalmon.csv"))
cat(paste0("Dataset:,pinkSalmon\n"), file=thisOutputFile)
cat("Region:,Asia and North America,\n", file=thisOutputFile, append=T)
cat("Source:,Ruggerone et al.\n", file=thisOutputFile, append=T)
cat("Contact:,\n", file=thisOutputFile, append=T)
suppressWarnings(write.table(pinkSalmonTotal, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))

# Updated data from Ruggerone; 28may25 e-mail
pinkSalmon2025 <- read.csv(pinkSalmon2025file)
pinkSalmonAsia2025 <- pinkSalmon2025 %>% select(Year, Asia) %>% rename(Value="Asia")
pinkSalmonNorthAmerica2025 <- pinkSalmon2025 %>% select(Year, North.America) %>% rename(Value="North.America")
pinkSalmonTotal2025 <- pinkSalmon2025 %>% mutate(Value=Asia + North.America) %>% select(Year, Value)

inds <- c("pinkSalmonAsia_2025", "pinkSalmonNorthAmerica_2025", "pinkSalmon_2025")
dfs <- list(pinkSalmonAsia2025, pinkSalmonNorthAmerica2025, pinkSalmonTotal2025)
for(i in 1:length(inds)) {
    thisInd <- inds[i]
    thisIndicators <- indicators %>% filter(indicator==thisInd)
    thisOutputFile <- file.path(pinkSalmon2025outputDir, paste0(thisIndicators$shortName, ".csv"))
    cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
    cat("Region:,Asia,\n", file=thisOutputFile, append=T)
    cat("Source:,", thisIndicators$source, "\n", file=thisOutputFile, append=T)
    cat("Contact:,\n", file=thisOutputFile, append=T)
    thisDat <- addMissingYears(dfs[[i]])
    suppressWarnings(write.table(dfs[[i]], file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))    
}
###########################################################################
# Obtain smolt-to-adult return ratio data
SAR <- read.csv(SARfile)
names(SAR) <- c("Year", "SAR")
write.csv(SAR, file.path(dataDir, basename(SARfile)), row.names=F)
scratch <- file.copy(newSARfile, file.path(dataDir, basename(newSARfile)))

###########################################################################
# Extract Lisa's data
for(f in LisaFiles) {
    dat <- read.csv(f)
    if("year" %in% names(dat)) {dat <- dat %>% rename(Year=year)}
    ind <- names(dat)
    ind <- ind[ind!="Year"]
    
    for(i in 1:length(ind)) {
        thisInd <- ind[i]
        thisDat <- dat[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        
        thisShortName <- indicators %>% filter(indicator==gsub(" ", "_", thisInd), dataset=="SEM_data_2024") %>% select(shortName) %>% pull(shortName)
        
        thisOutputFile <- file.path(LisaSEMoutputDir, paste0(thisShortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat("Source:,Lisa SEM\n", file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))
        
        # Create 2-year lead version of ssl.est.wholerange
        if(thisInd=="ssl.est.wholerange") {
            thisIndLead <- paste0(thisInd, "_2yrLead")
            thisDat <- dat[, c("Year", thisInd)]
            names(thisDat) <- c("Year", "Value")
            
            thisDat$Year <- thisDat$Year - 2
            
            thisShortName <- indicators %>% filter(indicator==gsub(" ", "_", thisIndLead), dataset=="SEM_data_2024") %>% select(shortName) %>% pull(shortName)
            
            thisOutputFile <- file.path(LisaSEMoutputDir, paste0(thisShortName, ".csv"))
            cat(paste0("Dataset:,", thisIndLead, "\n"), file=thisOutputFile)
            cat("Region:,\n", file=thisOutputFile, append=T)
            cat("Source:,Lisa SEM\n", file=thisOutputFile, append=T)
            cat("Contact:,\n", file=thisOutputFile, append=T)
            thisDat <- addMissingYears(thisDat)
            suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))
        }
    }
}

###########################################################################
# Extract and preprocess Bonneville pinniped data
BonnPinn <- read_excel(BonnPinnFile)

# Use the maximum daily count per spring (prior to July)
# see 30apr25 e-mail from Lisa
BonnPinn <- BonnPinn %>% filter(MONTH<=7) %>% group_by(YEAR) %>% summarize(EJU=max(`EJU Abundance`, na.rm=T),
                                                                        ZCA=max(`ZCA Abundance`, na.rm=T),
                                                                        PVI=max(`PVI Abundance`, na.rm=T)) %>% 
    mutate(Value=EJU + ZCA + PVI) %>% rename(Year=YEAR) %>% select(Year, Value)

BonnPinn <- addMissingYears(BonnPinn)

thisInd <- "Bonneville pinnipeds"
thisIndicators <- indicators %>% filter(indicator==thisInd)
thisOutputFile <- file.path(BonnPinnOutputDir, "BonnPinn.csv")
cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
cat("Region:,\n", file=thisOutputFile, append=T)
cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
cat("Contact:,\n", file=thisOutputFile, append=T)
suppressWarnings(write.table(BonnPinn, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))

###########################################################################
# Extract and preprocess Columbia River pinniped data
# Adaptation of Lisa's script, "Process mammal data for SEM.R"

#Prep Mammal counts for SEM modeling
#Notes: 3 spp (PV: harbor seals, EJ: California sea lions, and ZC: Steller sea lions) are counted at multiple locations by multiple agencies at and below Bonneville Dam
#In this script, I extract the information from just one location for harbor seals (Desdemona Sands)
#and two locations for both sea lion species -- Bonn Dam and East Mooring Basin
#I calculate the maximum daily count conducted by any agency within the spring Chinook smolt or adult migration season
#which I define as April 1 through June 30 for juveniles and hence harbor seals
#Ref for juv: https://www.cbr.washington.edu/dart/wrapper?type=php&fname=hrt_pit_1746042404_902.php
#Apr 1 to Jun 30 for adults and hence sea lions at Bonn
#Ref for adults: https://www.cbr.washington.edu/dart/wrapper?type=php&fname=hrt_pitadult_1746117777_949.php
#and Mar 1 to Jun 30 for adults and hence sea lions at EMB
#the extra month before Bonn passage time for adults acounts for aggregations at the mouth of the estuary for an unknown amount of time before they start directed upstream migration, and they would be vulnerable to EMB predation for that entire time.
#2025 field work was incomplete at the time of data transfer (April 2025)

#1. Sea lions at Bonneville dam, ZC and EJ, Apr 1-Jun 30 --------

#Data from Kyle Tidwell, Bonn counts
# EJU = Steller
# ZCA = California
# PVI = Harbor Seal
bon <- read_xlsx(bonFile, sheet=1)

bon <- bon %>% filter(MONTH>3, MONTH<7, YEAR<2025) %>% group_by(YEAR) %>%
    summarize(CSLbonn=max(`ZCA Abundance`, na.rm=T), SSLbonn=max(`EJU Abundance`, na.rm=T), AllPinnDayBonn=max(`Total Pinniped Abundance`, na.rm=T)) %>% 
    rename(Year=YEAR)

#2. Sea lions at East Mooring Basin, ZC and EJ, Mar 1-Jun 30    ----------
lcr <- read_xlsx(lcrFile, sheet=1)
lcr$Year <- as.numeric(substr(lcr$datemil, 1, 4))
lcr$month <- as.numeric(substr(lcr$datemil, 5, 6))
sl <- lcr %>% filter(month>2, month<7, spp!="PV", location=="COLUMBIA RIVER-EAST MOORING BASIN") 
sl <- tapply(sl$nonpup_total,list(sl$Year,sl$spp), max,na.rm=T)
sl[is.na(sl[,"EJ"]), "EJ"] <- 0
colnames(sl) <- c("SSLemb", "CSLemb")

sl<-as.data.frame(sl)
sl[,"Year"] <- as.numeric(as.character(rownames(sl)))

#3. Harbor seals at Desdemona Sands, PV, Apr 1-Jun 30 ------
# combine surveys from odfw and wdfw
hs.odfw <- lcr %>% filter(month>3, month<7, spp=="PV", location=="COLUMBIA RIVER-DESDEMONA SANDS") %>% group_by(Year) %>%
    summarize(odfw=max(nonpup_total, na.rm=T))

wdfw <- read.csv(wdfwFile, header=T, stringsAsFactors=T)
hs.wdfw <- wdfw %>%  filter(Sitecode==1.04)  %>% #Desdemona Sands
    filter(Julian>90, Julian<182) %>% #apr 1-july 1
    group_by(Year) %>% summarize(wdfw=max(Estimated.total))

hs <- full_join(hs.odfw,hs.wdfw, by = "Year") %>% rowwise() %>% mutate(HS.DesSands=max(odfw,wdfw,na.rm=T)) %>%
    select(Year,HS.DesSands)

# Put all times series in one wide dataframe -------
pinnipedsNCC <- full_join(hs, sl, by="Year") #%>%  
pinnipedsNCC <- full_join(pinnipedsNCC, bon, by="Year") #%>% 
pinnipedsNCC <- arrange(pinnipedsNCC, Year)

# Standardize names
pinnipedsNCC <- pinnipedsNCC %>% rename(HS_DesSands=HS.DesSands, SSL_EMB=SSLemb, CSL_EMB=CSLemb, CSL_Bonn=CSLbonn, SSL_Bonn=SSLbonn)

inds <- names(pinnipedsNCC)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    thisIndicators <- indicators %>% filter(indicator==thisInd)
    thisOutputFile <- file.path(ColumbiaRiverPinnOutputDir, paste0(thisIndicators$shortName, ".csv"))
    cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
    cat("Region:,\n", file=thisOutputFile, append=T)
    cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
    cat("Contact:,\n", file=thisOutputFile, append=T)
    thisDat <- pinnipedsNCC[, c("Year", thisInd)]
    names(thisDat) <- c("Year", "Value")
    thisDat <- addMissingYears(thisDat)
    suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))
}

###########################################################################
# Extract and preprocess new pinniped and salmon data from Lisa
pinnipeds2025 <- read.csv(file.path(Lisa2025dir, "pinnipedsNCC.datwide.csv"))

inds <- names(pinnipeds2025)
inds <- inds[inds!="Year" & inds!="X"]
for(thisInd in inds) {
    thisIndicators <- indicators %>% filter(indicator==thisInd)
    thisOutputFile <- file.path(SEM2025outputDir, paste0(thisIndicators$shortName, ".csv"))
    cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
    cat("Region:,\n", file=thisOutputFile, append=T)
    cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
    cat("Contact:,\n", file=thisOutputFile, append=T)
    thisDat <- pinnipeds2025[, c("Year", thisInd)]
    names(thisDat) <- c("Year", "Value")
    thisDat <- addMissingYears(thisDat)
    suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))
}

# Create 2-year lead versions of Harbor_seal_CR_2025, AllSeaLionsEMB_2025 and AllSeaLionsBonn_2025
inds <- c("Harbor_seal_CR_2025", "AllSeaLionsEMB_2025", "AllSeaLionsBonn_2025")
for(thisInd in inds) {
    thisIndLead <- gsub("_2025", "_2yrLead_2025", thisInd)
    thisIndicators <- indicators %>% filter(indicator==thisIndLead)
    thisOutputFile <- file.path(SEM2025outputDir, paste0(thisIndicators$shortName, ".csv"))
    cat(paste0("Dataset:,", thisIndLead, "\n"), file=thisOutputFile)
    cat("Region:,\n", file=thisOutputFile, append=T)
    cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
    cat("Contact:,\n", file=thisOutputFile, append=T)
    
    thisDat <- pinnipeds2025[, c("Year", thisInd)]
    names(thisDat) <- c("Year", "Value")
    
    thisDat$Year <- thisDat$Year - 2
    
    thisDat <- addMissingYears(thisDat)
    suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))
}

salmon2025 <- read.csv(file.path(Lisa2025dir, "salmon.datwide.csv"))
inds <- names(salmon2025)
inds <- inds[inds!="Year" & inds!="X"]
for(thisInd in inds) {
    thisIndicators <- indicators %>% filter(indicator==thisInd)
    thisOutputFile <- file.path(SEM2025outputDir, paste0(thisIndicators$shortName, ".csv"))
    cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
    cat("Region:,\n", file=thisOutputFile, append=T)
    cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
    cat("Contact:,\n", file=thisOutputFile, append=T)
    thisDat <- salmon2025[, c("Year", thisInd)]
    names(thisDat) <- c("Year", "Value")
    thisDat <- addMissingYears(thisDat)
    suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))
}

###########################################################################
# Extract and preprocess new versions of SEM 2024 data
hake2025 <- read.csv(file.path(SEM2025dir, "Age1Hake2024Assessment.csv"))
inds <- names(hake2025)
inds <- inds[inds!="Year"]
inds <- paste0(inds, "_2025")
names(hake2025) <- c("Year", inds)
for(thisInd in inds) {
    thisIndicators <- indicators %>% filter(indicator==thisInd)
    thisOutputFile <- file.path(SEM2025outputDir, paste0(thisIndicators$shortName, ".csv"))
    cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
    cat("Region:,\n", file=thisOutputFile, append=T)
    cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
    cat("Contact:,\n", file=thisOutputFile, append=T)
    thisDat <- hake2025[, c("Year", thisInd)]
    names(thisDat) <- c("Year", "Value")
    thisDat <- addMissingYears(thisDat)
    suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))
}

allSEM2025 <- read.csv(file.path(SEM2025dir, "AllDataNCC.csv"), skip=2)
names(allSEM2025) <- paste0(names(allSEM2025), "_2025")
allSEM2025 <- allSEM2025 %>% rename(Year="Metric_2025") %>% filter(!is.na(Year))

# Fill in killer whale data to match SEM_2024 version
allSEM2025$Killer.whales.BC_2025[is.na(allSEM2025$Killer.whales.BC_2025)] <- 261

# Create 2-year lead version of killer whale data
allSEM2025$Killer.whales.BC_2yrLead_2025 <- c(allSEM2025$Killer.whales.BC_2025[3:(nrow(allSEM2025))], NA, NA)

inds <- names(allSEM2025)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    if(thisInd %in% indicators$indicator) {
        thisIndicators <- indicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(SEM2025outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- allSEM2025[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
    }
}

###########################################################################
# Extract and preprocess new data, 2025
blackRockfish_2025 <- read_excel(blackRockfishFile, range="B4:J99")
blackRockfish_2025 <- blackRockfish_2025 %>% select(Year, `Total Biomass (mt)`) %>% rename(Value="Total Biomass (mt)")

inds <- c("blackRockfish_2025")
dfs <- list(blackRockfish_2025)
for(i in 1:length(inds)) {
    thisInd <- inds[i]
    thisIndicators <- indicators %>% filter(indicator==thisInd)
    thisOutputFile <- file.path(data2025outputDir, paste0(thisIndicators$shortName, ".csv"))
    cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
    cat("Region:,,\n", file=thisOutputFile, append=T)
    cat("Source:,", thisIndicators$source, "\n", file=thisOutputFile, append=T)
    cat("Contact:,\n", file=thisOutputFile, append=T)
    thisDat <- addMissingYears(dfs[[i]])
    suppressWarnings(write.table(dfs[[i]], file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))    
}

FRAM <- read.csv(FRAMfile)
for(thisInd in unique(FRAM$species)) {
    thisDat <- FRAM %>% filter(species==thisInd) %>% rename(Value=annualCPUE) %>% select(Year, Value)
    thisIndicators <- indicators %>% filter(indicator==thisInd)
    thisOutputFile <- file.path(data2025outputDir, paste0(thisIndicators$shortName, ".csv"))
    cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
    cat("Region:,,\n", file=thisOutputFile, append=T)
    cat("Source:,", thisIndicators$source, "\n", file=thisOutputFile, append=T)
    cat("Contact:,\n", file=thisOutputFile, append=T)
    thisDat <- addMissingYears(thisDat)
    suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))  
}

JSOESjellies <- read.csv(JSOESjelliesFile)

inds <- names(JSOESjellies)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    if(thisInd %in% indicators$indicator) {
        thisIndicators <- indicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(data2025outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- JSOESjellies[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
    }
}

bongo <- read.csv(bongoFile)

inds <- names(bongo)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    if(thisInd %in% indicators$indicator) {
        thisIndicators <- indicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(data2025outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- bongo[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
    }   
}

CPS <- read.csv(CPSfile)

inds <- unique(CPS$species)
for(thisInd in inds) {
    if(thisInd %in% indicators$indicator) {
        thisIndicators <- indicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(data2025outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- CPS %>% filter(species==thisInd) %>% rename(Year=year, Value=meanProbOccur) %>% select(Year, Value)
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
    }   
}

hakeASHOP <- read.csv(hakeASHOPfile)
inds <- names(hakeASHOP)
inds <- inds[inds!="year"]
for(thisInd in inds) {
    if(thisInd %in% indicators$indicator) {
        thisIndicators <- indicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(data2025outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- hakeASHOP[, c("year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
    }   
}

###########################################################################
# Extract and preprocess NHL data, 2025
planktonNHL <- read.csv(planktonNHLfile)
inds <- names(planktonNHL)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    if(thisInd %in% indicators$indicator) {
        thisIndicators <- indicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(NHL2025outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- planktonNHL[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))  
    }
}

###########################################################################
# Extract and preprocess planktonJuneNCC data, 2025
planktonJuneNCC <- read.csv(planktonJuneNCCfile)
inds <- names(planktonJuneNCC)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    if(thisInd %in% indicators$indicator) {
        thisIndicators <- indicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(planktonJuneNCCoutputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- planktonJuneNCC[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))  
    }
}

###########################################################################
# Extract and preprocess plankton_2026 data
plankton2026 <- read.csv(plankton2026file)
inds <- names(plankton2026)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    if(thisInd %in% indicators$indicator) {
        thisIndicators <- indicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(plankton2026outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- plankton2026[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))  
    }
}

###########################################################################
# Extract and preprocess krill NHL data, 2026
krill <- read.csv(krillFile)
inds <- names(krill)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    if(thisInd %in% indicators$indicator) {
        thisIndicators <- indicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(NHL2026outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- krill[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
    }   
}

###########################################################################
# Extract and preprocess JSOES data, 2026
bird <- read.csv(birdFile)
inds <- names(bird)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    thisDatasetIndicators <- indicators %>% filter(dataset=="JSOES_2026")
    if(thisInd %in% thisDatasetIndicators$indicator) {
        thisIndicators <- thisDatasetIndicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(JSOES2026outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- bird[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
    } 
}

###########################################################################
# Extract and preprocess Westport survey data, 2026
for(f in c(WS2026file, WS2026groupedFile)) {
    WS <- read.csv(f)
    inds <- names(WS)
    inds <- inds[inds!="year"]
    for(thisInd in inds) {
        thisDatasetIndicators <- indicators %>% filter(dataset=="Westport_2026")
        if(thisInd %in% thisDatasetIndicators$indicator) {
            thisIndicators <- thisDatasetIndicators %>% filter(indicator==thisInd)
            thisOutputFile <- file.path(WS2026outputDir, paste0(thisIndicators$shortName, ".csv"))
            cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
            cat("Region:,\n", file=thisOutputFile, append=T)
            cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
            cat("Contact:,\n", file=thisOutputFile, append=T)
            thisDat <- WS[, c("year", thisInd)]
            names(thisDat) <- c("Year", "Value")
            thisDat <- addMissingYears(thisDat)
            suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
        } 
    } 
}

# Create 2-year lead versions of Californian_s_l_WS_2026, Northern_f_s_WS_2026, Harbour_s_WS_2026,
# and Steller_s_l_WS_2026
WS <- read.csv(WS2026file)
inds <- c("Californian_s_l_WS_2026", "Northern_f_s_WS_2026", "Harbour_s_WS_2026", "Steller_s_l_WS_2026")
for(thisInd in inds) {
    thisIndLead <- gsub("WS_2026", "2yrLead_WS_2026", thisInd)
    thisIndicators <- indicators %>% filter(indicator==thisIndLead)
    thisOutputFile <- file.path(WS2026outputDir, paste0(thisIndicators$shortName, ".csv"))
    cat(paste0("Dataset:,", thisIndLead, "\n"), file=thisOutputFile)
    cat("Region:,\n", file=thisOutputFile, append=T)
    cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
    cat("Contact:,\n", file=thisOutputFile, append=T)
    
    thisDat <- WS[, c("year", thisInd)]
    names(thisDat) <- c("Year", "Value")
    
    thisDat$Year <- thisDat$Year - 2
    
    thisDat <- addMissingYears(thisDat)
    suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
}


###########################################################################
# Extract and preprocess WGOA DFA data, 2026
WGOA_DFA_2026 <- read.csv(WGOA_DFA_2026file)
inds <- names(WGOA_DFA_2026)
inds <- inds[inds!="year"]
for(thisInd in inds) {
    thisDatasetIndicators <- indicators %>% filter(dataset=="WGOA_DFA_2026")
    if(thisInd %in% thisDatasetIndicators$indicator) {
        thisIndicators <- thisDatasetIndicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(WGOA_DFA_2026outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- WGOA_DFA_2026[, c("year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
    } 
}

###########################################################################
# Extract and preprocess lingcod stock assessment data, 2026
lingcodSA_2026 <- read.csv(lingcodSA_2026file)
inds <- names(lingcodSA_2026)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    thisDatasetIndicators <- indicators %>% filter(dataset=="lingcodSA_2026")
    if(thisInd %in% thisDatasetIndicators$indicator) {
        thisIndicators <- thisDatasetIndicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(lingcodSA_2026outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- lingcodSA_2026[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
    } 
}

###########################################################################
# Extract and process AK predator data, 2026
predAK_2026 <- read.csv(predAK_2026file)
inds <- names(predAK_2026)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    thisDatasetIndicators <- indicators %>% filter(dataset=="predAK_2026")
    if(thisInd %in% thisDatasetIndicators$indicator) {
        thisIndicators <- thisDatasetIndicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(predAK_2026outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- predAK_2026[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
    } 
}

###########################################################################
# Orca data, 2026
orca_2026 <- read.csv(orca_2026file)

# Apply 2-year lead
orca_2026$YEAR <- orca_2026$YEAR-2

orca_2026 <- orca_2026 %>% rename(Year=YEAR, Killer.whales.NR.BC_2yrLead_2026=Size_best) %>% 
    select(Year, Killer.whales.NR.BC_2yrLead_2026)

inds <- names(orca_2026)
inds <- inds[inds!="Year"]
for(thisInd in inds) {
    thisDatasetIndicators <- indicators %>% filter(dataset=="orca_2026")
    if(thisInd %in% thisDatasetIndicators$indicator) {
        thisIndicators <- thisDatasetIndicators %>% filter(indicator==thisInd)
        thisOutputFile <- file.path(orca_2026outputDir, paste0(thisIndicators$shortName, ".csv"))
        cat(paste0("Dataset:,", thisInd, "\n"), file=thisOutputFile)
        cat("Region:,\n", file=thisOutputFile, append=T)
        cat(paste0("Source:,", thisIndicators$source, "\n"), file=thisOutputFile, append=T)
        cat("Contact:,\n", file=thisOutputFile, append=T)
        thisDat <- orca_2026[, c("Year", thisInd)]
        names(thisDat) <- c("Year", "Value")
        thisDat <- addMissingYears(thisDat)
        suppressWarnings(write.table(thisDat, file=thisOutputFile, sep=",", row.names=F, quote=F, append=T))   
    } 
}



