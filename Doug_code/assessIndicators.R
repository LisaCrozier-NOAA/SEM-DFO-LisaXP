# Assess indicators to see if and how they should be included in the analysis.
# 
# Doug Jackson
# doug@QEDAconsulting.com
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(imputeTS)
###########################################################################
# Constants
###########################################################################
workingDir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices"
indicatorsFile <- file.path(workingDir, "indicators.csv")

figWidth <- 6
figHeight <- 4

maxLag <- 10

###########################################################################
# Run
###########################################################################
setwd(workingDir)

source("functions.R")

backgroundTheme <- theme(panel.background=element_rect(fill="white"),
                         plot.background=element_rect(fill="white"))

dataDir <- file.path(workingDir, "data")
outputDir <-  file.path(workingDir, "output_assessIndicators")
dir.create(outputDir, showWarnings=F, recursive=T)
outputTimeSeriesDir <- file.path(workingDir, "output_assessIndicators", "outputTimeseries")
dir.create(outputTimeSeriesDir, showWarnings=F, recursive=T)

# Only include indicators with an associated category
indicators <- read.csv(indicatorsFile)
indicators <- indicators[indicators$category!="", ]

# Create a data frame containing all the indicators in each subdirectory
dataFileIndicatorList <- list()
for(subDir in subDirs) {
    thisDataFiles <- list.files(file.path(dataDir, subDir), pattern="\\.csv")
    thisDataset <- subDir
    
    for(dataFile in thisDataFiles) {
        thisFile <- file.path(dataDir, subDir, dataFile)
        thisHeader <- readLines(thisFile, n=4)
        thisIndicator <- trimws(unlist(strsplit(thisHeader[[1]], ","))[2])
        thisIndicator <- gsub('"', "", thisIndicator)
        
        thisDataFileIndicator <- data.frame(indicator=thisIndicator, dataset=subDir, file=basename(thisFile))
        dataFileIndicatorList[[length(dataFileIndicatorList)+1]] <- thisDataFileIndicator
    }
}
dataFileIndicator <- bind_rows(dataFileIndicatorList)
dataFileIndicator <- dataFileIndicator %>% arrange(indicator, dataset)
write.csv(dataFileIndicator, file=file.path(outputDir, "dataFileIndicator.csv"), row.names=F)

# Datetimes defining screening period
screenStartDatetime <- dmy(screenStartDate)
screenEndDatetime <- dmy(screenEndDate)
yearDF <- data.frame(year=year(seq(screenStartDatetime, screenEndDatetime, by="1 year")), ID=1)

fracCompleteList <- list()

# Read and plot data
for(subDir in subDirs) {
    thisDataFiles <- list.files(file.path(dataDir, subDir), pattern="\\.csv")
    thisDataset <- subDir
    
    for(dataFile in thisDataFiles) {
        thisFile <- file.path(dataDir, subDir, dataFile)
        thisHeader <- readLines(thisFile, n=4)
        thisIndicator <- trimws(unlist(strsplit(thisHeader[[1]], ","))[2])
        thisIndicator <- gsub('"', "", thisIndicator)
        
        if(thisIndicator %in% indicators$indicator) {
        
            cat("Processing", dataFile, "-", thisIndicator, "\n")
            fileName <- gsub(".csv", "", basename(thisFile))
            
            thisData <- read.csv(thisFile, skip=4, colClasses=c("character", "numeric"), na.strings=c("null", "NA"))
            
            # Convert time column to datetime
            if(names(thisData)[1]=="Year") {
                thisData$date <- ymd(thisData$Year, truncated=2L)
                thisData <- thisData[, c("date", "Value")]
            } else {
                thisData$date <- ymd(thisData$date)
            }
            names(thisData) <- c("date", "value")
            
            # Impute
            thisInd <- getInd(indicators, thisIndicator, thisDataset)
            thisData <- impute(thisInd, thisData)
            
            thisShortName <- thisInd$shortName
            
            thisScreenDat <- thisData %>% filter(date>=screenStartDatetime & date<=screenEndDatetime)
            
            # Check minimum completeness criterion
            # Note: this assumes that data are yearly
            thisScreenDat$year <- year(thisScreenDat$date)
            thisYearDat <- thisScreenDat %>% group_by(year) %>% summarize(value=mean(value, na.rm=T), imputed=mean(imputed, na.rm=T))
            thisYearDat <- full_join(thisYearDat, yearDF, by="year") %>% arrange(year)
            fracComplete <- sum(is.finite(thisYearDat$imputed))/nrow(thisYearDat)
            
            fracCompleteList[[length(fracCompleteList)+1]] <- data.frame(subDir=subDir, dataFile=dataFile,
                                                                         shortName=thisInd$shortName, fracComplete=fracComplete,
                                                                         minYear=min(thisData$date, na.rm=T), maxYear=max(thisData$date, na.rm=T))
            
            endY <- 0.05*max(thisData$value, na.rm=T)
            p <- ggplot() + geom_point(data=thisData, aes(x=date, y=value)) +
                geom_line(data=thisData, aes(x=date, y=value)) +
                geom_point(data=thisData, aes(x=date, y=imputed), color="red", size=0.5) + 
                geom_line(data=thisData, aes(x=date, y=imputed), color="red", linewidth=0.25) +
                geom_segment(data=thisData[is.na(thisData$value), ],
                             aes(x=date, xend=date, y=0, yend=endY), color="red", linewidth=0.25) +
                labs(title=thisIndicator, subtitle=subDir) +
                theme_light()
            
            # Calculate autocorrelation
            # See https://aosmith.rbind.io/2018/06/27/uneven-grouped-autocorrelation/
            png(filename=file.path(outputTimeSeriesDir, "scratch.png"))
            acfOut <- acf(thisData$value, na.action=na.pass, main=paste0(subDir, ": ", thisIndicator),
                          ylim=c(-2, 2))
            acfDF <- with(acfOut, data.frame(lag, acf))

            # Calculate confidence intervals given different number of occurrences of different lags
            thisData$row <- 1:nrow(thisData)
            thisDataNoNA <- thisData[complete.cases(thisData), ]
            N <- c()
            for(thisLag in 1:maxLag) {
                thisN <- 0
                for(innerLag in 1:maxLag) {
                    thisN <- thisN + sum(diff(thisDataNoNA$row, lag=innerLag)==thisLag)
                }
                N <- c(N, thisN)
            }
            lags <- data.frame(lag=1:maxLag, N=N)
            lags$lowerCI <- -qnorm(1-.025)/sqrt(lags$N)
            lags$upperCI <- qnorm(1-.025)/sqrt(lags$N)
            lines(1:maxLag, lags$lowerCI, lty = 2)
            lines(1:maxLag, lags$upperCI, lty = 2)
            dev.off()

            # Create a ggplot ACF plot
            pACF <- ggplot() + geom_bar(data=acfDF, aes(x=lag, y=acf), stat="identity", position="identity") +
                geom_line(data=lags, aes(x=lag, y=lowerCI), linetype="dashed") +
                geom_line(data=lags, aes(x=lag, y=upperCI), linetype="dashed") +
                theme_light() + backgroundTheme

            g <- arrangeGrob(grobs=list(p, pACF), nrow=2)
            ggsave(filename=file.path(outputTimeSeriesDir, paste0(thisDataset, "_", thisShortName, "_ts.png")), g,
                   width=8, height=8)
        }
    }
}

fracCompleteDF <- bind_rows(fracCompleteList)
write.csv(fracCompleteDF, file=file.path(outputTimeSeriesDir, "fracComplete.csv"), row.names=F)