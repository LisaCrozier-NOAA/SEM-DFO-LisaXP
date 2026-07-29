# Explore and analyze various indicators:
# Alaskan indicators from https://apps-afsc.fisheries.noaa.gov/refm/reem/ecoweb/Index.php?ID=9
# CCIEA indicators from https://www.integratedecosystemassessment.noaa.gov/regions/california-current/california-current-iea-indicators
# WCVI indicators from DFO
# SEM indicators from separate SEM analysis led by Noble
# 
# Doug Jackson
# doug@QEDAconsulting.com
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(fpc)
library(imputeTS)
library(dendextend)
library(corrplot)
library(ggcorrplot)
library(htmlwidgets)
library(DT)
library(dendextend)
library(flextable)
library(scales)
library(stringr)
library(MARSS)
###########################################################################
# Constants
###########################################################################
workingDir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices"
indicatorsFile <- file.path(workingDir, "indicators.csv")

scens <- c("incl2026")

# Optionally specify the number of clusters for each scenario
numClusters <- list() #list("clusCond1"=3, "clusPred"=6)

# Group time series according to pamk clusters (TRUE) or hierarchical clusters (FALSE)
timeSeriesGroupsPAMK <- TRUE

# Minimum year at which to truncate time series
minYears <- c(1998) #c(1990, 2000)

# Maximum year at which to truncate time series, e.g., to match SAR availability
maxYear <- 2021

SARfileName <- "sar.19982021.csv"

# Minimum fraction of coverage (not missing) during screening period
# Reduce minFracComplete from 0.79 to 0.5 to include Harbor_seal_CR_2025
minFracComplete <- 0.5 #0.79

figWidth <- 6
figHeight <- 4

maxLag <- 10

legendTextSize <- 6

# Base font size factor for correlations
baseFontSize <- 60

###########################################################################
# Functions
###########################################################################
clusterAnalysis <- function(dataWide, scen, plotTitle) {
    dataInd <- dataWide %>% select(-date)
    dataIndStd <- dataInd %>% mutate_all(~(scale(.) %>% as.vector))
    
    # Plot standardized values
    plotData <- dataIndStd
    plotData$date <- dataWide$date
    plotData <- plotData %>% pivot_longer(cols=-date, names_to="shortName")
    
    # Drop rows with all missing data
    notMissingAll <- rowSums(is.na(dataIndStd))!=ncol(dataIndStd)
    dataIndStd <- dataIndStd[notMissingAll, ]
    
    dataIndStdTrans <- t(dataIndStd)
    
    # Perform partitioning around medoids clustering to determine cluster assignments.
    # Optionally prescribe the number of clusters
    if(scen %in% names(numClusters)) {
        kRange <- numClusters[[scen]]:numClusters[[scen]]
        plotTitle <- paste0(plotTitle, ", prescribed number of clusters: ", numClusters[[scen]])
    } else {
        kRange <- 2:(nrow(dataIndStdTrans)-1)
    }
    dataPAM <- pamk(dataIndStdTrans, krange=kRange)
    dataPAMclusters <- data.frame(dataPAM$pamobject$clustering) %>% rownames_to_column() %>% 
        rename(shortName=rowname, group=dataPAM.pamobject.clustering)
    
    # Perform hierarchical clustering for dendrograms
    dataIndDist <- dist(dataIndStdTrans, method="euclidean")
    dataFit <- hclust(dataIndDist, method="ward.D")

    png(file=file.path(outputDir, paste0(scen, "Clusters.png")), width=8, height=8, units="in", res=300)
    par(mar=c(5 + 5,4,4,2) + 0.1)
    dend <- as.dendrogram(dataFit) %>% color_branches(k=dataPAM$nc) %>% set("labels_colors", k=dataPAM$nc) %>% 
        set("branches_lwd", 3) %>% plot(main=plotTitle)
    dev.off()
    par(mar=c(5,4,4,2) + 0.1)
    
    # For time series, assign groups according to either pamk or hierarchical clustering
    if(timeSeriesGroupsPAMK) {
        dataGroups <- dataPAMclusters
        plotTitle <- paste0(plotTitle, ", pamk clusters")
    } else {
        dataGroups <- data.frame(cutree(dataFit, k=dataPAM$nc)) %>% rownames_to_column()
        names(dataGroups) <- c("shortName", "group")
        plotTitle <- paste0(plotTitle, ", hierarchical clusters")
    }
    
    # Plot standardized values, with separate plots by group
    plotData <- left_join(plotData, dataGroups, by="shortName") %>% arrange(group, shortName, date)
    p <- ggplot(plotData) + geom_point(aes(x=date, y=value, color=shortName, shape=shortName, fill=shortName, group=shortName)) + 
        geom_line(aes(x=date, y=value, color=shortName, group=shortName)) + 
        scale_shape_manual(values=rep(21:25, 100)) +
        facet_grid(group~.) +
        labs(title=plotTitle, x="", y="") +
        theme_light() + theme(legend.text=element_text(size=10)) + backgroundTheme
    ggsave(filename=file.path(outputDir, paste0(scen, "Std.png")), width=10, height=6)
    
    return(dataPAMclusters)
}

plotCor <- function(dat, filename, plotTitle, figWidth, figHeight, useShortName=F, createPlot=T) {
    colNamesDF <- dat %>% select(shortName, plotName) %>% distinct() %>% arrange(plotName)
    if(useShortName) {
        dat <- dat %>% mutate(plotName=shortName)
        indicators <- indicators %>% mutate(plotName=shortName)
        filename <- paste0(filename, "_shortName")
        colNamesDF$plotName <- colNamesDF$shortName
    }
    colNames <- colNamesDF$plotName
    
    # Merge in trophic group colors
    colNamesDF <- left_join(colNamesDF, indicators[, c("plotName", "trophicGroupCol")], by="plotName")
    
    dat <- dat %>% select(plotName, date, finalVal) %>% 
        pivot_wider(id_cols=date, names_from=plotName, values_from=finalVal) %>% select(-date) %>% mutate_all(~(scale(.) %>% as.vector))
    
    dat <- dat[, colNames]
    datCor <- cor(dat, use="pairwise.complete.obs", method="spearman")
    datRes <- cor.mtest(dat, conf.level=0.95)
    
    if(createPlot) {
        png(file.path(outputDir, paste0(filename, ".png")), width=figWidth, height=figHeight, units="in", res=300)
        corrplot(datCor, p.mat=datRes$p, method="circle", type="lower", insig="blank",
                 order="original", diag=FALSE, title=plotTitle, mar=c(0, 0, 1, 0))$corrPos -> p1
        text(p1$x, p1$y, round(p1$corr, 2))
        dev.off()        
    }
    
    # Roll my own
    pMat <- cor_pmat(dat)
    datCorLong <- datCor  %>% as_tibble(rownames="varA")  %>% pivot_longer(-varA, names_to="varB", values_to="correlation")
    varNames <- unique(datCorLong$varA)
    datCorLong <- datCorLong %>% mutate(varA=factor(varA, levels=varNames), varB=factor(varB, levels=rev(varNames)))
    # Set redundant upper half to NA
    datCorLong <- datCorLong %>% mutate(lvlA=as.numeric(varA), lvlB=as.numeric(fct_rev(varB)), correlation=ifelse(lvlA<lvlB, correlation, NA))
    pMatLong <- pMat %>% as_tibble(rownames="varA")  %>% pivot_longer(-varA, names_to="varB", values_to="p")
    datCorLong <- full_join(datCorLong, pMatLong, by=c("varA", "varB"))
    datCorLong <- datCorLong %>% mutate(corrForCol=ifelse(p<0.05, correlation, NA))
    
    # Color axis labels according to trophicGroup
    trophicGroupColB <- colNamesDF$trophicGroupCol
    trophicGroupColA <- rev(trophicGroupColB)
    
    datCorLong <- datCorLong %>% mutate(varA=factor(varA, levels=varNames), varB=factor(varB, levels=rev(varNames))) %>% 
        arrange(varA, varB)
    
    if(createPlot) {
        p <- ggplot(data=datCorLong, aes(varA, varB)) + geom_tile(aes(fill=corrForCol), color="black") +
            scale_fill_gradient2(high="dodgerblue4", mid="white", low="firebrick2", limits=c(-1, 1), midpoint=0, na.value="white") +
            geom_text(aes(label=round(correlation, 2)), size=baseFontSize/ncol(dat)) +
            theme_minimal(base_size=16) +
            labs(x=element_blank(), y=element_blank(), fill="corr", title=plotTitle) +
            theme(axis.text.x=ggtext::element_markdown(angle=45, vjust=1, hjust=1, color=trophicGroupColB, face="bold", size=10), 
                  axis.text.y=ggtext::element_markdown(color=trophicGroupColA, face="bold", size=10), legend.position="none") +
            backgroundTheme
        ggsave(file.path(outputDir, paste0("customCorr_", filename, ".png")), width=figWidth, height=figHeight)        
    }

    return(datCorLong)
}

# Check if data set meets various screening criteria
checkCriteria <- function(dat) {
    # Datetimes defining screening period
    screenStartDatetime <- dmy(screenStartDate)
    screenEndDatetime <- dmy(screenEndDate)
    
    yearDF <- data.frame(year=year(seq(screenStartDatetime, screenEndDatetime, by="1 year")), ID=1)
    
    meetsCriteriaDFlist <- list()
    for(thisShortName in unique(dat$shortName)) {
        thisDat <- dat %>% filter(shortName==thisShortName)
        
        thisScreenDat <- thisDat %>% filter(date>=screenStartDatetime & date<=screenEndDatetime)
        
        # Check minimum completeness criterion
        # Note: this assumes that data are yearly
        thisScreenDat$year <- year(thisScreenDat$date)
        thisYearDat <- thisScreenDat %>% group_by(year) %>% summarize(value=mean(value, na.rm=T), finalVal=mean(finalVal, na.rm=T))
        thisYearDat <- full_join(thisYearDat, yearDF, by="year") %>% arrange(year)
        fracComplete <- sum(is.finite(thisYearDat$finalVal))/nrow(thisYearDat)
        passFracComplete <- fracComplete>=minFracComplete
        
        # Don't apply fraction complete criterion to SAR data
        if(str_starts(toupper(thisShortName), "SAR_")) {
            passFracComplete <- TRUE
            cat(thisShortName, " is an SAR variable => automatically passed fracComplete criterion.\n")
        }
        
        if(!passFracComplete) {cat(thisShortName, " failed fracComplete criterion.\n")}

        # Logical AND of all criteria
        meetsCriteria <- (passFracComplete)
        dat[dat$shortName==thisShortName, "meetsCriteria"] <- meetsCriteria
        
        meetsCriteriaDFlist[[length(meetsCriteriaDFlist)+1]] <- data.frame(shortName=thisShortName, 
                                                                           fracComplete=fracComplete, passFracComplete=passFracComplete, 
                                                                           meetsCriteria=meetsCriteria)
    }
    meetsCriteriaDF <- bind_rows(meetsCriteriaDFlist)
    
    return(list(dat=dat, meetsCriteriaDF=meetsCriteriaDF))
}

colorCell <- JS(
    "function(data, type, row, meta) {",
    "  if(type === 'display') {",
    "    return '<span style=\"color: ' + data + '\">' + data + '</span>';",
    "  } else {",
    "    return data;",
    "  }",
    "}"
)
###########################################################################
# Run
###########################################################################
setwd(workingDir)

source("functions.R")

backgroundTheme <- theme(panel.background=element_rect(fill="white"),
                         plot.background=element_rect(fill="white"))

dataDir <- file.path(workingDir, "data")

# Load smolt-to-adult return ratio
SAR <- read.csv(file.path(workingDir, "data", SARfileName))
SAR$date <- ymd(SAR$Year, truncated=2L)
SAR <- SAR[, c("date", "SAR")]

clusDataList <- list()
corrDataList <- list()
for(scen in scens) {
    for(minYear in minYears) {
        cat("===================================================================\n")
        cat(paste0("Running ", scen, ", minYear = ", minYear, "\n"))
        
        outputDir <- file.path(workingDir, "output", paste0(scen, "_", minYear))
        dir.create(outputDir, showWarnings=F, recursive=T)
        
        scratchDir <- file.path(workingDir, "scratch")
        dir.create(scratchDir, showWarnings=F, recursive=T)
        
        # Only include indicators that have an associated category
        indicators <- read.csv(indicatorsFile)
        indicators <- indicators[indicators$category!="", ]
        
        # Only include specified indicators
        indicators$include <- indicators[, scen]
        indicators <- indicators[indicators$include=="Y", ]
        
        # Create plotName from region_subregion_trophicGroup_ID
        indicators$plotName <- paste(indicators$region, indicators$subregion, indicators$trophicGroup, indicators$ID, sep="_")
        
        # Create colors for correlation matrix axis labels
        trophicGroupCols <- read.csv(file.path(workingDir, "trophicGroupCols.csv"))
        indicators <- left_join(indicators, trophicGroupCols, by="trophicGroup")
        
        # Read in data
        dataList <- list()
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
                    if(names(thisData)[2]=="Index") {
                        thisData$date <- ymd(thisData$Year, truncated=2L)
                        thisData <- thisData[, c("date", "Index")]
                    } else if(names(thisData)[1]=="Year") {
                        thisData$date <- ymd(thisData$Year, truncated=2L)
                        thisData <- thisData[, c("date", "Value")]
                    } else {
                        thisData$date <- ymd(thisData$date)
                    }
                    names(thisData) <- c("date", "value")
                    
                    # Impute
                    thisInd <- getInd(indicators, thisIndicator, thisDataset)
                    thisData <- impute(thisInd, thisData)
                    
                    # Log transform
                    out <- logTransform(thisInd, thisData)
                    thisData <- out$thisData
                    
                    # Use untransformed or transformed values as final values, as appropriate
                    if(out$transformed) {
                        thisData$finalVal <- thisData$logTransformed
                    } else {
                        thisData$finalVal <- thisData$imputed
                    }
                    
                    # Create smoothed version
                    thisVal <- thisData %>% select(finalVal) %>% mutate_all(~(scale(.) %>% as.vector))
                    thisVal <- t(thisVal)
                    
                    # Only smooth between the first year with data and the last year with data
                    notNAindices <- which(!is.na(thisVal))
                    firstYearIndex <- notNAindices[[1]]
                    lastYearIndex <- notNAindices[[length(notNAindices)]]
                    
                    fit <- MARSS(thisVal[firstYearIndex:lastYearIndex], fit=F)
                    fit$par <- fit$start
                    kfList <- MARSSkf(fit)
                    
                    smoothedDat <- thisVal
                    smoothedDat[firstYearIndex:lastYearIndex] <- as.numeric(t(kfList$xtT))
                    thisData$smoothed <- as.numeric(t(smoothedDat))

                    thisData$shortName <- thisInd$shortName
                    thisData$plotName <- thisInd$plotName
                    
                    thisData$dataset <- subDir
                    thisData$indicator <- thisIndicator
                    thisData$category <- thisInd$category
                    thisData$SEMlatent <- thisInd$SEMlatent
                    
                    thisData <- thisData[, c("indicator", "dataset", "shortName", "plotName", "SEMlatent", "date", 
                                             "value", "imputed", "finalVal", "smoothed")]
                    
                    dataList[[length(dataList)+1]] <- thisData
                } else {
                    cat(thisIndicator, "not found in indicators.csv. subDir =", subDir, ", dataFile = ", dataFile, "\n")
                }
            }
        }
        data <- bind_rows(dataList)

        saveRDS(data, file=file.path(workingDir, "output", "allData.rds"))
        
        # Drop years outside of specified minYear-maxYear range
        data <- data %>% filter(date>=ymd(minYear, truncated=2L), date<=ymd(maxYear, truncated=2L))
        
        # Check screening criteria
        output <- checkCriteria(data)
        data <- output$dat %>% filter(meetsCriteria)
        write.csv(output$meetsCriteriaDF, file.path(outputDir, "meetsCriteria.csv"), row.names=F)
        
        # Only include years for which we have SAR data
        # data <- data %>% filter(date %in% SAR$date)
        
        # Create html table of indicators
        indicatorsTable <- indicators %>% filter(indicator %in% data$indicator) %>% arrange(ID) %>% select(plotName, shortName, indicator, trophicGroup)
        indicatorsDT <- datatable(indicatorsTable, rownames=F) %>% formatStyle("trophicGroup", target="row", 
                                                                         color=styleEqual(trophicGroupCols$trophicGroup,
                                                                                                    trophicGroupCols$trophicGroupCol))
        saveWidget(indicatorsDT, file.path(outputDir, paste0("indicators_", scen, "_", minYear, ".html")))
        
        # Create wide version of data frame so I can easily see if years overlap
        data <- data %>% arrange(indicator, dataset, date)
        datWide <- data %>% pivot_wider(id_cols=c("date"), names_from=c("shortName"), values_from="finalVal") %>% arrange(date)
        write.csv(datWide, file.path(outputDir, "datWide.csv"), row.names=F)
        
        # Can't run clustering or correlation analysis unless there are multiple indicators
        if(ncol(datWide)<3) {
            cat("Only one indicator => skipping clustering and correlation analysis\n")
            
            if(startsWith(scen, "clus")) {
                # Save the single indicator as a single cluster
                thisClusData <- data %>% mutate(meetsCriteria=TRUE, group=1, scen=scen, minYear=minYear)
                clusDataList[[length(clusDataList)+1]] <- thisClusData
            }
            next
        }
        
        plotTitle <- paste0("minYear = ", minYear, ", ", scen)
        ###################################
        # Clustering analysis
        if(startsWith(scen, "clus") | startsWith(scen, "incl")) {
            dataPAMclusters <- clusterAnalysis(datWide, scen, plotTitle)
            data <- left_join(data, dataPAMclusters, by="shortName")
            
            thisClusData <- data %>% mutate(scen=scen, minYear=minYear)
            clusDataList[[length(clusDataList)+1]] <- thisClusData
            
            # Save table of clusters
            thisClusters <- thisClusData %>% rename(dendClusID=group) %>% select(SEMlatent, dendClusID, shortName) %>% distinct() %>% 
                arrange(SEMlatent, dendClusID, shortName) %>% select(dendClusID, shortName) %>% rename(cluster=dendClusID) %>% 
                mutate(cluster=as.factor(cluster)) %>% arrange(cluster, shortName)
            
            colorer <- col_factor("RdYlBu", domain=NULL)
            ft <- flextable(thisClusters) %>% bg(j="cluster", bg=colorer, part="body")
            save_as_image(ft, file.path(outputDir, paste0(scen, "ClustersTable.png")))
            
            # Correlation analysis within clusters
            thisGroups <- unique(thisClusData$group)
            for(thisGroup in thisGroups) {
                thisShortNames <- thisClusData %>% filter(group==thisGroup) %>% distinct(shortName)
                
                if(nrow(thisShortNames)>1) {
                    thisDat <- data %>% filter(shortName %in% thisShortNames$shortName)
                    thisCor <- plotCor(thisDat, "scratch", "", -999, -999, useShortName=T, createPlot=F)
                    thisCor <- thisCor %>% filter(!is.na(correlation)) %>% select(varA, varB, correlation, p)
                    thisCor$SEMnode <- tolower(gsub("clus", "", scen))
                    thisCor$group <- thisGroup
                    corrDataList[[length(corrDataList)+1]] <- thisCor
                }
            }
        }
        ###################################
        # Correlation analysis
        if(startsWith(scen, "corr")) {
            if(scen=="corrPreyComp") {
                figSize <- 20
            } else {
                figSize <- 8
            }
            scratch <- plotCor(data, scen, plotTitle, figWidth=figSize, figHeight=figSize)
            
            # Create another set of correlation plots using shortName as the plot name
            datCorLong <- plotCor(data, scen, plotTitle, figWidth=figSize, figHeight=figSize, useShortName=T)
            
            # Create a summary table of the correlations
            datCorLong <- datCorLong %>% filter(!is.na(correlation)) %>% select(varA, varB, correlation, p, corrForCol)
            thisInd <- thisInd %>% select(shortName, indicator, region, subregion, trophicGroup)
            datCorLong <- left_join(datCorLong, thisInd, by=c("varA"="shortName")) %>% rename(indicatorA=indicator, trophicGroupA=trophicGroup,
                                                                                                  regionA=region, subregionA=subregion)
            datCorLong <- left_join(datCorLong, thisInd, by=c("varB"="shortName")) %>% rename(indicatorB=indicator, trophicGroupB=trophicGroup,
                                                                                                  regionB=region, subregionB=subregion)
            
            datCorLong <- datCorLong %>% select(indicatorA, indicatorB, correlation, p, varA, trophicGroupA, regionA, subregionA, 
                                                varB, trophicGroupB, regionB, subregionB, corrForCol) %>% arrange(p) %>% 
                mutate(correlation=round(correlation, 4), p=round(p, 4))
            
            # Calculate colors
            datCorLong$col <- "#EEE9E9"
            thisIndices <- datCorLong$correlation<0 & !is.na(datCorLong$corrForCol)
            datCorLong[thisIndices, "col"] <- rgb(colorRamp(color=c("white", "firebrick2"))(-datCorLong[thisIndices, "corrForCol"][[1]]), maxColorValue=255)
            thisIndices <- datCorLong$correlation>0 & !is.na(datCorLong$corrForCol)
            datCorLong[thisIndices, "col"] <- rgb(colorRamp(color=c("white", "dodgerblue4"))(datCorLong[thisIndices, "corrForCol"][[1]]), maxColorValue=255)
            datCorLong$corrForCol <- NULL
            corrDT <- datatable(datCorLong, rownames=F) %>% formatStyle("col", target="row", color=styleValue())
            saveWidget(corrDT, file.path(outputDir, paste0("corr_", scen, "_", minYear, ".html")))
        }
    }
}
clusData <- bind_rows(clusDataList)
corrData <- bind_rows(corrDataList)

# Assign more meaningful names for SEM
clusData <- clusData %>% rename(dendID=scen, dendClusID=group)

saveRDS(clusData, file=file.path(workingDir, "output", "clusData.rds"))
saveRDS(corrData, file=file.path(workingDir, "output", "corrData.rds"))

cat("===================================================================\n")
cat("Cluster assignment\n")
# Create a summary data frame showing the cluster assignments
clusters <- clusData %>% select(SEMlatent, dendClusID, shortName) %>% distinct() %>% arrange(SEMlatent, dendClusID, shortName)
write.csv(clusters, file=file.path(workingDir, "output", "clusters.csv"), row.names=F)
print(clusters)

cat("===================================================================\n")
cat("Number of indicators per cluster\n")
# Calculate the number of indicators per cluster
numPerCluster <- clusters %>% group_by(SEMlatent, dendClusID) %>% summarize(numPerCluster=n(), .groups="drop")
write.csv(numPerCluster, file=file.path(workingDir, "output", "numPerCluster.csv"), row.names=F)
print(numPerCluster, n=100)

cat("===================================================================\n")
cat("Number of clusters\n")
# Calculate the number of clusters
numClusters <- numPerCluster %>% group_by(SEMlatent) %>% summarize(numClusters=n(), .groups="drop")
write.csv(numClusters, file=file.path(workingDir, "output", "numClusters.csv"), row.names=F)
print(numClusters)

cat("===================================================================\n")
cat("Maximum number of indicators per cluster\n")
# Calculate the maximum number of indicators the calculated clusters would imply:
maxNumPerCluster <- numPerCluster %>% group_by(SEMlatent) %>% summarize(maxNumPerCluster=max(numPerCluster), .groups="drop")
write.csv(maxNumPerCluster, file=file.path(workingDir, "output", "maxNumPerCluster.csv"), row.names=F)
print(maxNumPerCluster)

cat("===================================================================\n")
cat("Maximum number of indicators, north: ", 
    sum(maxNumPerCluster$maxNumPerCluster[maxNumPerCluster$SEMlatent!="predator"]), "\n")
cat("Maximum number of indicators, predator: ",
    sum(maxNumPerCluster$maxNumPerCluster[maxNumPerCluster$SEMlatent!="north"]), "\n")
