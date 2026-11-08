# Run Structural Equation Model (SEM) using various combinations of indicators
#
# Doug Jackson
# doug@QEDAconsulting.com
library(tidyverse)
library(lubridate)
library(lavaan)
library(performance)
library(evaluate)
library(igraph)
library(readxl)
library(doParallel)
###########################################################################
# Constants
###########################################################################
workingDir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices"

# Select the model structure
modStruct <- c("pred", "north", "guildsA", "DAG1A", "DAG1B", "DAG1C")[5]

# Latent variables to include
if(modStruct=="pred") {
    latentVars <- c("prey", "cond1", "cond2", "predator")
} else if(modStruct=="prey") {
    latentVars <- c("prey", "cond1", "cond2", "north")
} else if(modStruct=="guildsA") {
    latentVars <- c("PreyNCC", "Growth", "Abundance", "PreyAK")
} else if(modStruct=="DAG1A" | modStruct=="DAG1B" | modStruct=="DAG1C") {
    latentVars <- c("PreyNCC", "PredNCC", "Growth", "Abundance", "PreyAK", "PredAK")
}

# Indicate whether to exclude the north/pred -> cond2 link
excludeLinkToCond2 <- FALSE

# Minimum year at which to truncate time series (must be one of the options used in analyzeIndicators.R)
useMinYear <- 1998

# analysis type:
# short = analyze years when all time series have complete data
# long = specify range of years and exclude incomplete time series
analysisType <- c("short", "long")[1]
longStartYear <- 1998
longEndYear <- 2021

useDFA <- TRUE
clusDataFile <- ifelse(useDFA, file.path(workingDir, "output", "DFA", "clusDataDFA.rds"),
                       file.path(workingDir, "output", "clusData.rds"))

# This is only used if useDFA==FALSE
SARfile <- file.path(workingDir, "data", "sar.19982021.csv")

# Fix indicator and latent variable variances to one
fixVariances <- FALSE

# Specify number of top models to plot and metric used for ranking
# If rankMetric==NA, ranking will be by the highest ranking across all metrics
rankMetric <- "AIC"
numTopModels <- 4

numCores <- 30
###########################################################################
# Functions
###########################################################################
# Obtain data for a specified latent variable and cluster ID
obtainClus <- function(latentVar, clusID) {
    thisClus <- clusters %>% filter(SEMlatent==latentVar & latentClusID==clusID)
    dat <- clusData %>% 
        filter(SEMlatent==latentVar & dendID==thisClus$dendID & dendClusID==thisClus$dendClusID) %>% 
        dplyr::select(shortName, date, finalVal) %>% 
        pivot_wider(id_cols=date, names_from=shortName, values_from=finalVal)
    
    # Standardize
    thisDate <- dat$date
    dat <- dat %>% dplyr::select(-date) %>% mutate_all(~(scale(.) %>% as.vector)) %>% 
        mutate(date=thisDate) %>% dplyr::select(date, everything())
    
    indNames <- names(dat)
    indNames <- indNames[indNames!="date"]
    indNames <- paste(indNames, collapse=", ")
    
    return(list(dat=dat, indNames=indNames))
}

createModStr <- function(latentVar, dat) {
    indNames <- names(dat)
    indNames <- indNames[indNames!="date"]
    
    if(length(indNames)>1) {
        modStr <- paste(latentVar, "=~", paste(indNames, collapse=" + "))
    } else {
        modStr <- ""
        # Change the column name to equal the latent variable name
        names(dat) <- c("date", latentVar)
    }
    
    if(fixVariances) {
        if(length(indNames)>1) {
            modStr <- paste0(modStr, paste0("\n", indNames, "~~1*", indNames, collapse=" " ))
        } else {
            modStr <- paste0(modStr, "\n", latentVar, "~~1*", latentVar)
        }
    }
    
    return(list(modStr=modStr, dat=dat))
}

constructSEMnorth <- function(thisPrey, thisCond1, thisCond2, thisNorth, logSAR,
                              preyModStr, cond1modStr, cond2modStr, northModStr,
                              preyIndNames, cond1indNames, cond2indNames, northIndNames) {
    thisDat <- full_join(thisPrey, thisCond1, by="date")
    thisDat <- full_join(thisDat, thisCond2, by="date")
    thisDat <- full_join(thisDat, thisNorth, by="date")
    thisDat <- full_join(thisDat, logSAR, by="date")
    
    # Filter years with missing logSAR
    thisDat <- thisDat %>% filter(!is.na(logSAR)) %>% arrange(date)
    
    # Standardize indicators
    origDate <- thisDat$date
    origLogSAR <- thisDat$logSAR
    thisDat <- thisDat %>% mutate_all(~(scale(.)) %>% as.vector)
    thisDat$date <- origDate
    thisDat$logSAR <- origLogSAR
    
    if(excludeLinkToCond2) {
        modStr <- paste(c(cond1modStr, cond2modStr, preyModStr, northModStr, "cond1 ~ prey", "cond2 ~ cond1", "logSAR ~ cond2 + north"),
                        collapse="\n")
    } else {
        modStr <- paste(c(cond1modStr, cond2modStr, preyModStr, northModStr, "cond1 ~ prey", "cond2 ~ north + cond1", "logSAR ~ cond2 + north"),
                        collapse="\n")
    }
    
    if(fixVariances) {
        modStr <- paste0(modStr, "\nlogSAR~~1*logSAR")
    }
    
    return(list(thisDat=thisDat, modStr=modStr))
}

constructSEMpred <- function(thisPrey, thisCond1, thisCond2, thisPred, logSAR,
                              preyModStr, cond1modStr, cond2modStr, predModStr) {
    thisDat <- full_join(thisPrey, thisCond1, by="date")
    thisDat <- full_join(thisDat, thisCond2, by="date")
    thisDat <- full_join(thisDat, thisPred, by="date")
    thisDat <- full_join(thisDat, logSAR, by="date")
    
    # Filter years with missing logSAR
    thisDat <- thisDat %>% filter(!is.na(logSAR)) %>% arrange(date)
    
    # Standardize indicators
    origDate <- thisDat$date
    origLogSAR <- thisDat$logSAR
    thisDat <- thisDat %>% mutate_all(~(scale(.)) %>% as.vector)
    thisDat$date <- origDate
    thisDat$logSAR <- origLogSAR
    
    if(excludeLinkToCond2) {
        modStr <- paste(c(cond1modStr, cond2modStr, preyModStr, predModStr, "cond1 ~ prey", "cond2 ~ cond1", "logSAR ~ cond2 + predator"),
                        collapse="\n")
    } else {
        modStr <- paste(c(cond1modStr, cond2modStr, preyModStr, predModStr, "cond1 ~ prey", "cond2 ~ predator + cond1", "logSAR ~ cond2 + predator"),
                        collapse="\n")
    }
    
    if(fixVariances) {
        modStr <- paste0(modStr, "\nlogSAR~~1*logSAR")
    }
    
    return(list(thisDat=thisDat, modStr=modStr))
}

constructSEMguildsA <- function(thisPreyNCC, thisGrowth, thisAbundance, thisPreyAK, 
                                logSAR,
                                preyNCCmodStr, growthModStr, abundanceModStr, preyAKmodStr) {
    thisDat <- full_join(thisPreyNCC, thisGrowth, by="date")
    thisDat <- full_join(thisDat, thisAbundance, by="date")
    thisDat <- full_join(thisDat, thisPreyAK, by="date")
    thisDat <- full_join(thisDat, logSAR, by="date")
    
    # Filter years with missing logSAR
    thisDat <- thisDat %>% filter(!is.na(logSAR)) %>% arrange(date)
    
    # Standardize indicators
    origDate <- thisDat$date
    origLogSAR <- thisDat$logSAR
    thisDat <- thisDat %>% mutate_all(~(scale(.)) %>% as.vector)
    thisDat$date <- origDate
    thisDat$logSAR <- origLogSAR
    
    modStr <- paste(c(preyNCCmodStr, growthModStr, abundanceModStr, preyAKmodStr,
                      "Growth ~ PreyNCC", "Abundance ~ Growth", "logSAR ~ Abundance + PreyAK"),
                    collapse="\n")
    
    if(fixVariances) {
        modStr <- paste0(modStr, "\nlogSAR~~1*logSAR")
    }
    
    return(list(thisDat=thisDat, modStr=modStr))
}

constructSEM_DAG1A <- function(thisPreyNCC, thisPredNCC, thisGrowth, thisAbundance, thisPreyAK, thisPredAK, 
                               logSAR,
                               preyNCCmodStr, predNCCmodStr, growthModStr, abundanceModStr, preyAKmodStr, predAKmodStr) {
    thisDat <- full_join(thisPreyNCC, thisPredNCC, by="date")
    thisDat <- full_join(thisDat, thisGrowth, by="date")
    thisDat <- full_join(thisDat, thisAbundance, by="date")
    thisDat <- full_join(thisDat, thisPreyAK, by="date")
    thisDat <- full_join(thisDat, thisPredAK, by="date")
    thisDat <- full_join(thisDat, logSAR, by="date")
    
    # Filter years with missing logSAR
    thisDat <- thisDat %>% filter(!is.na(logSAR)) %>% arrange(date)
    
    # Standardize indicators
    origDate <- thisDat$date
    origLogSAR <- thisDat$logSAR
    thisDat <- thisDat %>% mutate_all(~(scale(.)) %>% as.vector)
    thisDat$date <- origDate
    thisDat$logSAR <- origLogSAR

    modStr <- paste(c(preyNCCmodStr, predNCCmodStr, growthModStr, abundanceModStr, preyAKmodStr, predAKmodStr,
                      "Growth ~ PreyNCC", "Abundance ~ Growth + PredNCC", "logSAR ~ Abundance + PreyAK + PredAK"),
                    collapse="\n")
    
    if(fixVariances) {
        modStr <- paste0(modStr, "\nlogSAR~~1*logSAR")
    }
    
    return(list(thisDat=thisDat, modStr=modStr))
}

constructSEM_DAG1B <- function(thisPreyNCC, thisPredNCC, thisGrowth, thisAbundance, thisPreyAK, thisPredAK, 
                                logSAR,
                                preyNCCmodStr, predNCCmodStr, growthModStr, abundanceModStr, preyAKmodStr, predAKmodStr) {
    thisDat <- full_join(thisPreyNCC, thisPredNCC, by="date")
    thisDat <- full_join(thisDat, thisGrowth, by="date")
    thisDat <- full_join(thisDat, thisAbundance, by="date")
    thisDat <- full_join(thisDat, thisPreyAK, by="date")
    thisDat <- full_join(thisDat, thisPredAK, by="date")
    thisDat <- full_join(thisDat, logSAR, by="date")
    
    # Filter years with missing logSAR
    thisDat <- thisDat %>% filter(!is.na(logSAR)) %>% arrange(date)
    
    # Standardize indicators
    origDate <- thisDat$date
    origLogSAR <- thisDat$logSAR
    thisDat <- thisDat %>% mutate_all(~(scale(.)) %>% as.vector)
    thisDat$date <- origDate
    thisDat$logSAR <- origLogSAR
    
    modStr <- paste(c(preyNCCmodStr, predNCCmodStr, growthModStr, abundanceModStr, preyAKmodStr, predAKmodStr,
                      "Growth ~ PreyNCC", "Abundance ~ PredNCC", "logSAR ~ Growth + Abundance + PreyAK + PredAK"),
                    collapse="\n")
    
    if(fixVariances) {
        modStr <- paste0(modStr, "\nlogSAR~~1*logSAR")
    }
    
    return(list(thisDat=thisDat, modStr=modStr))
}

constructSEM_DAG1C <- function(thisPreyNCC, thisPredNCC, thisGrowth, thisAbundance, thisPreyAK, thisPredAK, 
                               logSAR,
                               preyNCCmodStr, predNCCmodStr, growthModStr, abundanceModStr, preyAKmodStr, predAKmodStr) {
    thisDat <- full_join(thisPreyNCC, thisPredNCC, by="date")
    thisDat <- full_join(thisDat, thisGrowth, by="date")
    thisDat <- full_join(thisDat, thisAbundance, by="date")
    thisDat <- full_join(thisDat, thisPreyAK, by="date")
    thisDat <- full_join(thisDat, thisPredAK, by="date")
    thisDat <- full_join(thisDat, logSAR, by="date")
    
    # Filter years with missing logSAR
    thisDat <- thisDat %>% filter(!is.na(logSAR)) %>% arrange(date)
    
    # Standardize indicators
    origDate <- thisDat$date
    origLogSAR <- thisDat$logSAR
    thisDat <- thisDat %>% mutate_all(~(scale(.)) %>% as.vector)
    thisDat$date <- origDate
    thisDat$logSAR <- origLogSAR
    
    modStr <- paste(c(preyNCCmodStr, predNCCmodStr, growthModStr, abundanceModStr, preyAKmodStr, predAKmodStr,
                      "Abundance ~ PreyNCC + PredNCC", "logSAR ~ PreyNCC + Growth + Abundance + PredNCC + PreyAK + PredAK"),
                    collapse="\n")
    
    if(fixVariances) {
        modStr <- paste0(modStr, "\nlogSAR~~1*logSAR")
    }
    
    return(list(thisDat=thisDat, modStr=modStr))
}

evaluateSEM <- function(thisDat, modStr, thisIndNames) {
    # Fit the model
    rm(list=c("mod", "modPerf"))
    
    if(fixVariances) {
        fitOut <- evaluate::evaluate("mod <- sem(modStr, thisDat, std.lv=T)")
    } else {
        fitOut <- evaluate::evaluate("mod <- sem(modStr, thisDat)")
    }

    fitMessage <- ifelse(length(fitOut)>1, paste(unlist(fitOut[2:length(fitOut)]), collapse="; "), "")
    fitMessage <- gsub("\n", "", fitMessage)
    fitMessage <- gsub("\\s+", " ", fitMessage)
    
    if(exists("mod")) {
        perfOut <- evaluate::evaluate(paste0("modPerf <- model_performance(mod, metrics = c('", paste(metrics, collapse="', '"), "'))"))
    } else {
        mod <- NA
    }
    
    # Save the results
    thisResults <- thisComb
    if(exists("modPerf")) {
        for(metric in metrics) {
            thisResults[, metric] <- modPerf[[metric]]
        }        
    } else {
        for(metric in metrics) {
            thisResults[, metric] <- NA
        }     
    }
    
    if(modStruct=="pred") {
        thisResults$preyIndNames <- thisIndNames$prey
        thisResults$cond1indNames <- thisIndNames$cond1
        thisResults$cond2indNames <- thisIndNames$cond2
        thisResults$predIndNames <- thisIndNames$pred
    } else if(modStruct=="north") {
        thisResults$preyIndNames <- thisIndNames$prey
        thisResults$cond1indNames <- thisIndNames$cond1
        thisResults$cond2indNames <- thisIndNames$cond2
        thisResults$northIndNames <- thisIndNames$north
    } else if(modStruct=="guildsA") {
        thisResults$PreyNCCindNames <- thisIndNames$PreyNCC
        thisResults$GrowthIndNames <- thisIndNames$Growth
        thisResults$AbundanceIndNames <- thisIndNames$Abundance
        thisResults$PreyAKindNames <- thisIndNames$PreyAK
    } else if(modStruct=="DAG1A" | modStruct=="DAG1B" | modStruct=="DAG1C") {
        thisResults$PreyNCCindNames <- thisIndNames$PreyNCC
        thisResults$PredNCCindNames <- thisIndNames$PredNCC
        thisResults$GrowthIndNames <- thisIndNames$Growth
        thisResults$AbundanceIndNames <- thisIndNames$Abundance
        thisResults$PreyAKindNames <- thisIndNames$PreyAK
        thisResults$PredAKindNames <- thisIndNames$PredAK
    }
    thisResults$messages <- fitMessage
    
    return(list(thisResults=thisResults, mod=mod))
}

createPathFig <- function(mod, modNum) {
    
    # Path figure code from Noble, Eric, et al.
    parameters <- parameterEstimates(mod, standardized=T)
    
    # Prepare edges data for igraph
    edges <- parameters %>%
        filter(op %in% c('=~', '~')) %>%
        dplyr::select(lhs, rhs, estimate = est, se = se)
    
    g <- graph_from_data_frame(edges, directed = TRUE)# Create an igraph object
    layout <- layout_nicely(g)# Automatically layout the graph
    
    # Extract node positions and merge with node names
    node_positions <- data.frame(name=V(g)$name, x=layout[, 1], y=layout[, 2])
    
    # Merge node positions with edges to create complete plot data
    plot_data <- merge(edges, node_positions, by.x = "lhs", by.y = "name", all.x = TRUE)
    colnames(plot_data)[(ncol(plot_data)-1):ncol(plot_data)] <- c("lhs_x", "lhs_y")
    
    plot_data <- merge(plot_data, node_positions, by.x = "rhs", by.y = "name", all.x = TRUE)
    colnames(plot_data)[(ncol(plot_data)-1):ncol(plot_data)] <- c("rhs_x", "rhs_y")
    
    # Modify plot_data to include a new column for adjusted size
    plot_data$size <- 1 / plot_data$se
    plot_data$size[is.infinite(plot_data$size)] <- max(plot_data$size[!is.infinite(plot_data$size)], na.rm = TRUE)  # Replace Inf with max size
    plot_data$size[is.na(plot_data$size) | plot_data$size <= 0] <- min(plot_data$size[plot_data$size > 0], na.rm = TRUE) / 2  # Replace NA or non-positive with smallest positive size
    
    p <- ggplot(plot_data, aes(x = lhs_x, y = lhs_y, xend = rhs_x, yend = rhs_y)) +
        geom_segment(aes(size = size, alpha = size),  # Use adjusted size for both size and alpha
                     arrow = arrow(length = unit(0.02, "npc")), color = "grey") +
        geom_text(aes(x = (lhs_x + rhs_x) / 2, y = (lhs_y + rhs_y) / 2, label = round(estimate, 2)), vjust = 1, color = "blue") +
        geom_point(aes(x = lhs_x, y = lhs_y), color = "red", size = 3) +
        geom_point(aes(x = rhs_x, y = rhs_y), color = "red", size = 3) +
        geom_text(aes(x = lhs_x, y = lhs_y, label = lhs), vjust = -1, color = "red", size = 0.5) +
        geom_text(aes(x = rhs_x, y = rhs_y, label = rhs), vjust = -1, color = "red") +
        scale_size_continuous(range = c(0.5, 3)) +  # Control size range for better visualization (0.5,3)
        scale_alpha_continuous(range = c(0.5, 1)) +  # Control transparency range
        theme_minimal() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.position = "right") +
        labs(title=paste("model", modNum)) +
        scale_x_continuous(expand = expansion(mult = 0.1)) + scale_y_continuous(expand = expansion(mult = 0.1)) +
        backgroundTheme
    
    ggsave(file.path(figOutputDir, paste0("pathFig_mod_", modNum, ".png")))
    
}

###########################################################################
# Run
###########################################################################
setwd(workingDir)

outputDir <- file.path(workingDir, "output", "SEM")
dir.create(outputDir, showWarnings=F)

figOutputDir <- file.path(outputDir, "figs")
dir.create(figOutputDir, showWarnings=F)

doParOutputDir <- file.path(outputDir, "doPar")
success <- unlink(doParOutputDir, recursive=T)
if(success==1) {
    cat("Failed to delete doPar output directory,", doParOutputDir, "\n")
}
dir.create(doParOutputDir, showWarnings=F)

# Set up parallel processing
registerDoParallel(cores=numCores)

# Load clusData that was generated by analyzeIndicators.R
clusData <- readRDS(clusDataFile)

# Load smolt-to-adult return ratios
if(useDFA) {
    # SAR in DFA was already log-transformed
    logSAR <- clusData %>% filter(SEMlatent=="SAR") %>% select(date, finalVal) %>% rename(logSAR=finalVal) %>% 
        dplyr::select(date, logSAR)
} else {
    SAR <- read.csv(SARfile)
    SAR$date <- ymd(SAR$Year, truncated=2L)
    SAR <- SAR[, c("date", "SAR")]
    logSAR <- SAR %>% mutate(logSAR=log(SAR)) %>% dplyr::select(date, logSAR)
}

clusData <- clusData %>% filter(minYear==useMinYear, SEMlatent %in% latentVars)

# Create completeness data frame
completeness <- clusData %>% select(shortName, date, finalVal) %>% 
    pivot_wider(id_cols=date, values_from=finalVal, names_from=shortName) %>% arrange(date)
completeness$complete <- complete.cases(completeness)
write.csv(completeness, file.path(outputDir, "completeness.csv"), row.names=F)

if(analysisType=="short") {
    completeDates <- completeness[completeness$complete, "date"]
    clusData <- clusData %>% filter(date %in% completeDates$date)
} else if(analysisType=="long") {
    completeness$year <- year(ymd(completeness$date))
    completeness <- completeness %>% filter(year>=longStartYear, year<=longEndYear) %>% arrange(year)
    complete <- sapply(completeness, function(col) all(complete.cases(col)))
    completeInds <- names(complete[complete])
    clusData <- clusData %>% mutate(year=year(ymd(date))) %>% 
                                        filter(year>=longStartYear, year<=longEndYear, shortName %in% completeInds)
}

# Create completeness data frame (post-filter)
completenessFiltered <- clusData %>% select(shortName, date, finalVal) %>% 
    pivot_wider(id_cols=date, values_from=finalVal, names_from=shortName)
completenessFiltered$complete <- complete.cases(completenessFiltered)
write.csv(completenessFiltered, file.path(outputDir, paste0("completeness_", analysisType, ".csv")), row.names=F)

metrics <- c('AIC', 'p_Chi2', 'AGFI', 'CFI', 'SRMR', 'RMSEA')

backgroundTheme <- theme(panel.background=element_rect(fill="white"),
                         plot.background=element_rect(fill="white"))

###########################################################################
# Run different combinations of clusters

# Identify unique clusters
clusters <- clusData %>% dplyr::select(SEMlatent, dendID, dendClusID) %>% distinct() %>% 
    arrange(SEMlatent, dendID, dendClusID) %>% group_by(SEMlatent) %>% mutate(latentClusID=row_number())

# Calculate total number of clusters per latent variable
numClus <- clusters %>% group_by(SEMlatent) %>% summarize(numClus=n(), .groups="drop")

# Generate all possible combinations of clusters
sequences <- lapply(numClus$numClus, function(x) seq_len(x))
names(sequences) <- numClus$SEMlatent
combinations <- expand.grid(sequences)

logFile <- file.path(outputDir, "runSEM.log")
cat("\nTrack progress in logfile:", logFile, "\n")
cat("Number of cores to use (numCores):", numCores, "\n")
cat("Check progress in", doParOutputDir, "\n")
cat("runSEM logfile\n", file=logFile)

# Run SEM on all combinations
numCombinations <- nrow(combinations)
startTime <- Sys.time()
r <- foreach(i=1:numCombinations, .packages=c("tidyverse", "lubridate", "lavaan", "dplyr", "performance")) %dopar% { 
    cat("Running SEM model", i, "of", nrow(combinations), "\n")
    cat("Running SEM model", i, "of", nrow(combinations), "\n", file=logFile, append=T)
    thisComb <- combinations[i, ]
    
    thisDat <- list()
    thisIndNames <- list()
    thisModStr <- list()
    numInd <- 0
    
    # Assemble the data. If there's only one indicator, the column name
    # will be set to the latent variable name.
    for(latentVar in latentVars) {
        out <- obtainClus(latentVar, thisComb[, latentVar])
        thisDat[[latentVar]] <- out$dat
        thisIndNames[[latentVar]] <- out$indNames
        numInd <- numInd + length(str_split(thisIndNames[[latentVar]], ",")[[1]])
        out <- createModStr(latentVar, thisDat[[latentVar]])
        thisModStr[[latentVar]] <- out$modStr
        thisDat[[latentVar]] <- out$dat
    }

    if(modStruct=="north") {
        thisMod <- constructSEMnorth(thisDat$prey, thisDat$cond1, thisDat$cond2, thisDat$north, logSAR,
                                     thisModStr$prey, thisModStr$cond1, thisModStr$cond2, thisModStr$north)
        thisResults <- evaluateSEM(thisMod$thisDat, thisMod$modStr, thisIndNames)
    } else if(modStruct=="pred") {
        thisMod <- constructSEMpred(thisDat$prey, thisDat$cond1, thisDat$cond2, thisDat$pred, logSAR,
                                    thisModStr$prey, thisModStr$cond1, thisModStr$cond2, thisModStr$pred)
        thisResults <- evaluateSEM(thisMod$thisDat, thisMod$modStr, thisIndNames)
    } else if(modStruct=="guildsA") {
        thisMod <- constructSEMguildsA(thisDat$PreyNCC, thisDat$Growth, thisDat$Abundance, thisDat$PreyAK,
                                      logSAR,
                                      thisModStr$PreyNCC, thisModStr$Growth, thisModStr$Abundance, thisModStr$PreyAK)
        thisResults <- evaluateSEM(thisMod$thisDat, thisMod$modStr, thisIndNames)        
    } else if(modStruct=="DAG1A") {
        thisMod <- constructSEM_DAG1A(thisDat$PreyNCC, thisDat$PredNCC, thisDat$Growth, thisDat$Abundance, thisDat$PreyAK, thisDat$PredAK, 
                    logSAR,
                    thisModStr$PreyNCC, thisModStr$PredNCC, thisModStr$Growth, thisModStr$Abundance, thisModStr$PreyAK, thisModStr$PredAK)
        thisResults <- evaluateSEM(thisMod$thisDat, thisMod$modStr, thisIndNames)
    } else if(modStruct=="DAG1B") {
        thisMod <- constructSEM_DAG1B(thisDat$PreyNCC, thisDat$PredNCC, thisDat$Growth, thisDat$Abundance, thisDat$PreyAK, thisDat$PredAK, 
                                       logSAR,
                                       thisModStr$PreyNCC, thisModStr$PredNCC, thisModStr$Growth, thisModStr$Abundance, thisModStr$PreyAK, thisModStr$PredAK)
        thisResults <- evaluateSEM(thisMod$thisDat, thisMod$modStr, thisIndNames)
    } else if(modStruct=="DAG1C") {
        thisMod <- constructSEM_DAG1C(thisDat$PreyNCC, thisDat$PredNCC, thisDat$Growth, thisDat$Abundance, thisDat$PreyAK, thisDat$PredAK, 
                                      logSAR,
                                      thisModStr$PreyNCC, thisModStr$PredNCC, thisModStr$Growth, thisModStr$Abundance, thisModStr$PreyAK, thisModStr$PredAK)
        thisResults <- evaluateSEM(thisMod$thisDat, thisMod$modStr, thisIndNames)
    }
    
    thisResults$thisResults$numObs <- ifelse(is.na(thisResults$mod), NA, thisResults$mod@SampleStats@nobs[[1]])
    thisResults$thisResults$numInd <- numInd
    thisResults$thisResults$modNum <- i
    thisResults$thisResults$numIterations <- thisResults$mod@optim$iterations
    
    # Save results to output files
    saveRDS(thisResults$thisResults, file=file.path(doParOutputDir, paste0("resultsList_", i, ".rds")))
    saveRDS(thisResults$mod, file=file.path(doParOutputDir, paste0("modList_", i, ".rds")))
    if(!is.na(thisResults$mod)) {
        thisParamEst <- parameterEstimates(thisResults$mod)
        thisParamEst$modNum <- i
        saveRDS(thisParamEst, file=file.path(doParOutputDir, paste0("paramEstList_", i, ".rds")))
    }

}

# Read the doPar outputs
cat("Reading doParallel outputs.\n")
resultsList <- list()
modList <- list()
paramEstList <- list()
for(i in 1:numCombinations) {
    cat("Loading results for model", i, "of", numCombinations, "\n", file=logFile, append=T)
    thisResultsList <- readRDS(file.path(doParOutputDir, paste0("resultsList_", i, ".rds")))
    thisModList <- readRDS(file.path(doParOutputDir, paste0("modList_", i, ".rds")))
    
    resultsList[[length(resultsList)+1]] <- thisResultsList
    modList[[length(modList)+1]] <- thisModList
    
    if(!is.na(thisModList)) {
        thisParamEst <- parameterEstimates(thisModList)
        thisParamEst$modNum <- i
        paramEstList[[length(paramEstList)+1]] <- thisParamEst
    }
}

# # Save list of models to an RDS
# cat("Saving list of models to RDS.\n")
# saveRDS(modList, file=file.path(outputDir, "modList.rds"))

resultsByClus <- bind_rows(resultsList)

if(modStruct=="north") {
    resultsByClus <- resultsByClus %>% dplyr::rename(cond1cluster=cond1, cond2cluster=cond2, preyCluster=prey, northCluster=north) %>% 
        relocate(preyCluster) %>% relocate(modNum)
} else if(modStruct=="pred") {
    resultsByClus <- resultsByClus %>% dplyr::rename(cond1cluster=cond1, cond2cluster=cond2, preyCluster=prey, predCluster=predator) %>% 
        relocate(preyCluster) %>% relocate(modNum)
} else if(modStruct=="guildsA") {
    resultsByClus <- resultsByClus %>% dplyr::rename(PreyNCCclus=PreyNCC, GrowthClus=Growth, AbundanceClus=Abundance, PreyAKclus=PreyAK) %>% 
        relocate(PreyAKclus) %>% relocate(AbundanceClus) %>% relocate(GrowthClus) %>% relocate(PreyNCCclus) %>% relocate(modNum)    
} else if(modStruct=="DAG1A" | modStruct=="DAG1B" | modStruct=="DAG1C") {
    resultsByClus <- resultsByClus %>% dplyr::rename(PreyNCCclus=PreyNCC, PredNCCclus=PredNCC, GrowthClus=Growth, AbundanceClus=Abundance,
                                                     PreyAKclus=PreyAK, PredAKclus=PredAK) %>% 
        relocate(PreyAKclus) %>% relocate(AbundanceClus) %>% relocate(GrowthClus) %>% relocate(PredNCCclus) %>% relocate(PreyNCCclus) %>% relocate(modNum)
}

# Rank by each metric and then calculate the highest rank
resultsByClus <- resultsByClus %>% mutate(rankAIC=dense_rank(AIC), rankChi2=dense_rank(desc(p_Chi2)), 
                                          rankAGFI=dense_rank(desc(AGFI)), rankCFI=dense_rank(desc(CFI)),
                                          rankSRMR=dense_rank(SRMR), rankRMSEA=dense_rank(RMSEA)) %>% 
    mutate(highestRank=do.call(pmin, across(starts_with("rank")))) %>% arrange(highestRank)

if(!is.na(rankMetric)) {
    resultsByClus <- resultsByClus[order(resultsByClus[, "AIC"]), ]
}

write.csv(resultsByClus, file=file.path(outputDir, "SEMresultsByClus.csv"), row.names=F)

# Compile and save parameter estimates
paramEst <- bind_rows(paramEstList)
paramEst <- paramEst %>% relocate(modNum)
write.csv(paramEst, file=file.path(outputDir, "parameterEstimates.csv"), row.names=F)

cat("Runtime for by-cluster runs:\n")
print(Sys.time()-startTime)

###########################################################################
# Create path figures for top models
for(i in 1:numTopModels) {
    modNum <- resultsByClus$modNum[i]
    createPathFig(modList[[modNum]], modNum)   
}
###########################################################################
# Create ranked plots of metrics
rankedMetrics <- resultsByClus[, c("modNum", metrics)]
rankedMetrics <- rankedMetrics[complete.cases(rankedMetrics), ]
rankedMetrics <- rankedMetrics %>% pivot_longer(cols=-modNum, names_to="metric") %>% arrange(value) %>% 
    group_by(metric) %>% mutate(rank=row_number())
p <- ggplot(rankedMetrics) + geom_line(aes(x=rank, y=value)) + facet_wrap(~metric, ncol=1, scales="free_y") +
    theme_light() + backgroundTheme
ggsave(file.path(figOutputDir, "rankedMetrics.png"), width=6, height=6)

###########################################################################
# Plot AIC by indicator
for(latentVar in latentVars) {
    thisDat <- resultsByClus
    if(paste0(latentVar, "indNames") %in% names(thisDat)) {
        thisDat$ind <- thisDat[, paste0(latentVar, "indNames")]
    } else {
        thisDat$ind <- thisDat[, paste0(latentVar, "IndNames")]
    }
    
    meanAIC <- thisDat %>% group_by(ind) %>% summarize(meanAIC=mean(AIC), .groups="drop")
    p <- ggplot() + geom_histogram(data=thisDat, aes(x=AIC), bins=40) + facet_wrap(~ind, scales="free_y", ncol=1) +
        geom_vline(data=meanAIC, aes(xintercept=meanAIC, group=ind), col="red", size=1) +
        theme_light() + backgroundTheme
    ggsave(file=file.path(figOutputDir, paste0("histAIC_", latentVar, ".pdf")), width=8, height=10)
}

###########################################################################
# Calculate variable importance (see 20may25 Obsidian)
resultsByClus$deltaAIC <- resultsByClus$AIC - min(resultsByClus$AIC, na.rm=T)
resultsByClus$modelWeight <- exp(-0.5*resultsByClus$deltaAIC)/sum(exp(-0.5*resultsByClus$deltaAIC), na.rm=T)

allVars <- resultsByClus %>% select(ends_with("indNames", ignore.case=T)) %>% rowwise() %>%
    mutate(vars=paste(c_across(everything()), collapse=",")) %>%
    ungroup()
resultsByClus$allVars <- allVars$vars

varImportance <- data.frame(var=unique(clusData$shortName))
for(i in 1:nrow(varImportance)) {
    v <- varImportance$var[i]
    thisResults <- resultsByClus
    thisResults$inModel <- grepl(v, thisResults$allVars)
    varImportance[i, "importance"] <- sum(thisResults[thisResults$inModel, "modelWeight"])
}
varImportance <- varImportance %>% arrange(desc(importance))

# Add indicator information if using DFAs
if(useDFA) {
    guilds <- read.csv(file.path(workingDir, "output", "DFA", "guildsDFA.csv"))
    loadings <- read_excel(file.path(workingDir, "output", "DFA", "loadings.xlsx"))
    
    for(i in 1:nrow(varImportance)) {
        if(grepl("smoothed", varImportance[i, "var"])) {
            guildName <- gsub("_smoothed", "", varImportance[i, "var"])
            varImportance[i, "indicators"] <- guilds[guilds$guild==guildName, "shortName"]
            varImportance[i, "SEMnode"] <- guilds[guilds$guild==guildName, "guildSEMnode"]
        } else if(grepl("DFA1", varImportance[i, "var"])) {
            guildName <- gsub("_DFA1", "", varImportance[i, "var"])
            thisLoadings <- loadings[loadings$guild==guildName, ]
            thisLoadings <- thisLoadings %>% arrange(desc(abs(Z_est)))
            varImportance[i, "indicators"] <- paste0(thisLoadings$indicator, collapse=", ")
            varImportance[i, "SEMnode"] <- unique(thisLoadings$guildSEMnode)
        }
    }
}
write.csv(varImportance, file.path(outputDir, "varImportance.csv"), row.names=F)

# Make a copy of runSEM.R for future reference
file.copy(file.path(workingDir, "runSEM.R"), file.path(outputDir, "runSEM.R"))

###########################################################################
# Create version of completeness with indicators in column names
# NOTE: since the inds lookup relies on varImportance, only the "short" version of this 
# file will contain inds for all columns
completenessNames <- names(completeness)
completenessWithInds <- completeness
for(i in 1:length(completenessNames)) {
    
    thisName <- gsub("^X", "", completenessNames[i])
    
    thisInds <- varImportance$indicators[varImportance$var==thisName]
    
    names(completenessWithInds)[i] <- paste0(c(thisName, thisInds), collapse=": ")
}

write.csv(completenessWithInds, file.path(outputDir, "completenessWithInds.csv"), row.names=F)

cat("Total runtime:\n")
print(Sys.time()-startTime)