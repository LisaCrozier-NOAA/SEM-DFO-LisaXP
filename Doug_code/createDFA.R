# Create a new clusDataDFA data frame by applying DFA to each cluster in clusData.
# 
# Doug Jackson
# doug@QEDAconsulting.com
library(tidyverse)
library(MARSS)
library(lubridate)
#library(writexl)
# An alternative to writexl to prevent the following error:
# Error: package or namespace load failed for ‘writexl’ in inDL(x, as.logical(local), as.logical(now), ...):
#     unable to load shared object 'C:/Users/dougj/AppData/Local/R/win-library/4.5/writexl/libs/x64/writexl.dll':
#     LoadLibrary failure:  An Application Control policy has blocked this file.
library(openxlsx2)

# install.packages("devtools")
# devtools::install_github("thomasp85/patchwork")
library(patchwork)
###########################################################################
# Constants
###########################################################################
workingDir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices"

clusDataFile <- file.path(workingDir, "output", "clusData.rds")
allDataFile <- file.path(workingDir, "output", "allData.rds")
specifyDFAfile <- file.path(workingDir, "specifyDFA.csv")

guildFile <- file.path(workingDir, "guilds.csv")

# Column in guild file to use
guildCol <- c("shiftLisa_step2_26jun26")

# Number of dynamic factor analysis (DFA) factors
numFactors <- 1

# Minimum |loading| threshold
minLoading <- 0.2

clusMethod <- c("pamk", "guild")[2]
###########################################################################
# Run
###########################################################################
setwd(workingDir)

outputDir <- file.path(workingDir, "output", "DFA")
dir.create(outputDir, showWarnings=F)

clusData <- readRDS(clusDataFile)
allData <- readRDS(allDataFile)

clusData$year <- year(ymd(clusData$date))

backgroundTheme <- theme(panel.background=element_rect(fill="white"),
                         plot.background=element_rect(fill="white"))

###########################################################################
# Perform DFAs on clusters created by pamk
clusDataDFAlist <- list()
fitList <- list()
if(clusMethod=="pamk") {
    specifyDFA <- read.csv(specifyDFAfile)
    file.copy(specifyDFAfile, file.path(outputDir, basename(specifyDFAfile)))
    
    # Ensure that useDFA1 and useDFA2 are either Y or N
    for(col in c("useDFA1", "useDFA2")) {
        specifyDFA[is.na(specifyDFA[, col]), col] <- "N"
        specifyDFA[specifyDFA[, col]!="Y", col] <- "N"
    }
    
    clusters <- clusData %>% distinct(dendID, dendClusID)
    
    for(i in 1:nrow(clusters)) {
        
        thisCluster <- clusters[i, ]
        thisDendID <- thisCluster$dendID
        thisDendClusID <- thisCluster$dendClusID
        thisClusData <- clusData %>% filter(dendID==thisDendID, dendClusID==thisDendClusID)
        
        thisSpecifyDFA <- specifyDFA %>% filter(dendID==thisDendID, dendClusID==thisDendClusID)
        
        # Run DFA and create plots whether we end up using it in clusDataDFA or not
        if(length(unique(thisClusData$shortName))>1) {
            thisClusData <- thisClusData %>% select(shortName, date, finalVal) %>% 
                pivot_wider(id_cols=date, names_from=shortName, values_from=finalVal) %>% arrange(date)
            thisDate <- thisClusData$date
            thisClusData <- thisClusData %>% select(-date) %>% mutate_all(~(scale(.) %>% as.vector))
            
            # Following example from here: https://cran.r-project.org/web//packages/dsem/vignettes/dynamic_factor_analysis.html
            # Transpose so years are in columns and time series are in rows
            thisClusData <- t(thisClusData)
            fit <- MARSS(thisClusData, model=list(m=numFactors), form="dfa", method="BFGS")
            
            # Plot states using all data
            png(file.path(outputDir, paste0("factors_dendID_", thisCluster$dendID, "_cluster_", thisCluster$dendClusID, ".png")))
            plot(fit, plot.type="xtT")
            dev.off()
            
            # Plot expectation for data using all data
            png(file.path(outputDir, paste0("expectation_dendID_", thisCluster$dendID, "_cluster_", thisCluster$dendClusID, ".png")))
            plot(fit, plot.type="fitted.ytT")
            dev.off()
            
            # Following along with Noble's code (makeDFAs.R):
            Z_est <- coef(fit, type="matrix")$Z
            # Get the inverse of the rotation matrix
            H_inv <- varimax(Z_est)$rotmat
            # Rotate factor loadings
            Z_rot = Z_est %*% H_inv
            # Rotate processes
            processes = solve(H_inv) %*% fit$states
            
            # Save the loadings
            save(file=file.path(outputDir, paste0("loadings_dendID_", thisCluster$dendID, "_cluster_", thisCluster$dendClusID, ".Rdata")),
                 list=c("Z_est", "Z_rot"))
            
            # Construct the clusData entries
            thisClusData <- clusData %>% filter(dendID==thisCluster$dendID, dendClusID==thisCluster$dendClusID)
            thisSEMlatent <- unique(thisClusData$SEMlatent)
            thisMinYear <- unique(thisClusData$minYear)
            
            for(j in 1:numFactors) {
                # Don't add this DFA to clusData if it's not selected
                if(j==1 & thisSpecifyDFA$useDFA1=="N") {next}
                if(j==2 & thisSpecifyDFA$useDFA2=="N") {next}
                
                thisDFA <- data.frame(SEMlatent=thisSEMlatent, dendID=thisDendID, dendClusID=thisDendClusID,
                                      shortName=paste(thisDendID, thisDendClusID, paste0("DFA", j), sep="_"), date=thisDate, finalVal=processes[j, ],
                                      minYear=thisMinYear)
                clusDataDFAlist[[length(clusDataDFAlist)+1]] <- thisDFA
            }
        }
        
        # If neither DFA is marked for inclusion in clusDataDFA, include the original indicators
        if(length(unique(thisClusData$shortName))==1 | (thisSpecifyDFA$useDFA1=="N" & thisSpecifyDFA$useDFA2=="N")) {
            thisDFA <- thisClusData %>% select(SEMlatent, dendID, dendClusID, shortName, date, finalVal, minYear)
            clusDataDFAlist[[length(clusDataDFAlist)+1]] <- thisDFA
        }
    }
} else if (clusMethod=="guild") {
    guilds <- read.csv(guildFile)
    
    # Only include indicators that are in clusData
    guilds$guild <- guilds[, guildCol]
    guilds <- guilds %>% filter(shortName %in% clusData$shortName) %>% rename(guildSEMnode=SEMnode.for.AIC.analysis) %>% 
        select(shortName, guild, guildSEMnode)
    write.csv(guilds, file.path(outputDir, "guildsDFA.csv"), row.names=F)
    
    clusData <- left_join(clusData, guilds, by="shortName")
    
    # Create artificial dendID and dendClusID
    guilds <- guilds %>% distinct(guild, guildSEMnode) %>% filter(guild!="", guildSEMnode!="") %>% arrange(guildSEMnode, guild) %>% 
        group_by(guildSEMnode) %>% mutate(dendID=cur_group_id(), dendClusID=row_number()) %>% arrange(guild)
    
    plotList <- list()
    loadingsList <- list()
    smoothedIndicatorsList <- list()
    patchworkPlotList <- list()
    patchworkPlotIndex <- 0
    prevStraggler <- ""
    groupList <- list()
    for(guildRow in 1:nrow(guilds)) {
        
        thisGuildRow <- guilds[guildRow, ]
        thisGuild <- thisGuildRow$guild
        thisGuildSEMnode <- thisGuildRow$guildSEMnode
        thisDendID <- thisGuildRow$dendID
        thisDendClusID <- thisGuildRow$dendClusID
        
        thisClusData <- clusData %>% filter(guild==thisGuild, guildSEMnode==thisGuildSEMnode)
        thisMinYear <- unique(thisClusData$minYear)
        
        # Infer indicator ranking from order in guilds file
        ranks <- read.csv(guildFile)
        ranks <- ranks %>% filter(!!sym(guildCol)==thisGuild, SEMnode.for.AIC.analysis==thisGuildSEMnode) %>% 
            select(shortName) %>% mutate(rank=row_number())
        
        # Arrange by rank
        thisClusData <- left_join(thisClusData, ranks, by="shortName")
        thisClusData <- thisClusData %>% arrange(rank, year)
        
        # Run DFA and create plots whether we end up using it in clusDataDFA or not
        if(length(unique(thisClusData$shortName))>1) {
            thisGroup <- data.frame(guildSEMnode=thisGuildRow$guildSEMnode, guild=thisGuild, shortName=unique(thisClusData$shortName),
                                    type="DFA")
            groupList[[length(groupList)+1]] <- thisGroup
            
            thisClusData <- thisClusData %>% dplyr::select(shortName, date, finalVal) %>% 
                pivot_wider(id_cols=date, names_from=shortName, values_from=finalVal) %>% arrange(date)
            allMissing <- thisClusData %>% select(-date)
            allMissing$allMissing <- rowSums(is.na(allMissing))==ncol(allMissing)
            thisClusData <- thisClusData[!allMissing$allMissing, ]
            thisDate <- thisClusData$date
            
            thisClusData <- thisClusData %>% select(-date) %>% mutate_all(~(scale(.) %>% as.vector))
            
            # Following example from here: https://cran.r-project.org/web//packages/dsem/vignettes/dynamic_factor_analysis.html
            # Transpose so years are in columns and time series are in rows
            thisClusData <- t(thisClusData)
            fit <- MARSS(thisClusData, model=list(m=numFactors), form="dfa", method="BFGS")
            
            # Ensure that the sign of the DFA matches the indicator with the highest loading
            factor <- 1
            Z_matrix <- fit$par$Z
            maxIndex <- which.max(abs(Z_matrix[, factor]))
            
            if (Z_matrix[maxIndex, factor]<0) {
                fit$par$Z[, factor] <- -fit$par$Z[, factor]
                fit$states[factor, ] <- -fit$states[factor, ]
            }
            
            # Custom states plot
            ap <- autoplot(fit, plot.type="xtT")
            plotDF <- ap$data
            plotDF$datetime <- ymd(thisDate)
            p1 <- ggplot(plotDF, aes(x=datetime, y=.estimate)) + geom_ribbon(aes(ymin=.conf.low, ymax=.conf.up), fill="lightgray") + 
                geom_line(color="black", linewidth=1) +
                labs(title=paste0("factors, guildSEMnode = ", thisGuildSEMnode, ", guild = ", thisGuild), x="", y="Estimate") +
                theme_light()
            plotList[[length(plotList)+1]] <- p1
            
            # Plot expectation for data using all data
            ap <- autoplot(fit, plot.type="fitted.ytT")
            plotDF <- ap$data %>% rename(indName=.rownames, origDat=y, est=.fitted, confLow=.conf.low, confUp=.conf.up) %>% 
                select(indName, origDat, est, confLow, confUp)
            inds <- thisGroup$shortName
            numInds <- length(inds)
            plotDF$indName <- factor(plotDF$indName, levels=inds[order(abs(Z_matrix), decreasing=T)])
            plotDF$datetime <- rep(ymd(thisDate), numInds) 
            p2 <- ggplot(plotDF, aes(x=datetime, y=est)) + geom_ribbon(aes(ymin=confLow, ymax=confUp), fill="lightgray") + 
                geom_point(aes(x=datetime, y=origDat), color="blue") +
                geom_line(color="black", linewidth=1) +
                facet_wrap(~indName, ncol=ceiling(sqrt(numInds))) +
                labs(title=paste0("expectation, guildSEMnode = ", thisGuildSEMnode, ", guild = ", thisGuild), x="", y="Estimate") +
                theme_light() + theme(strip.text=element_text(size=6), axis.text.x=element_text(size=6))
            plotList[[length(plotList)+1]] <- p2
            
            patchworkPlotIndex <- patchworkPlotIndex + 1
            patchworkPlotList[[patchworkPlotIndex]] <- p1 / p2
            
            # Following along with Noble's code (makeDFAs.R):
            Z_est <- coef(fit, type="matrix")$Z
            if(numFactors>1) {
                # Get the inverse of the rotation matrix
                H_inv <- varimax(Z_est)$rotmat
                # Rotate factor loadings
                Z_rot = Z_est %*% H_inv
                # Rotate processes
                processes = solve(H_inv) %*% fit$states
                
                # Save the loadings
                loadingsDF <- data.frame(indicator=row.names(thisClusData), Z_est=as.numeric(Z_est), Z_rot=as.numeric(Z_rot))
            } else {
                processes <- fit$states
                
                # Save the loadings
                loadingsDF <- data.frame(indicator=row.names(thisClusData), Z_est=as.numeric(Z_est))
            }
            
            loadingsDF$guild <- thisGuild
            loadingsDF$guildSEMnode <- thisGuildSEMnode
            loadingsDF <- loadingsDF[, c("guild", "guildSEMnode", "indicator", "Z_est")]
            
            # Assign time series with |loading|<threshold to new guild
            loadingsDF$newGuild <- loadingsDF$guild
            newGuildIndex <- 0
            for(i in 1:nrow(loadingsDF)) {
                if(abs(loadingsDF[i, "Z_est"])<minLoading) {
                    loadingsDF[i, "newGuild"] <- paste0(loadingsDF[i, "guild"], "_", newGuildIndex)
                    newGuildIndex <- newGuildIndex + 1
                }
            }
            
            loadingsList[[length(loadingsList)+1]] <- loadingsDF

            for(j in 1:numFactors) {
                thisDFA <- data.frame(SEMlatent=thisGuildSEMnode, dendID=thisDendID, dendClusID=thisDendClusID, guild=thisGuild,
                                      shortName=paste(thisGuild, paste0("DFA", j), sep="_"), date=thisDate, finalVal=processes[j, ],
                                      minYear=thisMinYear)
                clusDataDFAlist[[length(clusDataDFAlist)+1]] <- thisDFA
            }
            
            fit$guildSEMnode <- thisGuildSEMnode
            fit$guild <- thisGuild
            fitList[[length(fitList)+1]] <- fit
            
        } else {
            thisStraggler <- gsub("_\\d+", "", thisGuild)
            thisGroup <- data.frame(guildSEMnode=thisGuildRow$guildSEMnode, guild=thisStraggler, shortName=unique(thisClusData$shortName),
                                    type="straggler")
            groupList[[length(groupList)+1]] <- thisGroup
            
            # Use pre-calculated smoothed data
            thisDate <- thisClusData$date
            thisDFA <- data.frame(SEMlatent=thisGuildSEMnode, dendID=thisDendID, dendClusID=thisDendClusID, guild=thisGuild,
                                  shortName=paste(thisGuild, "smoothed", sep="_"), date=thisDate, finalVal=thisClusData$smoothed,
                                  minYear=thisMinYear)
            clusDataDFAlist[[length(clusDataDFAlist)+1]] <- thisDFA
            
            # Plot smoothed time series
            thisClusData <- thisClusData %>% dplyr::select(shortName, date, finalVal) %>%
                pivot_wider(id_cols=date, names_from=shortName, values_from=finalVal)
            thisClusData <- thisClusData %>% select(-date) %>% mutate_all(~(scale(.) %>% as.vector))
            
            shortName <- clusData %>% filter(guild==thisGuild, guildSEMnode==thisGuildSEMnode) %>% distinct(shortName) %>% select(shortName) %>% pull(shortName)
            smoothedIndicatorsList[[length(smoothedIndicatorsList)+1]] <- data.frame(DFAname=paste(thisGuild, "smoothed", sep="_"), rankedIndicators=shortName)
            plotDF <- thisDFA[, c("date", "finalVal")]
            plotDF$raw <- as.numeric(unlist(thisClusData))
            
            p <- ggplot(plotDF) + geom_point(aes(x=date, y=raw)) + 
                geom_line(aes(x=date, y=finalVal)) + 
                theme_light() + backgroundTheme + theme(plot.subtitle=element_text(size=8))
            #ggsave(file=file.path(outputDir, paste0("smoothed_", shortName, ".png")), width=6, height=4)
            plotList[[length(plotList)+1]] <- p
            
            if(thisStraggler!=prevStraggler) {
                p <- p + labs(title=paste0("smoothed, guildSEMnode = ", thisGuildSEMnode, ", guild = ", thisStraggler), 
                              subtitle=shortName, x="", y="")
                patchworkPlotIndex <- patchworkPlotIndex + 1
                patchworkPlotList[[patchworkPlotIndex]] <- p + plot_layout(ncol=3)
            } else {
                p <- p + labs(title="", subtitle=shortName, x="", y="")
                patchworkPlotList[[patchworkPlotIndex]] <- patchworkPlotList[[patchworkPlotIndex]] + p
            }
            prevStraggler <- thisStraggler
        }
    }
    loadings <- bind_rows(loadingsList)
    write_xlsx(loadings, file.path(outputDir, "loadings.xlsx"))
    
    smoothedIndicators <- bind_rows(smoothedIndicatorsList)

    groups <- bind_rows(groupList)
    
    # Plot summary of number of indicators per DFA, number of stragglers, etc.
    groups <- groups %>% group_by(guildSEMnode, guild, type) %>% summarize(count=n(), .groups="drop") %>% 
        arrange(guild)
    p <- ggplot(groups, aes(x=guild, y=count, fill=type, label=count)) + geom_bar(stat="identity") + 
        geom_text(position=position_stack(vjust=0.5)) +
        scale_fill_brewer(palette="Set2") +
        labs(title="Guild splits", y="number of indicators", x="") +
        theme_light() + theme(axis.text.x=element_text(angle=90, hjust=1))
    patchworkPlotList[[length(patchworkPlotList)+1]] <- p
    
    thisTheme <- theme_light() + backgroundTheme + theme(plot.title=element_text(size=6))
    p <- autoplot(fitList[[1]], plot.type="xtT", fig.notes=F) + 
        labs(title=paste0(fitList[[1]]$guildSEMnode, ", ", fitList[[1]]$guild), x="", y="") + thisTheme
    if(length(fitList)>1) {
        for(i in 2:length(fitList)) {
            p <- p + autoplot(fitList[[i]], plot.type="xtT", fig.notes=F) + 
                labs(title=paste0(fitList[[i]]$guildSEMnode, ", ", fitList[[i]]$guild), x="", y="") + thisTheme
        }
    }
    patchworkPlotList[[length(patchworkPlotList)+1]] <- p
        
    pdf(file.path(outputDir, "allDFAs.pdf"))
    invisible(lapply(patchworkPlotList, print))
    dev.off()
}

clusDataDFA <- bind_rows(clusDataDFAlist)
saveRDS(clusDataDFA, file=file.path(outputDir, "clusDataDFA.rds"))

saveRDS(plotList, file=file.path(outputDir, "plotList.rds"))

# Create a new guilds file
guilds <- read.csv(guildFile)

newGuilds <- loadings %>% select(indicator, newGuild) %>% rename(shortName=indicator)
guilds <- left_join(guilds, newGuilds, by="shortName")
write.csv(guilds, file.path(outputDir, "newGuilds.csv"), row.names=F)

###########################################################################
# Ensure that the direction of the DFA matches the indicator with the highest loading
for(thisGuild in unique(loadings$guild)) {
    thisLoadings <- loadings %>% filter(guild==thisGuild)
    
    thisSign <- sign(thisLoadings$Z_est[which.max(abs(thisLoadings$Z_est))])
    cat(thisGuild, " sign:", thisSign, "\n")
}

# Create combined names from concatenated list of names ordered by absolute value of loading
DFAnames <- loadings %>% mutate(DFAname=paste0(guild, "_DFA1")) %>% arrange(guild, desc(abs(Z_est))) %>%
    group_by(guild, DFAname) %>% summarize(rankedIndicators=paste(indicator, collapse = " - "), .groups="drop") %>% 
    select(DFAname, rankedIndicators) %>% mutate(type="DFA")
smoothedIndicators$type <- "smoothed"
DFAnames <- bind_rows(DFAnames, smoothedIndicators) %>% arrange(DFAname)
write.csv(DFAnames, file.path(outputDir, "rankedIndicators.csv"), row.names=F)




