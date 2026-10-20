
# library(gridExtra)
# library(grid)
# source("R code/multiplot.fxn.r")
library(readxl)
library(dplyr)
 library(ggplot2)

#1. Growth -- IGF and Length from JSOES cruise
#2. Abundance -- cpue June from JSOES cruise
#3. SAR -- pit tag data from DART 

#outputs:
#mu.igf -- May, Jun, H, W, IntF, IntS all separate
#mu.igf.final -- grand mean (after being grouped: May, Jun, H, W, IntF, IntS)
#mu.length -- both Length and Fulton's K --  May, Jun, H, W, IntS all separate
#mu.length.final -- grand mean (after being grouped: May, Jun, H, W, IntS)
#abundance -- sdmTMB output to be added to existing time series until that is updated
#sar.UCsprCh -- this is a new SAR time series to be added to existing through dfa, another set of populations


rootDir <- c("C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices","C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices")[1]
dataDir <- file.path(rootDir, "LisaDataProcessScripts 2025", "salmon")
outputDir <- file.path(rootDir, "LisaDataProcessScripts 2025")

igfFilePath <- file.path(dataDir, "2000_2022_IGF_MayJune_Subs_Year_Chin_Year_Coho_7.5.23.xlsx")
len0filePath <- file.path(dataDir, "1998_to_2022_Chinook_All_7.27.23.xlsx")
abundanceFilePath <- file.path(dataDir, "markus_interior_spring_chinook_yearlings_index.csv")
cpue0filePath <- file.path(dataDir, "Markus_Min_Trawl_CTD_Chl_Nuts_Thorson_Scheuerell_5.15.25_FINAL.xlsx")
newSARfilePath <- file.path(dataDir, "SARfromDART.UCspringCh.pitHW_poolbyyear20022024.csv")
oldSARfilePath <- file.path(dataDir, "sar.19982021.csv")
stomachFullnessPath <- file.path(rootDir, "LisaDataProcessScripts 2025", "raw_data", "StomachFullnesCKInterior_Crozier082125.xlsx")

#1a. Growth -- IGF ----------
    igf.orig<-read_xlsx(igfFilePath);summary(igf) #input
      #interior spring and fall yearling chinook
      names(igf.orig)
      
      igf<-igf.orig %>% 
        filter (FINAL_GSI == "Interior_Sp", FINAL_GSI_PROB>=0.8,Age_ClassByLength=="yearling")%>%
        select (Month,Year,Tag,FINAL_GSI,LabLength_mm,FultonCondition,revised_IGF)
      
      fit<-lm(revised_IGF~Year*Tag,data=igf);summary(fit) #significant linear trend in year, but Tag and YxTag NS
      fit<-lm(revised_IGF~Year*Month  ,data=igf);summary(fit) #significant linear trend in year, but Month and YxMonth NS

      
      mu.igf<-igf %>%
        group_by(Year,Month) %>%
        summarize(igf=mean(revised_IGF),n=n())
      
      mu.igf.final<-igf %>%
        group_by(Year) %>%
        summarize(IGF_mu_2025=mean(revised_IGF),n=n())
      print(mu.igf.final,n=Inf)
      
      p1<-
        ggplot(data=mu,aes(x=Year,y=igf,color=Tag))+
        geom_point() +
        geom_smooth(method= "gam",formula = y ~ s(x, bs = "cs"), aes(fill=Tag) )+
        ggtitle("IGF ~ HW, Fall v Spring") +
        facet_wrap(~FINAL_GSI)
      p1

      p2<-
        ggplot(data=mu.igf.final,aes(x=Year,y=IGF_mu_2025))+
        geom_point() +
        geom_smooth(method= "gam",formula = y ~ s(x, bs = "cs"))


      p2


#1b. Growth -- Length ----------
      len0<-read_xlsx(len0filePath) %>%
        rename(FishWt_gm="Fish Wt (gm)");len0 #input
      len0[len0$Tag=="unknown","Tag"]<-"YES" #assume "unknown" fish were hatchery origin
      
      len<-len0 %>% 
        filter (FINAL_GSI %in% c("Interior_Sp"),FINAL_GSI_PROB>=0.8,Age_ClassByLength=="yearling",Month !="September")%>%
        select (Month,Year,Tag,LabLength_mm,FishWt_gm,FultonCondition)
      
      fit<-lm(LabLength_mm~Year*Tag,data=len);summary(fit) #significant trend in year, Tag and YxTag
      fit<-lm(LabLength_mm~Year*Month*Tag  ,data=len);summary(fit) #all factors sign in general
      fit<-lm(LabLength_mm~Year  ,data=len %>% filter(Tag=="YES",Month=="May"));summary(fit) #Y sign
      fit<-lm(LabLength_mm~Year  ,data=len %>% filter(Tag=="NO",Month=="June"));summary(fit) #Y sign
      
      mu.length<-len %>%
        group_by(Year,Month,Tag) %>%
        summarize(Lmu=mean(LabLength_mm,na.rm=T),Wtmu=mean(FishWt_gm,na.rm=T),Kmu=mean(FultonCondition,na.rm=T),n=n())
      mu.length
      
      L_H_May<-mu.length %>% filter(Tag=="YES",Month=="May");print(L_H_May,n=Inf)
      L_W_May<-mu.length %>% filter(Tag=="NO",Month=="May");print(L_W_May,n=Inf) #several years small sample sizes
      L_H_Jun<-mu.length %>% filter(Tag=="YES",Month=="June");print(L_H_Jun,n=Inf)
      L_W_Jun<-mu.length %>% filter(Tag=="NO",Month=="June");print(L_W_Jun,n=Inf) #VERY SMALL SAMPLE SIZES --EXCLUDE
      
      mu.length %>% filter(n>10) 
      years_per_group <- mu.length %>%
        filter(n >= 10) %>%
        group_by(Tag, Month) %>%
        summarize(
          n_years = n_distinct(Year),
          .groups = "drop"
        )
      
      years_per_group
      #   Tag   Month n_years
      # 1 NO    June        9
      # 2 NO    May        14
      # 3 YES   June       20
      # 4 YES   May        20
       
      
      mu.length.final.1<-mu.length %>%
        filter (Month =="June", Tag=="YES")%>%
        group_by(Year) %>%
        summarize(Lmu_IntSprJunH_2025=mean(Lmu,na.rm=T))
      
      mu.length.final.2<-mu.length %>%
        filter (Month =="May", Tag=="NO")%>%
        group_by(Year) %>%
        summarize(Lmu_IntSprMayW_2025=mean(Lmu,na.rm=T))
      mu.length.final.2
      
      salmongrowth<-  full_join(mu.length.final.1,mu.length.final.2, by = "Year") #%>%  
      salmongrowth<-  full_join(salmongrowth,mu.igf.final, by = "Year") #%>% 
      salmongrowth

      p1.L<-
        ggplot(data=mu.length,aes(x=Year,y=Lmu,color=Tag))+
        geom_point() +
        geom_smooth(method= "gam",formula = y ~ s(x, bs = "cs"), aes(fill=Tag) )+
        ggtitle("Length ~ HW, Spring") +
        facet_wrap(~Month)
      p1.L
      # 
      # p2.L<-
      #   ggplot(data=mu.length.final,aes(x=Year,y=Lmu))+
      #   ggtitle("Length ~ HW + May/Jun avg") +
      #   geom_point() + 
      #   geom_smooth(method= "gam",formula = y ~ s(x, bs = "cs"))
      # p2.L
      # 
      # p1.K<-
      #   ggplot(data=mu.length,aes(x=Year,y=Kmu,color=Tag))+
      #   geom_point() + 
      #   geom_smooth(method= "gam",formula = y ~ s(x, bs = "cs"), aes(fill=Tag) )+
      #   ggtitle("Fulton condition ~ HW, Spring") +
      #   facet_wrap(~Month)
      # p1.K
      # 
      # p2.K<-
      #   ggplot(data=mu.length.final,aes(x=Year,y=Kmu))+
      #   ggtitle("Fulton condition ~  HW + May/Jun avg") +
      #   geom_point() + 
      #   geom_smooth(method= "gam",formula = y ~ s(x, bs = "cs"))
      # p2.K
      
#1c. Growth -- Stomach Fullness from Elizabeth Daly ----------
      stomachFullness <- readxl::read_excel(stomachFullnessPath, sheet=1,skip=1) %>%
        rename(StomFull.May="Average of %BW...2",StomFull.Jun="Average of %BW...5") %>%
        select(Year,StomFull.May,StomFull.Jun)
      
      # Follow the naming convention of the other data generated by this script
      names(stomachFullness)[2:3] <- paste0(names(stomachFullness)[2:3], "_2025")
      

#2. Condition 2 -----------------------
      #2a. sdmTMB estimate from Markus Min May 8, 2025
       abundance<-read.csv(file=abundanceFilePath, row.names = 1)%>%
         rename(Year=year,abund.sdmTMB.junIntSpr_2025=est) %>%
        select(Year,abund.sdmTMB.junIntSpr_2025)
      abundance
      
      
      #2b. raw cpue from survey data
       cpue0<-read_xlsx(cpue0filePath);cpue0 
       cpue0 <- cpue0 %>% 
         filter (Month == "June")%>%
         select (Year,`Chinook salmon_yearling_Interior_Sp`)
       cpue0
       
        cpue<-aggregate.data.frame(cpue0$`Chinook salmon_yearling_Interior_Sp`,list(cpue0$Year),mean,na.rm=T)
       names(cpue)<-c("Year","cpue_IntSprJunHW_2025")
       cpue
       
        salmoncpue<-  full_join(cpue,abundance, by = "Year") #%>%  
        salmon<-  full_join(salmongrowth,salmoncpue, by = "Year") #%>%  
        salmon
       
       
#3. SAR----------
       #SAR from DART for all Upper Columbia populations 2002-2024
       sar.UCsprCh<-read.csv(file=newSARfilePath, nrows = 22,row.names = NULL)%>% 
         select(year,meanSAR)%>%
         rename(Year=year,sar.UCsprCh_2025=meanSAR)
       sar.UCsprCh
       
       sar.old<-read.csv(file=oldSARfilePath, row.names = NULL) %>%
         rename(SAR_2025=SAR)
 
       salmonsar<-  full_join(sar.old,sar.UCsprCh, by = "Year") #%>%  
       salmon<-  full_join(salmon,salmonsar, by = "Year") #%>%  
       
       # Bring in stomach fullness
       salmon <- full_join(salmon, stomachFullness, by="Year")
       
       write.csv(salmon,file=file.path(outputDir, "salmon.datwide.csv"))
       