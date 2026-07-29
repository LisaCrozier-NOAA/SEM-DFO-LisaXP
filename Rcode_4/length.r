rootDir <- "C:/Users/dougj/Documents/QEDA/NWFSC/ECOTRAN/programs/analyzeAKindices"
rootDir <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/analyzeAKindices"
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
igf<-read_xlsx(igfFilePath);summary(igf) #input
#interior spring and fall yearling chinook
names(igf)

igf<-igf %>% 
  filter (FINAL_GSI %in% c("Interior_Sp","Interior_F"),FINAL_GSI_PROB>=0.8,Age_ClassByLength=="yearling")%>%
  select (Month,Year,Tag,FINAL_GSI,LabLength_mm,FultonCondition,revised_IGF)

mu.igf<-mu<-igf %>%
  group_by(Year,Month,Tag,FINAL_GSI) %>%
  summarize(igf=mean(revised_IGF),n=n())

mu.igf.final<-mu.final<-mu %>%
  group_by(Year) %>%
  summarize(IGF_mu_2025=mean(igf))
mu.igf.final

p1<-
  ggplot(data=mu.igf,aes(x=Year,y=igf,color=Tag))+
  geom_point() +
  geom_smooth(method= "gam",formula = y ~ s(x, bs = "cs"), aes(fill=Tag) )+
  ggtitle("IGF ~ HW, Fall v Spring") +
  facet_wrap(~FINAL_GSI)
p1

dat<-mu.igf %>%filter(FINAL_GSI=="Interior_Sp")
p1<-
  ggplot(data=dat,aes(x=Year,y=igf,color=Tag))+
  geom_point() +
  geom_smooth(method= "gam",formula = y ~ s(x, bs = "cs"), aes(fill=Tag) )+
  ggtitle("IGF ~ HW, Fall v Spring") +
  facet_wrap(~Month)
p1

p2<-
  ggplot(data=mu.igf.final,aes(x=Year,y=IGF_mu_2025))+
  geom_point() +
  geom_smooth(method= "gam",formula = y ~ s(x, bs = "cs"))


p2


#1b. Growth -- Length ----------
len0<-read_xlsx(len0filePath);len0 #input
len0[len0$Tag=="unknown","Tag"]<-"YES" #assume "unknown" fish were hatchery origin

len<-len0 %>% 
  filter (FINAL_GSI %in% c("Interior_Sp"),FINAL_GSI_PROB>=0.8,Age_ClassByLength=="yearling",Month !="September")%>%
  select (Month,Year,Tag,LabLength_mm,FultonCondition)

mu.length<-len %>%
  group_by(Year,Month,Tag) %>%
  summarize(Lmu=mean(LabLength_mm,na.rm=T),Kmu=mean(FultonCondition,na.rm=T),n=n())
mu.length

mayW<-mu.length %>% filter(Month=="May",Tag=="NO")
mayH<-mu.length %>% filter(Month=="May",Tag=="YES")
junH<-mu.length %>% filter(Month=="June",Tag=="YES")
junH.big<-mu.length %>% filter(Month=="June",Tag=="YES",n>10)


mayW.plot <- ggplot(data = mayW, aes(x = Year, y = Lmu)) +
  ggtitle("Length ~ HW + May/Jun avg") +
  
  # --- ALL DATA (Blue Theme) ---
  geom_line(color = "grey70", alpha = 0.5) +
  geom_point(aes(size = n), color = "black", alpha = 0.4) +
  geom_smooth(
    method = "gam", formula = y ~ s(x, bs = "cs"), aes(weight = n),
    color = "darkblue", fill = "lightblue", alpha = 0.2
  ) +
  
  # --- HIGH SAMPLE DATA ONLY (Red Theme) ---
  geom_line(data = filter(mayW, n > 10), color = "firebrick", linewidth = 0.8) +
  geom_point(data = filter(mayW, n > 10), aes(size = n), color = "firebrick", alpha = 0.9) +
  geom_smooth(
    data = filter(mayW, n > 10), # <-- Drops low-n years entirely from this model
    method = "gam", formula = y ~ s(x, bs = "cs"),
    color = "firebrick", fill = "firebrick", alpha = 0.1
  ) +
  
  labs(size = "Sample Size (n)")

mayW.plot

mayH.plot <- ggplot(data = mayH, aes(x = Year, y = Lmu)) +
  ggtitle("Length ~ HW + May/Jun avg") +
  
  # --- ALL DATA (Blue Theme) ---
  geom_line(color = "grey70", alpha = 0.5) +
  geom_point(aes(size = n), color = "black", alpha = 0.4) +
  geom_smooth(
    method = "gam", formula = y ~ s(x, bs = "cs"), aes(weight = n),
    color = "darkblue", fill = "lightblue", alpha = 0.2
  ) +
  
  # --- HIGH SAMPLE DATA ONLY (Red Theme) ---
  geom_line(data = filter(mayH, n > 10), color = "firebrick", linewidth = 0.8) +
  geom_point(data = filter(mayH, n > 10), aes(size = n), color = "firebrick", alpha = 0.9) +
  geom_smooth(
    data = filter(mayH, n > 10), # <-- Drops low-n years entirely from this model
    method = "gam", formula = y ~ s(x, bs = "cs"),
    color = "firebrick", fill = "firebrick", alpha = 0.1
  ) +
  
  labs(size = "Sample Size (n)")

mayH.plot


junH.big<-mu.length %>% filter(Month=="June",Tag=="YES",n>10)
junH<-mu.length %>% filter(Month=="June",Tag=="YES")
dat<-junH
junH.plot <- 
  ggplot(data = dat, aes(x = Year, y = Lmu)) +
  ggtitle("Length ~ HW + May/Jun avg") +
  
  # --- ALL DATA (Blue Theme) ---
  geom_line(color = "grey70", alpha = 0.5) +
  geom_point(aes(size = n), color = "black", alpha = 0.4) +
  geom_smooth(
    method = "gam", formula = y ~ s(x, bs = "cs"), aes(weight = n),
    color = "darkblue", fill = "lightblue", alpha = 0.2
  ) +
  
  # --- HIGH SAMPLE DATA ONLY (Red Theme) ---
  geom_line(data = filter(dat, n > 10), color = "firebrick", linewidth = 0.8) +
  geom_point(data = filter(dat, n > 10), aes(size = n), color = "firebrick", alpha = 0.9) +
  geom_smooth(
    data = filter(dat, n > 10), # <-- Drops low-n years entirely from this model
    method = "gam", formula = y ~ s(x, bs = "cs"),
    color = "firebrick", fill = "firebrick", alpha = 0.1
  ) +
  
  labs(size = "Sample Size (n)")

junH.plot

junW<-mu.length %>% filter(Month=="June",Tag=="NO")
dat<-junW
junW.plot <- 
  ggplot(data = dat, aes(x = Year, y = Lmu)) +
  ggtitle("Length ~ HW + May/Jun avg") +
  
  # --- ALL DATA (Blue Theme) ---
  geom_line(color = "grey70", alpha = 0.5) +
  geom_point(aes(size = n), color = "black", alpha = 0.4) +
  geom_smooth(
    method = "gam", formula = y ~ s(x, bs = "cs"), aes(weight = n),
    color = "darkblue", fill = "lightblue", alpha = 0.2
  ) +
  
  # --- HIGH SAMPLE DATA ONLY (Red Theme) ---
  geom_line(data = filter(dat, n > 10), color = "firebrick", linewidth = 0.8) +
  geom_point(data = filter(dat, n > 10), aes(size = n), color = "firebrick", alpha = 0.9) +
  geom_smooth(
    data = filter(dat, n > 10), # <-- Drops low-n years entirely from this model
    method = "gam", formula = y ~ s(x, bs = "cs"),
    color = "firebrick", fill = "firebrick", alpha = 0.1
  ) +
  
  labs(size = "Sample Size (n)")

junW.plot


library(patchwork)

# 1. Update the title on the second plot so you can easily tell them apart!
mayW.plot <- mayW.plot + ggtitle("May (Tag: NO) - Length ~ W")
junH.plot <- junH.plot + ggtitle("June (Tag: YES) - Length ~ H")
mayH.plot <- mayH.plot + ggtitle("May (Tag: YES) - Length ~ H")
junW.plot <- junW.plot + ggtitle("June (Tag: NO) - Length ~ W")

# 2. Put them side-by-side and collect the legends so it looks clean
combined_plots <- (mayW.plot + mayH.plot + junW.plot + junH.plot ) + 
  plot_layout(guides = "collect")

# 3. Print the final layout
length.plots<-combined_plots

ggsave("LisaXP/outputs_4/length.plots.png",length.plots)


#older---------
p1.L<-
  ggplot(data=junH.big,aes(x=Year,y=Lmu,color=Tag))+
  geom_point() +
  geom_smooth(method= "gam",formula = y ~ s(x, bs = "cs"), aes(fill=Tag) )+
  ggtitle("Length ~ HW, Spring") +
  facet_wrap(~Month)
p1.L
# 
mayW<-
  ggplot(data=mayW,aes(x=Year,y=Lmu))+
  ggtitle("Length ~ HW + May/Jun avg") +
  geom_point(aes(size = n), alpha = 0.7) +
  geom_smooth(
    method = "gam", 
    formula = y ~ s(x, bs = "cs"), 
    aes(weight = n) # <-- This applies the weight to the model
  ) +
  labs(size = "Sample Size (n)")

mayW

mayW <- ggplot(data = mayW, aes(x = Year, y = Lmu)) +
  ggtitle("Length ~ HW + May/Jun avg") +
  # 1. Connect the points chronologically with a thin line
  geom_line(color = "grey50", alpha = 0.7) +
  # 2. Size the points by sample size
  geom_point(aes(size = n), alpha = 0.7) + 
  # 3. Weight the GAM curve by sample size
  geom_smooth(
    method = "gam", 
    formula = y ~ s(x, bs = "cs"), 
    aes(weight = n)
  ) +
  labs(size = "Sample Size (n)") 

mayW


junH<-
  ggplot(data=junH,aes(x=Year,y=Lmu))+
  ggtitle("Length ~ HW + May/Jun avg") +
  geom_line(color = "grey50", alpha = 0.7) +
  geom_point(aes(size = n), alpha = 0.7) +
  geom_smooth(
    method = "gam", 
    formula = y ~ s(x, bs = "cs"), 
    aes(weight = n) # <-- This applies the weight to the model
  ) +
  labs(size = "Sample Size (n)")

junH
