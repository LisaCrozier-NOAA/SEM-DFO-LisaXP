#Goal: generate a predator exposure index for salmon sharks and SSL in GOA

#1. Salmon sharks are opportunistic pelagic hunters. 
#When salmon are scarce or spread thin, salmon sharks switch to other abundant midwater targets:
#Gonatid Squid (Berryteuthis anonychus & Gonatus spp.): 
#The primary alternative prey for salmon sharks in the central GOA and offshore waters. 
#Squid populations fluctuate wildly with oceanographic conditions.
#Adult Pollock & Pacific Herring: Large schooling forage in the upper-to-mid water column.
#Sablefish (Anoplopoma fimbria): Juvenile and subadult sablefish occupying pelagic/slope waters.
#The risk term for salmon sharks ($R_{\text{SalmonShark}}$) can be structured as:$$R_{\text{SalmonShark}, t} = \frac{\text{SalmonShark\_CPUE}_t}{\text{Squid}_t + \text{Pollock}_t + \text{Herring}_t + \epsilon} \times \text{ThermalCompression}_t$$

#wSS=SS/pollock (or residual from the fishery?)
#HCI * wSS/squid + herring + sablefish (juv)

#HCI is the thermal modifier for increased metabolic rate
#However, temperature dictates where in the water column they encounter Chinook:
#In Cold/Normal Years: Salmon sharks forage across a wider depth band, and 
#juvenile/subadult salmon have more vertical space to evade them.
#In Warm/Heatwave Years: Surface waters exceed comfortable foraging temperatures (12-14C). 
#Both subadult Chinook and salmon sharks are squeezed into a narrow thermocline layer (40-80m), 
#dramatically increasing encounter rates.

#Salmon Sharks: almost entirely caught in the walleye pollock fishery

sem_master_data<-sem_master_data %>%
relocate(year)

png("LisaXP/data/sharks_raw.png")        
      n=which(grepl("Shark|shark",names(sem_master_data)));n
      matplot(sem_master_data$year,sem_master_data[,n],type='b',col=1:length(n))
      legend("topright",legend=paste(1:length(n),names(sem_master_data)[n]),col=1:n,bty='n',cex=0.75)
      mtext(side=3,"Sharks")
dev.off()

shark.dat<-sem_master_data[,which(grepl("year|Shark|shark|pollock",names(sem_master_data)))];head(shark.dat)
pollock

png("LisaXP/data/pollock_raw.png")        
        n=which(grepl("pollock",names(sem_master_data)));n
        matplot(sem_master_data$year,sem_master_data[,n],type='b',col=1:length(n))
        legend("topleft",legend=paste(1:length(n),names(sem_master_data)[n]),col=1:n,bty='n',cex=0.75)
        mtext(side=3,"Walleye pollock")
dev.off()


library(readxl)
catch<-read_excel(path="LisaXP/data/GOApollock.table1.23_2024assess.xlsx",skip=2) %>%
  rename(catch=`Catch\r\n(t)`,
         tb="3+ total\r\nbiomass\r\n(kt)...2") %>%
  select(Year,catch,tb) %>%
  mutate(across(everything(), parse_number))

head(catch)

head(catch)

shark<-left_join(shark.dat,catch,join_by(year==Year)) %>%
  mutate(across(-year, ~ as.vector(scale(.x))))
head(shark)

matplot(shark$year,shark%>% select(X15_salmonSharkGoA_predAK,X15_sharkCatchGoA_predAK,catch,tb))
legend("topleft",legend=paste(1:4,c("X15_salmonSharkGoA_predAK","X15_sharkCatchGoA_predAK","pollock catch","pollock total biomass")),col=1:4,bty='n')
       