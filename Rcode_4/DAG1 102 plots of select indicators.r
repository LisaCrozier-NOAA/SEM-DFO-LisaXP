#Explore PreyAK for insights

# Van Doornik compared SEAK/WAOR cpue 2011-2015, and found that generall strongly correlated but 2011 lower in SEAK
# They blamed low productivity in SEAK: high salps that ate all the phyto, low chl-a, lots of picoplankton, 
# low Pcod and walleye pollock larvae collected in GOA
library(readxl)

#Look at loadings
path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_06_29_SEM_AKPred/shiftLisa_step3_26jun26/"
model_names <- c("DAG1A_long", "DAG1A_short", "DAG1B_long", "DAG1B_short", "DAG1C_long", "DAG1C_short")
#loadings<-read_xls(paste0(path, "DFA/loadings.xlsx"),sheet=1);loadings
loadings<-read.csv("LisaXP/Rcode_4/loadings_shiftLisa_step3_26jun26.csv",row.names=NULL);loadings
egoa.krill_amph.csv<-read.csv("~/Marine survival/Ferris-DFAIndicators-goa/egoa.krill_amph.csv",row.names=NULL);head(egoa.krill_amph.csv)

# var_lookup_NCC_AK<-read.csv("LisaXP/outputs_4/var_lookup_NCC_AK.SAR.csv",row.names=1);names(var_lookup_NCC_AK)
# guild.dfas1<-read.csv("LisaXP/outputs_4/guild.dfa.NCC.AK.csv",row.names=1);head(guild.dfas1)
# sem_master_data<-read.csv("LisaXP/outputs_4/sem_master_data.csv",row.names=1);head(sem_master_data)
# guild.dfas1$sarSR<-sem_master_data$sarSR

#Krill----------
#First plot in R, then ggplot for saving
          plot(1998:2021,sem_master_data$X02.krill,type='l',ylim=c(-3,3))
          lines(1998:2021,scale(sem_master_data$X12_DFA_biomassEuphShelfSum),col=2)
          lines(1998:2021,scale(egoa.krill_amph.csv$SECM.EUPH.DENS),col=4)
          legend("bottomleft",legend=c("NPL krill DFA","WCVI krill_amph DFA","EGOA Krill"),col=c(1,2,4),lty=1,bty="n")
          
          lines(1998:2021,sem_master_data$X09_DFA_HakeAge5Plus,col=2)
          
          lines(1998:2021,scale(sem_master_data$X16_SAR),col=8)
          lines(1998:2021,-1*scale(sem_master_data$X12_DFA_biomassEuphShelfSum),col=3)
          lines(1998:2021,scale(egoa.krill_amph.csv$SECM.HYP.AMPH),col=5)
          
          
          library(tidyverse)

        # 1. Combine vectors into a single data frame
        krill_data <- tibble(
          Year = 1998:2021,
          `NPL krill DFA`       = sem_master_data$X02.krill,
          `WCVI krill_amph DFA` = as.numeric(scale(sem_master_data$X12_DFA_biomassEuphShelfSum)),
          `EGOA Krill_Ferris2025`          = as.numeric(scale(egoa.krill_amph.csv$SECM.EUPH.DENS))
        )
        
        write.csv(krill_data,"LisaXP/data/krill_NPL_DFO_AK.csv")
        
        # 2. Pivot to long format and plot
        krill_plot <- krill_data %>%
          pivot_longer(
            cols = -Year, 
            names_to = "Indicator", 
            values_to = "Value"
          ) %>%
          ggplot(aes(x = Year, y = Value, color = Indicator)) +
          geom_line(linewidth = 1) +
          scale_color_manual(
            values = c(
              "NPL krill DFA"       = "black",
              "WCVI krill_amph DFA" = "red",
              "EGOA Krill"          = "blue"
            )
          ) +
          scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
          coord_cartesian(ylim = c(-3, 3)) +
          labs(
            x = "Year",
            y = "Standardized Anomaly",
            color = NULL
          ) +
          theme_classic() +
          theme(
            legend.position = c(0.2, 0.15),
            legend.background = element_blank(),
            legend.key = element_blank()
          )
        
        # Display the plot
        print(krill_plot)
        
        ggsave("LisaXP/data/krill_NPL_DFO_AK.png", krill_plot)

#Amphipods----------
        
        amphipod_data <- tibble(
          Year = 1998:2021,
          `JSOES sumPreyofPrey DFA`       = sem_master_data$X01_DFA_sumPreyOfPrey_planktonJun,
          `WCVI krill_amph DFA` = as.numeric(scale(sem_master_data$X12_DFA_biomassEuphShelfSum)),
          `EGOA amphipod_Ferris2025`          = as.numeric(scale(egoa.krill_amph.csv$SECM.HYP.AMPH))
        )
        
        amphipod_plot <- amphipod_data %>%
          pivot_longer(
            cols = -Year, 
            names_to = "Indicator", 
            values_to = "Value"
          ) %>%
          ggplot(aes(x = Year, y = Value, color = Indicator)) +
          geom_line(linewidth = 1) +
          scale_color_manual(
            values = c(
              "JSOES sumPreyofPrey DFA"       = "black",
              "WCVI krill_amph DFA" = "red",
              "EGOA amphipod_Ferris2025"          = "blue"
            )
          ) +
          scale_x_continuous(breaks = seq(1998, 2021, by = 2)) +
          coord_cartesian(ylim = c(-3, 3)) +
          labs(
            x = "Year",
            y = "Standardized Anomaly",
            color = NULL
          ) +
          theme_classic() +
          theme(
            legend.position = c(0.2, 0.75),
            legend.background = element_blank(),
            legend.key = element_blank()
          )
        
        # Display the plot
        print(amphipod_plot)
        
        write.csv(amphipod_data,"LisaXP/data/amphipods_JSOES_DFO_AK.csv")
        ggsave("LisaXP/data/amphipods_JSOES_DFO_AK.png", amphipod_plot)
        
        