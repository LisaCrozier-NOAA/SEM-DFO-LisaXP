library(dplyr)
library(tibble)



proj_dir     <- paste0(getwd(),"/Rcode_for_paper")
meta_dir     <- file.path(proj_dir, "metadata")
doug_dir     <- "C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/2026_10_23_SEM/shiftLisa_1_Lisa_22oct26"
crosswalk_names<-read.csv(file.path(meta_dir, "crosswalk_names.csv"), row.names = NULL)

path_doug<-"C:/Users/Lisa.Crozier/Documents/Marine survival/Doug results/"
doug_complete<-read.csv(file.path(path_doug,"2026_10_23_SEM/shiftLisa_1_Lisa_22oct26/DAG1A_long/completeness.csv")) %>%
  rename(year=date) %>% mutate(year=year(year));
head(doug_complete[,1:5])

doug_complete %>% select(contains("X08."))

#to change names in doug's file
# 1. Create a named vector: names are old column names (dfa_cols), values are new column names (LisaName)
name_map <- deframe(crosswalk_names %>% select(dfa_cols, LisaName))

# 2. Rename columns dynamically
doug_dat <- doug_complete %>%
  rename_with(~ ifelse(.x %in% names(name_map), name_map[.x], .x))

doug_dat %>% select(contains("X08"))


dat<-doug_dat %>% select(year,contains("X02"),contains("X04"),contains("X08"))
n<-ncol(dat)-1
matplot(dat$year,dat[,-1],type='b',col=1:n,lty=1:n)
legend('bottom',legend=paste0(1:n,"_",names(dat)[-1]),,col=1:n,lty=1:n,bty='n')

#seanettles replaced plankton in short model
      dat<-doug_dat %>% select(year,contains("X02"),contains("X01"))
      n<-ncol(dat)-1
      matplot(dat$year,dat[,-1],type='b',col=1:n,lty=1:n)
      legend('bottom',legend=paste0(1:n,"_",names(dat)[-1]),,col=1:n,lty=1:n,bty='n')
      
      cor(doug_dat$X02_DFA_SeaNettle,doug_dat$X01_DFA_sumPreyOfPrey_planktonJun,use="pair") #-0.635511


#marketsquid replaced forage fish in long model
      dat<-doug_dat %>% select(year,contains("X02"),contains("X01"))
      n<-ncol(dat)-1
      matplot(dat$year,dat[,-1],type='b',col=1:n,lty=1:n)
      legend('bottom',legend=paste0(1:n,"_",names(dat)[-1]),,col=1:n,lty=1:n,bty='n')
      
      cor(doug_dat$X02_DFA_SeaNettle,doug_dat$X01_DFA_sumPreyOfPrey_planktonJun,use="pair") #-0.635511

#birds gone completely
      
      