dat<-mu.igf %>%filter(FINAL_GSI=="Interior_Sp")
p1<-
  ggplot(data=dat,aes(x=Year,y=igf,color=Tag))+
  geom_point() +
  geom_smooth(method= "gam",formula = y ~ s(x, bs = "cs"), aes(fill=Tag) )+
  ggtitle("IGF ~ HW, Fall v Spring") +
  facet_wrap(~Month)
p1


dat0<-mu.igf %>%filter(FINAL_GSI=="Interior_Sp") %>%
  group_by(Year,Month,Tag) %>%
  summarise(igf=mean(igf),n=sum(n))


#dat<-junH.big
dat<-dat0 %>% filter
igf.plot <- 
  ggplot(data = dat, aes(x = Year, y = igf,color=Month)) +
  ggtitle("IGF ~ HW + May/Jun avg") +
  
  # --- ALL DATA (Blue Theme) ---
  geom_line(color = "grey20", alpha = 0.9) +
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
  facet_wrap(Month~Tag) +
  
  labs(size = "Sample Size (n)")

igf.plot
