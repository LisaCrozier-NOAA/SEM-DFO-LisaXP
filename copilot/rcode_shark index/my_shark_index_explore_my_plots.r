shark_dat<-read.csv(file.path("data_Lisa/shark_wide.csv"))%>% 
  clean_names() %>%
  select(year,goa_pacific_sleeper_shark,goa_salmon_shark)
names(shark_dat)

ishark<-function(Nshark, dT,Tscale){
  Mt=2^(dT/10)
  Ot=plogis(Tscale) * 2
  ishark=Nshark*Mt*Ot
  return(ishark)
}



sst_dat<-read.csv(file.path("data_Lisa/goa_prey_clim_raw_trends_avg.csv"),row.names = 1) %>% 
  clean_names() %>%
  select(year,sst_egoa_coastwatch_junjulaug,sst_wgoa_coastwatch_junjulaug,swln_temp_fall_176to226m) %>%
  mutate(dT_sst_egoa=sst_egoa_coastwatch_junjulaug-mean(sst_egoa_coastwatch_junjulaug),
         dT_sst_wgoa=sst_wgoa_coastwatch_junjulaug-mean(sst_wgoa_coastwatch_junjulaug),
         dT_SewardLine_176to226m=swln_temp_fall_176to226m-mean(swln_temp_fall_176to226m),
         sc_sst_egoa=scale(sst_egoa_coastwatch_junjulaug),
         sc_sst_wgoa=scale(sst_wgoa_coastwatch_junjulaug),
         sc_SewardLine_176to226m=scale(swln_temp_fall_176to226m)        
         )

alldat<-full_join(shark_dat,sst_dat,by="year") %>%
  mutate(
        ishark_sst_egoa=ishark((goa_pacific_sleeper_shark) ,dT_sst_egoa, sc_sst_egoa),
        ishark_sst_wgoa=ishark((goa_pacific_sleeper_shark) ,dT_sst_wgoa, sc_sst_wgoa),
        ishark_SewardLine_176to226m=ishark((goa_pacific_sleeper_shark) ,dT_SewardLine_176to226m, sc_SewardLine_176to226m)
  )%>%
  mutate(
    ln_ishark_sst_egoa=ishark(log(goa_pacific_sleeper_shark) ,dT_sst_egoa, sc_sst_egoa),
    ln_ishark_sst_wgoa=ishark(log(goa_pacific_sleeper_shark) ,dT_sst_wgoa, sc_sst_wgoa),
    ln_ishark_SewardLine_176to226m=ishark(log(goa_pacific_sleeper_shark) ,dT_SewardLine_176to226m, sc_SewardLine_176to226m)
  )

ishark_dat<-alldat %>% select(goa_pacific_sleeper_shark,sst_wgoa_coastwatch_junjulaug)

names(sst_dat)
(alldat %>% select(year,goa_pacific_sleeper_shark ,ishark_sst_egoa,ishark_sst_wgoa,ishark_SewardLine_176to226m))
round(cor(sst_dat[,2:4]),3)

par(mfrow=c(1,1))

#Traw---------
plot(sst_egoa_coastwatch_junjulaug~year,data=sst_dat, lwd=2,col=2,type='l',ylim=c(0,15),bty='l')
    lines(sst_wgoa_coastwatch_junjulaug~year,data=sst_dat,lwd=2,col=3)
    lines(swln_temp_fall_176to226m~year,data=sst_dat,lwd=2,col=4)
    legend("bottomleft",legend=c("sst_egoa_coastwatch_junjulaug","sst_wgoa_coastwatch_junjulaug","swln_temp_fall_176to226m"),
           lty=1,lwd=2,col=c(2:4),bty='n')
    mtext(side=3,outer=T,"Traw",line=-3)
    mtext(side=3,outer=T,"Corr: 0.3 - 0.9",line=-5)
    

#dT   ------------- 
plot(dT_sst_egoa~year,data=sst_dat,ylab="Raw temperature differential",
     lwd=2,col=2,type='l',ylim=c(-3,2),bty='l')
  mtext(side=3,outer=T,"dT for M(t)",line=-3)
  mtext(side=3,outer=T,"Range: +/- 1.5 (sst) or \n +/- 0.75 (at depth)",line=-5)
  lines(dT_sst_wgoa~year,data=sst_dat,lwd=2,col=3)
    lines(dT_SewardLine_176to226m~year,data=sst_dat,lwd=2,col=4)
    legend("bottomleft",legend=c("sst_egoa_coastwatch_junjulaug","sst_wgoa_coastwatch_junjulaug","swln_temp_fall_176to226m"),
           lty=1,lwd=2,col=c(2:4),bty='n')

#M(t) ---------
    mt<-function(dT){
      mt=2^(dT/10)
    }
    plot(mt(dT_sst_egoa)~year,data=sst_dat,ylab="Mt",
         lwd=2,col=2,type='l',ylim=c(0.8,1.2),bty='l')
    mtext(side=3,outer=T,"M(t)",line=-3)
    mtext(side=3,outer=T,"Range: +/- 0.9-1.1",line=-5)
    lines(mt(dT_sst_wgoa)~year,data=sst_dat,lwd=2,col=3)
    lines(mt(dT_SewardLine_176to226m)~year,data=sst_dat,lwd=2,col=4)
    legend("bottomleft",legend=c("sst_egoa_coastwatch_junjulaug","sst_wgoa_coastwatch_junjulaug","swln_temp_fall_176to226m"),
           lty=1,lwd=2,col=c(2:4),bty='n')
    
    lty=1,lwd=2,col=c(2:4),bty='n')

#O(t)------- 
    Ot<-function(Tscale){
      Ot=plogis(Tscale) * 2
      return(Ot)
    }
plot(ot(sc_sst_egoa)~year,data=sst_dat,ylab="Spatial overlap",
     lwd=2,col=2,type='l',ylim=c(-0.2,2),bty='l')
      mtext(side=3,outer=T,"O(t)",line=-3)
      mtext(side=3,outer=T,"Range: 0.2 - 1.7",line=-5)
      lines(ot(sc_sst_wgoa)~year,data=sst_dat,lwd=2,col=3)
      lines(ot(sc_SewardLine_176to226m)~year,data=sst_dat,lwd=2,col=4)
      legend("bottomleft",legend=c("sst_egoa_coastwatch_junjulaug","sst_wgoa_coastwatch_junjulaug","swln_temp_fall_176to226m"),
             lty=1,lwd=2,col=c(2:4),bty='n')
#M(t) * O(t)--------
plot(ot(sc_sst_egoa)*mt(dT_sst_egoa)~year,data=sst_dat,ylab="Spatial overlap",
           lwd=2,col=2,type='l',ylim=c(-0.2,2),bty='l')
      mtext(side=3,outer=T,"M(t) * O(t)",line=-3)
      mtext(side=3,outer=T,"Range: 0 - 2",line=-5)
      lines(ot(sc_sst_wgoa)*mt(dT_sst_wgoa)~year,data=sst_dat,lwd=2,col=3)
      lines((ot(sc_SewardLine_176to226m) * mt(dT_SewardLine_176to226m))~year,data=sst_dat,lwd=2,col=4)
      legend("bottomleft",legend=c("sst_egoa_coastwatch_junjulaug","sst_wgoa_coastwatch_junjulaug","swln_temp_fall_176to226m"),
             lty=1,lwd=2,col=c(2:4),bty='n')
      
      
            
#Ishark-------
      ishark<-function(Nshark, dT,Tscale){
        Mt=2^(dT/10)
        Ot=plogis(Tscale) * 2
        ishark=Nshark * Mt * Ot
        return(ishark)
      }
      
      
plot(scale(goa_pacific_sleeper_shark) ~year,data=alldat,ylab="Shark index",
           lwd=2,col=1,type='l',ylim=c(-1,4),bty='l')
      lines(scale(ishark_sst_egoa)~year,data=alldat,lwd=2,col=2)
      lines(scale(ishark_sst_wgoa)~year,data=alldat,lwd=2,col=3)
      lines(scale(ishark_SewardLine_176to226m)~year,data=alldat,lwd=2,col=4)
      legend("topright",legend=c("Nshark","ishark_sst_egoa","ishark_sst_wgoa","ishark_swln_temp_fall_176to226m"),
             lty=1,lwd=2,col=c(1:4),bty='n')
      mtext(side=3,outer=T,"Raw shark count, and 3 indices",line=-3,cex=2)
      mtext(side=3,outer=T,"Range 0.5 - 2 prior to scaling",line=-4)
      
#Ishark_log-----
      
      plot(scale(log(goa_pacific_sleeper_shark)) ~year,data=alldat,ylab="Shark index",
           lwd=2,col=1,type='l',ylim=c(-2,4),bty='l')
   #   lines(scale(ln_ishark_sst_egoa)~year,data=alldat,lwd=2,col=2)
      lines(scale(ln_ishark_sst_wgoa)~year,data=alldat,lwd=2,col=3)
  #    lines(scale(ln_ishark_SewardLine_176to226m)~year,data=alldat,lwd=2,col=4)
      legend("topright",legend=c("ln_Nshark","ln_ishark_sst_wgoa"),
#             legend("topright",legend=c("ln_Nshark","ln_ishark_sst_egoa","ln_ishark_sst_wgoa","ln_ishark_swln_temp_fall_176to226m"),
                    lty=1,lwd=2,col=c(1:2),bty='n')
      mtext(side=3,outer=T,"Log transformed raw shark count, and best index",line=-3,cex=2)
      
      
      
      
    #Range: 
     round(apply(sst_dat,2,range,na.rm=T),2)
    # year sst_egoa_coastwatch_junjulaug sst_wgoa_coastwatch_junjulaug swln_temp_fall_176to226m
    # [1,] 1998                         10.88                          9.51                     4.74
    # [2,] 2021                         13.95                         12.33                     6.02
    # dT_sst_egoa dT_sst_wgoa dT_SewardLine_176to226m sc_sst_egoa sc_sst_wgoa sc_SewardLine_176to226m
    # [1,]       -1.63       -1.21                   -0.72       -2.15       -1.44                   -2.19
    # [2,]        1.43        1.62                    0.56        1.89        1.93                    1.69
    
    #Plot base index and all variants

tt=seq(-2,2,by=0.1)

ishark<-function(Nshark, dT,Tscale){
  Mt=2^(dT/10)
  Ot=plogis(Tscale) * 2
  ishark=Nshark*Mt*Ot
  return(ishark)
}

plot(tt,ishark(1,tt,tt),type='l',xlab="Temp")
abline(h=1)
legend("topleft",legend=c("dT=Tscale",
                          ), lty=1,col=1:)

lines(tt,ishark(1,tt,tt))