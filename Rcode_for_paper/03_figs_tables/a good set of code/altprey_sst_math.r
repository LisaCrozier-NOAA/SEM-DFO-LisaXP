#the scaled max/min fxn only ranges from a-b, e.g., from 0.5-1 instead of 0-1 impact on predator
min_max_scale <- function(x, a = 0, b = 1) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep((a + b) / 2, length(x)))
  a + (b - a) * ((x - rng[1]) / (rng[2] - rng[1]))
}


set.seed(100)
n=10
dat_raw<-data.frame(year=1998:(1997+n),ssl=rnorm(n),herr=rnorm(n),cap=rnorm(n),sst=rnorm(n)) 
dat_maxmin<- dat_raw %>%
  mutate(ssl=min_max_scale(ssl,a=0,b=1))%>%
  mutate(sst=min_max_scale(sst,a=0,b=1))%>%
  mutate(herr=min_max_scale(herr,a=0,b=1))%>%
  mutate(cap=min_max_scale(cap,a=0,b=1))%>%
  mutate(sst_sc_0.5wt=min_max_scale(sst,a=0.5,b=1))%>%
  mutate(herrwt=1-min_max_scale(herr,a=0,b=1))%>%
  mutate(capwt=1-min_max_scale(cap,a=0,b=1))%>%
  mutate(herr_sc_0.5wt=1-min_max_scale(herr,a=0,b=1))%>%
  mutate(cap_sc_0.5wt=1-min_max_scale(cap,a=0,b=1))
  
  head(dat_maxmin)

# min_max_scale <- function(x) {
#   rng <- range(x, na.rm = TRUE)
#   if (rng[1] == rng[2]) return(rep(0.5, length(x)))
#   (x - rng[1]) / (rng[2] - rng[1])
# }


# dat_maxmin <- dat_raw %>%
#   mutate(across(where(is.numeric) & !matches("Year|year"), min_max_scale,a=0.5,b=1))


#let wt=1
dat_maxmin<-dat_maxmin %>%
  mutate(peff_herr=ssl*(herrwt),
         peff_cap=ssl*(capwt),
         peff_sst=ssl*(sst_sc_0.5wt),
         peff_sst_herr=ssl*(herrwt)*(sst_sc_0.5wt),
         peff_cap_herr=ssl*(herrwt)*(capwt)
  )

#maybe 0-1 is too extreme
#sst        
dat<-dat_maxmin %>% select(year,ssl,sst,peff_sst) %>% head(10)
        nn=ncol(dat)-1
matplot(dat$year,dat[,-1],lty=1:nn,col=1:nn,type='b')
        legend("topright",legend=paste(1:nn,names(dat)[-1]),lty=1:nn,col=1:nn,bty='n')
abline(v=2007,col=2)     


#raw conditions
dat<-dat_maxmin %>% select(year,ssl,sst,herr,cap)
nn=ncol(dat)-1
matplot(dat$year,dat[,-1],lty=1:nn,col=1:nn,type='b',lwd=2,main="raw conditions")
legend("topright",legend=paste(1:nn,names(dat)[-1]),lty=1:nn,col=1:nn,bty='n')



mytitle="all 3 effect sizes"
dat<-dat_maxmin %>% select(year,ssl,peff_sst,peff_herr,peff_cap)
nn=ncol(dat)-1
matplot(dat$year,dat[,-1],lty=1:nn,col=1:nn,type='b',lwd=2,main=mytitle)
legend("topright",legend=paste(1:nn,names(dat)[-1]),lty=1:nn,col=1:nn,bty='n')



#herr + sst
dat<-dat_maxmin %>% select(year,ssl,sst,herr,cap,peff_sst,peff_herr,peff_sst_herr)
nn=ncol(dat)-1
matplot(dat$year,dat[,-1],lty=1:nn,col=1:nn,type='b',lwd=2)
        legend("topright",legend=paste(1:nn,names(dat)[-1]),lty=1:nn,col=1:nn,bty='n')
        
        
        #herr
        dat<-dat_maxmin %>% select(year,ssl,herr,peff_herr)
        nn=ncol(dat)-1
        matplot(dat$year,dat[,-1],lty=1:nn,col=1:nn,type='b')
        legend("topright",legend=paste(1:nn,names(dat)[-1]),lty=1:nn,col=1:nn,bty='n')
        
        
#look at scale of each effect, regardless of ssl
        dat<-dat_maxmin %>% select(year,sst_sc_0.5 ,herr,cap) %>% mutate(herr=1-herr, cap=1-cap)
        nn=ncol(dat)-1
        matplot(dat$year,dat[,-1],lty=1:nn,col=1:nn,type='b',lwd=2)
        legend("topright",legend=paste(1:nn,names(dat)[-1]),lty=1:nn,col=1:nn,bty='n')
        
        
        #herr
        dat<-dat_maxmin %>% select(year,ssl,herr,peff_herr)
        nn=ncol(dat)-1
        matplot(dat$year,dat[,-1],lty=1:nn,col=1:nn,type='b')
        legend("topright",legend=paste(1:nn,names(dat)[-1]),lty=1:nn,col=1:nn,bty='n')
        
        
   
