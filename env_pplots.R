###chla, sst, ssh per station###
# libraries####
library(lubridate)
library(rerddap)
library(rerddapXtracto)
library(ggplot2)
library(viridis)
library(dplyr)
library(mgcv)
library(lunar)
#station data#####
#coordinate conversion##
ik<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/IKMT Log - Sheet1.csv")
#convert ship coordinates into DMS then to DD
ik<-ik %>% filter(Max.Depth>10) #eliminates neustons and s10s
ik<-ik %>% mutate(lat_deg= as.numeric(substr(Lat.In,1,2)))#http://127.0.0.1:40641/graphics/4ee7ccc4-d91e-43e0-8b31-d427019a6d2c.png
ik<-ik %>% mutate(lon_deg= as.numeric(substr(Lon.In,1,4)))
lat_minsec_string<-substr(ik$Lat.In,3,10)
ik <- ik %>%mutate(lat_minsec = as.numeric(gsub(" ", "", lat_minsec_string)))
long_minsec_string <- substr(ik$Lon.In, 5, 10)
ik <- ik %>%mutate(long_minsec = as.numeric(gsub(" ", "", long_minsec_string)))
ik<-ik %>% mutate(lat_dd= lat_deg+(lat_minsec/60))
ik<-ik %>% mutate(lon_dd= lon_deg+(long_minsec/60))
ik<-ik %>% mutate(date= ymd(Date,tz="HST"))
ik<-ik %>% mutate(tow_duration= Time.Out.of.Water-Time.In..Local.)


#sst data extraction####
legend_title<-"SST (°C)"
swchlInfo <- as.info("noaacrwsstDaily", url="https://coastwatch.noaa.gov/erddap/")

SST_Match <- rxtracto(swchlInfo, parameter = 'analysed_sst', 
                      xcoord = ik$lon_dd, ycoord = ik$lat_dd, tcoord =ik$date , 
                      xlen = .2, ylen = .2, progress_bar = TRUE)
ik$sst<-SST_Match$`mean analysed_sst`

#chla####
legend_title <-expression("Chlorophyll Concentration" * "\n" * "(mg" ~ m^{-3} ~ ")")

swchlInfo <- as.info("noaacwN20VIIRSchlaDaily", url="https://coastwatch.noaa.gov/erddap/")

chla_Match <- rxtracto(swchlInfo, parameter = 'chlor_a', 
                      xcoord = ik$lon_dd, ycoord = ik$lat_dd, tcoord =ik$date, 
                      zcoord = rep(0, length(ik$lon_dd)),
                      zName = "altitude",
                      xlen = .2, ylen = .2, progress_bar = TRUE)
ik$chla_daily<-chla_Match$`mean chlor_a`

#ssh########
legend_title <- "SSH Anomaly (m)"

swchlInfo <- as.info("noaacwBLENDEDsshDaily", url="https://coastwatch.noaa.gov/erddap/")
sla_Match <- rxtracto(swchlInfo, parameter = 'sla', 
                      xcoord = ik$lon_dd, ycoord = ik$lat_dd, tcoord =ik$date , 
                      xlen = .2, ylen = .2, progress_bar = TRUE)
ik$sla<-sla_Match$`mean sla`

#add in fish info#####
larv<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/SE2303 ID Log - ID&Length.csv")
ik2<-left_join(ik, larv,by=join_by("Station"=="station"))
se23<-ik2%>%
  mutate(lon_dd=ifelse(lon_dd < 0,lon_dd+ 360, lon_dd))%>%
  mutate(depth_m=as.numeric(Max.Depth))%>%
  mutate(mean_flow=((Flowmeter.1.End-Flowmeter.1.Start)+
                      (Flowmeter.2.End-Flowmeter.2.Start))/2,
         na.rm=T)%>%
  filter(depth_m>30)%>%
  mutate(tow_duration=as.duration(Time.Out.of.Water-Time.In..Local.))%>%
  filter(Gear.x=="6' IKMT")%>%
  mutate(month=month(date))%>%
  mutate(toyos=as.numeric(Tow.yo.Count))%>%
  mutate(total=as.numeric(total))%>%
  mutate(density=(total/mean_flow))%>%
  dplyr::select(Station,lat_dd,lon_dd,mean_flow, month, depth_m,total, toyos,tow_duration,
                sst,chla_daily,FLAG, family, subfamily,genus, species,density)

