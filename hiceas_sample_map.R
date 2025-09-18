#HICEAS sampling map
#wd#####
setwd("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/")
#libraries#########
library(dplyr)
library(tidyverse)
library(stringr)
library(kableExtra)
library(ggOceanMaps)
library(marmap)
library(ggmap)
library(shape)
library(sf)
library(gpkg)
library(dplyr)
library(tidyr)
library(TSP)
library(readxl)
library(sp) #already in raster
library(raster)# error
library(ggplot2)
library(scales)
#library(rgdal)
library(marmap)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
#library(rgeos)
library(viridis)
library(raster)
library(ggnewscale)
library(tidyverse)
library(ggrepel)
#summary table#####
larv<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/SE2303 ID Log - ID&Length.csv")
larv2<-larv%>%
  filter(SORTING.PROCEDURE==1)%>%
  mutate(total_n=as.numeric(total))%>%
  group_by(station)%>%
  summarize(sum(total_n, na.rm=T))

kable(larv2,
      col.names =c("station", "total larvae per station")) %>%#,caption = "Total number of larvae from stations where all fish were sorted to family level and counted")
     kable_styling(latex_options = "striped",table.env='table*')


larv_uku<-larv%>%
  filter(genus=="Aprion")%>%
  mutate(total_n=as.numeric(total))%>%
  group_by(station)%>%
  summarize(sum(total_n, na.rm=T))

kable(larv_uku,
      col.names =c("station", "total larvae per station")) %>%#,caption = "Total number of larvae from stations where all fish were sorted to family level and counted")
  kable_styling(latex_options = "striped",table.env='table*')

#coordinate conversion####
ik<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/IKMT Log - Sheet1.csv")
#convert ship coordinates into DMS then to DD
ik<-ik %>% mutate(lat_deg= as.numeric(substr(Lat.In,1,2)))#http://127.0.0.1:40641/graphics/4ee7ccc4-d91e-43e0-8b31-d427019a6d2c.png
ik<-ik %>% mutate(lon_deg= as.numeric(substr(Lon.In,1,4)))
ik<-ik %>% mutate(lat_minsec= as.numeric(substr(Lat.In,3,10)))
ik<-ik %>% mutate(long_minsec= as.numeric(substr(Lon.In,5,10)))
ik<-ik %>% mutate(lat_dd= lat_deg+(lat_minsec/60))
ik<-ik %>% mutate(lon_dd= lon_deg+(long_minsec/60))
inventory<-read_excel("Data/SE2303 ID Log.xlsx", sheet="Inventory and Log")
ik2<-left_join(ik, inventory, by="Station")
#add in eDNA sites so far for funsies#
edan<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/eDNA Log - Cast_info (1).csv")

#cw's mapping codefor base map####
world<-ne_countries(scale="medium", returnclass = "sf")
oahu_raster <- raster(file.path("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/hi_eez_extract_grid.tiff"))
#from here: https://www.ncei.noaa.gov/maps/bathymetry/
#settings ETPO_2022(Bedrock, 15 arc seconds)
df <- fortify(as.bathy(oahu_raster))
#convert to 360
oahu_df<-df%>%mutate(x=ifelse(df$x < 0, df$x + 360, df$x))
load("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Hawaii.RData")
oahu_df <- fortify(Hawaii)
oahu_map <- ggplot(data = world) +
  geom_raster(data = oahu_df, aes(x = x, y = y, fill = z)) +
  labs(fill = "Depth (m)") +
  scale_fill_gradient(high = "lightskyblue1", low = "cornflowerblue", limits = c(-10000, 1000)) +# scale_fill_continuous(labels = scales::label_number(scale = 1000, suffix = "k")) +
  geom_sf() +
  coord_sf(xlim = c(206, 177), ylim = c(18, 32), label_graticule = "SW") + #THE FIX
  theme_bw() +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")
#http://127.0.0.1:19395/graphics/plot_zoom_png?width=1920&height=1027
#WIP: add in monument boundary######
library(sf)
PMNM<-read_sf("Data/PMNM_Shp_File/hi_noaa_nwhi_papahanaumokuakea.shp")
EEZ<-PMNM$geometry
#for some reason it is split into east and west of dateline, so need to combine 
EEZ_Plot1<-as.data.frame(EEZ[[1]][[1]])
EEZ_Plot2<-as.data.frame(EEZ[[1]][[2]])
EEZ_Plot<-rbind(EEZ_Plot1, EEZ_Plot2)
EEZ_Plot$Lon<-NULL
EEZ_Plot$Lon<-ifelse(EEZ_Plot$X1<0, EEZ_Plot$X1+360, EEZ_Plot$X1)

#connecting the dots in an efficient manner
xytsp <- ETSP(data.frame(EEZ_Plot$Lon, EEZ_Plot$X2))
colnames(xytsp) <- c("Lon", "Lat")
xytour <- solve_TSP(xytsp)
re_ordered_EEZ <- EEZ_Plot[xytour, ]
merp<-oahu_map+geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)

#mapping######
#FISH TIME#

#QCing tows and making a smaller, easier to work with data frame
ik_lite<-ik2%>%
  mutate(depth_m=as.numeric(Max.Depth))%>%
  mutate(mean_flow=((Flowmeter.1.End-Flowmeter.1.Start)+
                      (Flowmeter.2.End-Flowmeter.2.Start))/2,
         na.rm=T)%>%
  filter(depth_m>30)%>%
  filter(Gear.x=="6' IKMT")%>%
  mutate(datey=ymd(Date))%>%
  mutate(dec_date=julian(datey))%>%
  mutate(month=month(datey))%>%
  dplyr::select(Station,lat_dd,lon_dd,mean_flow, month, datey, dec_date)
#structuring fish ids
SE2303_ID<-read_excel("Data/SE2303 ID Log.xlsx", sheet = "ID&Length")
SE2303_ID$total<-as.numeric(SE2303_ID$total)
SE2303_ID$`1_over_split_fraction`<-as.numeric(SE2303_ID$`1_over_split_fraction`)

#unite dfs
se23<-left_join(SE2303_ID, ik_lite, join_by("station"=="Station"))

#density
se23<-se23%>%
  mutate(density=(total/mean_flow))

#plot fish####
##billfish####
fish_taxa="Istiophoridae"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)

##density
metric<-"density"
ggplot()+geom_point(data=df, mapping=aes(x=dec_date, y=density))+
  labs(title =paste(fish_taxa,metric,"by Julian day"))+
  geom_smooth(method="loess")

plot<-ggplot()+geom_point(data=df, mapping=aes(x=datey, y=density))+
  labs(title =paste(fish_taxa,metric,"by Date"))
#ggsave("plot",filename=paste("Plot Figures/IKMT_Tows_Plot_HICEAS_",fish_taxa,"_",metric,".png"), height=5, width=8, units="in", res=300)

ggplot()+geom_density(data=df, aes(x=datey),alpha=0.3)+
  labs(title =paste(fish_taxa,metric,"by Date"))

#Pomfrets#######
fish_taxa="Bramidae"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)

##density
metric<-"density"
ggplot()+geom_point(data=df, mapping=aes(x=dec_date, y=density))+
  labs(title =paste(fish_taxa,metric,"by Julian day"))+
  geom_smooth(method="loess")

plot<-ggplot()+geom_point(data=df, mapping=aes(x=datey, y=density))+
  labs(title =paste(fish_taxa,metric,"by Date"))
#ggsave("plot",filename=paste("Plot Figures/IKMT_Tows_Plot_HICEAS_",fish_taxa,"_",metric,".png"), height=5, width=8, units="in", res=300)

ggplot()+geom_density(data=df, aes(x=datey),alpha=0.3)+
  labs(title =paste(fish_taxa,metric,"by Date"))
###mackerel######
fish_taxa="Gempylidae"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)

##density
metric<-"density"
ggplot()+geom_point(data=df, mapping=aes(x=dec_date, y=density))+
  labs(title =paste(fish_taxa,metric,"by Julian day"))+
  geom_smooth(method="loess")

plot<-ggplot()+geom_point(data=df, mapping=aes(x=datey, y=density))+
  labs(title =paste(fish_taxa,metric,"by Date"))
#ggsave("plot",filename=paste("Plot Figures/IKMT_Tows_Plot_HICEAS_",fish_taxa,"_",metric,".png"), height=5, width=8, units="in", res=300)

ggplot()+geom_density(data=df, aes(x=datey),alpha=0.3)+
  labs(title =paste(fish_taxa,metric,"by Date"))
###mahi#######
fish_taxa="Coryphaenidae"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)

##density
metric<-"density"
ggplot()+geom_point(data=df, mapping=aes(x=dec_date, y=density))+
  labs(title =paste(fish_taxa,metric,"by Julian day"))+
  geom_smooth(method="loess")

plot<-ggplot()+geom_point(data=df, mapping=aes(x=datey, y=density))+
  labs(title =paste(fish_taxa,metric,"by Date"))
#ggsave("plot",filename=paste("Plot Figures/IKMT_Tows_Plot_HICEAS_",fish_taxa,"_",metric,".png"), height=5, width=8, units="in", res=300)

ggplot()+geom_density(data=df, aes(x=datey),alpha=0.3)+
  labs(title =paste(fish_taxa,metric,"by Date"))

###nehu###########
fish_taxa="Engraulidae"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)

##density
metric<-"density"
ggplot()+geom_point(data=df, mapping=aes(x=dec_date, y=density))+
  labs(title =paste(fish_taxa,metric,"by Julian day"))+
  geom_smooth(method="loess")

plot<-ggplot()+geom_point(data=df, mapping=aes(x=datey, y=density))+
  labs(title =paste(fish_taxa,metric,"by Date"))
#ggsave("plot",filename=paste("Plot Figures/IKMT_Tows_Plot_HICEAS_",fish_taxa,"_",metric,".png"), height=5, width=8, units="in", res=300)

ggplot()+geom_density(data=df, aes(x=datey),alpha=0.3)+
  labs(title =paste(fish_taxa,metric,"by Date"))
#snappers######
fish_taxa="Lutjanidae"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)

##density
metric<-"density"
ggplot()+geom_point(data=df, mapping=aes(x=dec_date, y=density))+
  labs(title =paste(fish_taxa,metric,"by Julian day"))+
  geom_smooth(method="loess")

plot<-ggplot()+geom_point(data=df, mapping=aes(x=datey, y=density))+
  labs(title =paste(fish_taxa,metric,"by Date"))
#ggsave("plot",filename=paste("Plot Figures/IKMT_Tows_Plot_HICEAS_",fish_taxa,"_",metric,".png"), height=5, width=8, units="in", res=300)

#Etelinae
fish_taxa="Etelinae"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(subfamily==fish_taxa)
ggplot()+geom_point(data=df, mapping=aes(x=dec_date, y=density))+
  labs(title =paste(fish_taxa,metric,"by Julian day"))+
  geom_smooth(method="loess")

ggplot()+geom_density(data=df, aes(x=datey),alpha=0.3)+
  labs(title =paste(fish_taxa,metric,"by Date"))

ggplot()+geom_density(data=df, aes(x=datey, fill=genus),alpha=0.3)+
  labs(title =paste(fish_taxa,metric,"by Date"))

#uku##
#Etelinae
fish_taxa="Aprion"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(genus==fish_taxa)
ggplot()+geom_point(data=df, mapping=aes(x=dec_date, y=density))+
  labs(title =paste(fish_taxa,metric,"by Julian day"))+
  geom_smooth(method="loess")
#Tunas##########
fish_taxa="Scombridae"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)

##density#
metric<-"density"
ggplot()+geom_point(data=df, aes(x=dec_date, y=density))+
  labs(title =paste(fish_taxa,metric,"by Julian day"))

ggplot()+geom_density(data=df, aes(x=datey, fill=genus),alpha=0.3)+
  labs(title =paste(fish_taxa, "genus",metric,"by Date"))
#ggsave("plot",filename=paste("Plot Figures/IKMT_Tows_Plot_HICEAS_",fish_taxa,"_",metric,".png"), height=5, width=8, units="in", res=300)


##map fish#######
#plot ALL possible stations#
iks_map<-oahu_map+geom_point(data=ik_lite, mapping=aes(y=lat_dd, x=lon_dd), shape=1, size=3, color="red")

df<-se23%>%
  filter(IS_MUS==T)
iks_map+geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density, color=family), shape=19)+

fish_taxa="Istiophoridae"
metric="density"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)
iks_map+geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density), shape=19, color="red")+
  labs(title =paste(fish_taxa,metric))

fish_taxa="Scombridae"
metric="density"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)
iks_map+geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density), shape=19, color="red")+
  labs(title =paste(fish_taxa,metric))+facet_grid(~genus)

fish_taxa="Coryphaenidae"
metric="density"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)
iks_map+geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density), shape=19, color="red")+
  labs(title =paste(fish_taxa,metric))

fish_taxa="Bramidae"
metric="density"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)
iks_map+geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density), shape=19, color="red")+
  labs(title =paste(fish_taxa,metric))

fish_taxa="Engraulidae"
metric="density"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)
iks_map+geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density), shape=19, color="red")+
  labs(title =paste(fish_taxa,metric))

fish_taxa="Gempylidae"
metric="density"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(family==fish_taxa)
iks_map+geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density), shape=19, color="red")+
  labs(title =paste(fish_taxa,metric))
