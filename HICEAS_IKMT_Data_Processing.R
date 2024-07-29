#plotting the IKMT tows


#plotting HICEAS sampling Stations 
library(ncdf4)
library(gsheet)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)
library(ggplot2)
library(ggnewscale)
library(sp)
library(sf)
library(lubridate)
library(readxl)
library(dplyr)
library(plyr)
library(viridis)
library(tidyr)

setwd("C:/JJS_Old_Lptp/HICEAS/IKMT_Processing")

IKMT_Log<-read_excel("IKMT Log.xlsx", sheet=1)

for (i in 1:nrow(IKMT_Log)){
  IKMT_Log$Lat_Minutes[i]<-substr(IKMT_Log$`Lat In`[i], 4,9)
  IKMT_Log$Lon_Minutes[i]<- substr(IKMT_Log$`Lon In`[i], 6,12)
  IKMT_Log$Lat_Deg[i]<-substr(IKMT_Log$`Lat In`[i], 1,2)}


for (i in 1:nrow(IKMT_Log)){
  IKMT_Log$Lon_Deg[i]<-ifelse(IKMT_Log$`Lon In`[i]<0, substr(IKMT_Log$`Lon In`[i], 1,4), substr(IKMT_Log$`Lon In`[i], 1,3))
}

IKMT_Log$Lat_DD<-as.numeric(IKMT_Log$Lat_Deg)+(as.numeric(IKMT_Log$Lat_Minutes)/60)

for (i in 1:nrow(IKMT_Log)){
  IKMT_Log$Lon_DD[i]<-ifelse(as.numeric(IKMT_Log$Lon_Deg[i])<0, as.numeric(IKMT_Log$Lon_Deg[i])+((as.numeric(IKMT_Log$Lon_Minutes[i])/60)*-1), as.numeric(IKMT_Log$Lon_Deg[i])+((as.numeric(IKMT_Log$Lon_Minutes[i])/60)))
}

for (i in 1:nrow(IKMT_Log)){
  IKMT_Log$Lon_East[i]<-ifelse(as.numeric(IKMT_Log$Lon_DD[i])<0, (360+IKMT_Log$Lon_DD[i]), IKMT_Log$Lon_DD[i])
}

# this is a lazy way, I'll fix this later
for (i in 1:nrow(IKMT_Log)){
  IKMT_Log$Lon_West[i]<-ifelse(as.numeric(IKMT_Log$Lon_DD[i])<0, IKMT_Log$Lon_DD[i], -180)
}

IKMT_Log$Leg<-as.factor(IKMT_Log$Leg)

#####################################################

SE2303_ID_Log<-read_excel("SE2303 ID Log.xlsx",sheet=4)

SE2303_ID_Log$station = as.factor(SE2303_ID_Log$station)

SE2303_Genus_Level<- SE2303_ID_Log %>%
  dplyr::group_by(genus,station, .drop=F) %>%
  dplyr::summarise(Total = sum(as.numeric(total), na.rm=T), Flag=sum(as.numeric(FLAG), na.rm=T))

merge_data<-merge(IKMT_Log,SE2303_Genus_Level,by.x = "Station",by.y = "station",all.y=TRUE)
genera<-unique(merge_data$genus)

##################################################
merge_data$`Flowmeter 1 End`<-as.numeric(merge_data$`Flowmeter 1 End`)
merge_data$`Flowmeter 2 End`<-as.numeric(merge_data$`Flowmeter 2 End`)
for (i in 1:nrow(merge_data)){
  merge_data$Maxflowcount[i]<-max(merge_data$`Flowmeter 1 End`[i],merge_data$`Flowmeter 2 End`[i],na.rm = T)
}
merge_data$volume<-(merge_data$Maxflowcount*.245)*2.94
merge_data$Abund_num_1000m3<-(merge_data$Total)/(merge_data$volume/1000)

#################################################
world <- ne_countries(scale=10,returnclass = "sf")#generate high res coastlines 

setwd("C:/JJS_Old_Lptp/HICEAS/IKMT_Processing/Genus_Level_Plots")
for(i in 1: length(genera)){
Genus_count<-merge_data[merge_data$genus==genera[i],]

png(paste0(genera[i],"_IKMT_Plots.png"), height=5, width=7, units="in", res=300)
p<-ggplot() +geom_point( data=Genus_count,aes(x=Lon_West, y=Lat_DD, size = Abund_num_1000m3, color = Abund_num_1000m3))+ geom_sf(data = world)+coord_sf(xlim = c(-180, -152), ylim = c(18,32))+scale_color_viridis()+ylab("Latitude")+xlab("Longitude")+ggtitle(genera[i])
print(p)
dev.off()}
 


