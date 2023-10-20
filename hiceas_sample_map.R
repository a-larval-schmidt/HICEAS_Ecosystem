#HICEAS sampling map
#map libraries#####
###=---Cam's map---
library(sp) #already in raster
library(raster)# error
library(ggplot2)
library(scales)
library(rgdal)
library(marmap)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(viridis)
library(ggplot2)
library(raster)
library(ggnewscale)
world<-ne_countries(scale="medium", returnclass = "sf")

#ik data libraries####
library(tidyverse)
ik<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/IKMT Log - Sheet1.csv")
#convert ship coordinates into DMS then to DD
ik<-ik %>% mutate(lat_deg= as.numeric(substr(Lat.In,1,2)))
ik<-ik %>% mutate(lon_deg= as.numeric(substr(Lon.In,1,4)))
ik<-ik %>% mutate(lat_minsec= as.numeric(substr(Lat.In,3,10)))
ik<-ik %>% mutate(long_minsec= as.numeric(substr(Lon.In,5,10)))
ik<-ik %>% mutate(lat_dd= lat_deg+(lat_minsec/60))
ik<-ik %>% mutate(lon_dd= lon_deg+(long_minsec/60))

#cam's mapping codefor base map####
oahu_raster <- raster(file.path("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/hi_eez_extract_grid.tiff"))
oahu_df <- fortify(as.bathy(oahu_raster))
str(oahu_df)
oahu_map <- ggplot(data=world) +
  geom_raster(data = oahu_df, aes(x = x, y = y, fill = z)) +labs(fill = "Depth (m)")+
  scale_fill_gradient(high = "lightskyblue1", low = "cornflowerblue",limits=c(-1000,1000))+
  new_scale_fill()+
  scale_fill_continuous(labels = scales::label_number(scale = 1000, suffix = "k"))+
  theme_bw()+geom_sf()+coord_sf(xlim=c(-154.9,-179.9999), ylim=c(19, 31.7))


oahu_map
#map with points####
oahu_map+geom_point(data=ik, mapping=aes(y=lat_dd, x=lon_dd))
                    
                    
                    
  scale_color_viridis_d(option="B")+scale_shape_manual(values=shapes)+
  theme(axis.title = element_blank(),legend.position = "bottom",legend.box = "horizonal",
        axis.ticks = element_blank(),axis.text = element_blank(),
        legend.title=element_text(size=16),legend.text=element_text(size=20))+
  labs(color="Sampling Effort by Cruise Identifier", shape="Larva present (1) or absent(0)", size=element_blank())+
  guides(color = guide_legend(override.aes = list(size = 5), order=1),
         shape=guide_legend(override.aes = list(size = 5), order=2),
         fill= guide_colorbar(barwidth = 20, barheight = 10, order=3),
         size="none")
