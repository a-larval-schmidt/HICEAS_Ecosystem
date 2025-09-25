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
ik<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/SE2303 ID Log - ikmt_metadata.csv")
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
inventory<-read_excel("Data/SE2303 ID Log.xlsx", sheet="Inventory and Log")
ik2<-left_join(ik, inventory, by="Station")
#add in eDNA sites so far for funsies#
edan<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/eDNA Log - Cast_info (1).csv")

#FISH POINTS########

#QCing tows and making a smaller, easier to work with data frame
ik_lite<-ik2%>%
  mutate(lon_dd=ifelse(lon_dd < 0,lon_dd+ 360, lon_dd))%>%
  mutate(depth_m=as.numeric(Max.Depth))%>%
  mutate(mean_flow=((Flowmeter.1.End-Flowmeter.1.Start)+
                      (Flowmeter.2.End-Flowmeter.2.Start))/2,
         na.rm=T)%>%
  filter(depth_m>30)%>%
  mutate(month=month(date.x))%>%
  mutate(dec_date=julian(date.x))%>%
  mutate(month = factor(month, levels = c(7, 8, 10,11), labels = c("July", "August", "October","November"))) %>%
  dplyr::select(Station,lat_dd,lon_dd,mean_flow, month, date.x, dec_date)
#structuring fish ids
SE2303_ID<-read_excel("Data/SE2303 ID Log.xlsx", sheet = "ID&Length")
SE2303_ID$total<-as.numeric(SE2303_ID$total)
SE2303_ID$`1_over_split_fraction`<-as.numeric(SE2303_ID$`1_over_split_fraction`)

#unite dfs
se23<-left_join(SE2303_ID, ik_lite, join_by("station"=="Station"))

#density
se23<-se23%>%
  filter(Trawl=="6' IKMT")%>%
  filter(is.na(month)==F)%>%
  mutate(density=(total/mean_flow))

#old base map data####
world<-ne_countries(scale="medium", returnclass = "sf")
oahu_raster <- raster(file.path("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/hi_eez_extract_grid.tiff"))
#from here: https://www.ncei.noaa.gov/maps/bathymetry/
#settings ETPO_2022(Bedrock, 15 arc seconds)
df <- fortify(as.bathy(oahu_raster))
#convert to 360
oahu_df<-df%>%mutate(x=ifelse(df$x < 0, df$x + 360, df$x))
load("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Hawaii.RData")
#oahu_df <- fortify(Hawaii)
Hawaii_df <- fortify(Hawaii)
#oahu_map <- ggplot(data = world) +geom_raster(data = oahu_df, aes(x = x, y = y, fill = z)) +labs(fill = "Depth (m)") +scale_fill_gradient(high = "lightskyblue1", low = "cornflowerblue", limits = c(-10000, 1000)) +geom_sf() +coord_sf(xlim = c(206, 177), ylim = c(18, 32), label_graticule = "SW") + theme_bw() +theme(axis.ticks.length = unit(0.25, "cm")) +ylab("Latitude") +xlab("Longitude")
#http://127.0.0.1:19395/graphics/plot_zoom_png?width=1920&height=1027
#add in monument boundary######
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

#jessie's Base map rendering/colors######
#first download the bathymetry data for Hawaiʻi
getNOAA.bathy(lon1 = 177,lon2 = -154,lat1 = 18,lat2 =32,resolution = 1, antimeridian = TRUE) -> Hawaii
Hawaii_df <- fortify(Hawaii)

# make similar to Marianas 2025 map #EEE
base_map<-ggplot() +
  geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +labs(fill = "Depth (m)")+
  scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4",limits=c(-10000,0), na.value="black",guide="none")+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  geom_point(data=se23,mapping=aes(y=lat_dd, x=lon_dd, color=factor(month)), size=2.2) +
  scale_color_manual(values = c("chartreuse4","#800074","#c99b38","orange4"),labels = c("July","Aug","Oct","Nov"),name="") +
  #geom_label(data=se23,mapping=aes(y=lat_dd, x=lon_dd, label=station))+
  coord_sf(xlim=c(178.3,204.3), ylim=c(18.6, 31.4)) +
  scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1)) +
  theme_bw()+theme(axis.ticks.length = unit(0.25, "cm"),panel.grid = element_blank(),text = element_text(size=16))+
  labs(y="Latitude", x="Longitude")#,title= "PIRIS 2023:\nHawaiʻi Pae ʻĀina")

ggsave("HICEAS_Distribution_Map.png", plot = base_map, width = 6, height = 5,path = "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Map Figures/")
#justins base quantiles PLUS gemini#####
# Separate presence and absence data
df_pres <- df %>% filter(total > 0)
df_abs <- anti_df %>% filter(total == 0)

# Check if there is presence data to calculate quantiles
if (nrow(df_pres) > 0) {
  # Calculate quantiles for total counts
  quant_levels <- quantile(df_pres$total, probs = seq(0, 1, by = 0.25), na.rm = TRUE)
  
  # Add a new column to the presence data for point sizes based on quantiles
  df_pres$Number <- 1
  df_pres$Number[df_pres$total > 0 & df_pres$total < quant_levels[2]] <- 1.5
  df_pres$Number[df_pres$total >= quant_levels[2] & df_pres$total < quant_levels[3]] <- 2
  df_pres$Number[df_pres$total >= quant_levels[3] & df_pres$total < quant_levels[4]] <- 2.5
  df_pres$Number[df_pres$total >= quant_levels[4]] <- 3
  
  # Define the legend labels and sizes for both absence and presence
  legend_labels <- c("0",
                     paste0("<", quant_levels[2]),
                     paste0(quant_levels[2], "-", quant_levels[3]),
                     paste0(quant_levels[3], "-", quant_levels[4]),
                     paste0(">", quant_levels[4]))
  
  cex_size_list <- c(1, 1.5, 2, 2.5, 3)
  pch_list <- c(1, 19, 19, 19, 19)
  
} else {
  # If no presence data, all points are absence points.
  df_pres <- data.frame() # Create an empty data frame to avoid errors
  df_abs <- se23 %>% filter(total == 0) # All points are absence points
  legend_labels <- "0"
  cex_size_list <- 1
  pch_list <- 1
}

# Now, create the plot using Base R
# Now, create the plot using Base R
png(paste0("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/", fish_taxa, "_map.png"), height=5, width=8, units="in", res=300)

par(mar = c(5, 5, 4, 2) + 0.1, cex = 1.5, bg = "white") 
# Sets margins, text size, and plot background

# This is the Base R equivalent of geom_raster and scale_fill_gradient
# `xlim` and `ylim` are set here.
plot(Hawaii,
     image = TRUE,
     land = TRUE,
     lwd = 0.1,
     # Corrected bpal to ensure colors are applied to appropriate depth ranges
     bpal = list(c(0, max(Hawaii), "lightskyblue1"), # For depths >= 0 (surface water)
                 c(min(Hawaii), 0, "lightsteelblue4")), # For depths < 0 (deeper water)
     xlim = c(178.3, 204.6),
     ylim = c(18.6, 31.4),
     xaxt = 'n', 
     yaxt = 'n',
     xlab = "Longitude", 
     ylab = "Latitude",
     main = paste(fish_taxa, "Density Locations"))

# Add the white EEZ boundary line
graphics::lines(re_ordered_EEZ$Lon, re_ordered_EEZ$X2, lwd=3, col="white")

# Plot the absence points
points(df_abs$lon_dd, df_abs$lat_dd, pch = 1, cex = 2.2, lwd = 1, col = "black")

# Plot the presence points (only if there are any)
if (nrow(df_pres) > 0) {
  points(df_pres$lon_dd, df_pres$lat_dd, pch = 19, cex = df_pres$Number, col = "black")
}

# Add axes and tick marks
axis(1, at = round(seq(180, 205, by = 5), 1))
axis(2, at = round(seq(18, 32, by = 2), 1))

# Add the legend
# Define the legend labels and sizes for both absence and presence
if (nrow(df_pres) > 0) {
  legend_labels <- c("0",
                     paste0("<", quant_levels[2]),
                     paste0(quant_levels[2], "-", quant_levels[3]),
                     paste0(quant_levels[3], "-", quant_levels[4]),
                     paste0(">", quant_levels[4]))
  
  cex_size_list <- c(2.2, 1.5, 2, 2.5, 3)
  pch_list <- c(1, 19, 19, 19, 19)
  
} else {
  legend_labels <- "0"
  cex_size_list <- 2.2
  pch_list <- 1
}

legend("topright",
       legend = legend_labels,
       pch = pch_list,
       pt.cex = cex_size_list,
       title = "Density (Counts)")

dev.off()
#Make sure to handle the case where anti_df might be empty
# map_plot <- ggplot() +
#   geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +
#   labs(fill = "Depth (m)") +
#   scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4", limits = c(-10000, 0), na.value = "black", guide = "none") +
#   geom_point(data = df, aes(y = lat_dd, x = lon_dd, size = density), shape = 19, color="black") +
#   scale_size_continuous(range = c(1.5, 5)) +
#   geom_point(data = anti_df, aes(y = lat_dd, x = lon_dd), shape = 1, size = 2.2) +
#   coord_sf(xlim = c(178.3, 204.6), ylim = c(18.6, 31.4)) +
#   scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5), 1)) +
#   theme_bw() +
#   theme(
#     axis.ticks.length = unit(0.25, "cm"),
#     panel.grid = element_blank(),
#     text = element_text(size = 16)
#   ) +
#   ylab("Latitude") +
#   xlab("Longitude") +
#   labs(title = paste(fish_taxa, "Density Locations")) +
#   geom_path(data = re_ordered_EEZ, aes(x = Lon, y = X2), color = "white", linewidth = 1)

#ggsave(paste0(fish_taxa, "_map.png"), plot = map_plot, path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/",width = 8, height = 6)




#plot fish####

#uku##############

fish_taxa="Aprion"
df<-paste("se23_",fish_taxa)
df<-se23%>%
  filter(genus==fish_taxa)
##plots
metric=density
yaxistitle="Larval density (larva/m^3 seawater)"
ggplot()+geom_point(data=df, mapping=aes(x=datey, y=density))+
  labs(title =paste(fish_taxa), x="Date", y=yaxistitle)+
  geom_smooth(method="loess")

anti_df<-se23%>%filter(genus!=fish_taxa)
merp+geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density, color=(month)), shape=19)+
  geom_point(data=anti_df,mapping=aes(y=lat_dd, x=lon_dd), shape=1)+labs(title =paste(fish_taxa,metric))
#snappers generally######
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
