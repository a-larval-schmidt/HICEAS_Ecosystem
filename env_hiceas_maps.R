#temp, chla, SSH means over each leg
# libraries####
library(RColorBrewer)
library(ncdf4)
library(httr)
library(tidyverse)
library(ggOceanMaps)
library(marmap)
library(ggmap)
library(shape)
library(sf)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(gpkg)
library(TSP)
library(readxl)
library(patchwork)
library(stringr)
library(ggtext)
library(raster)
library(terra)
library(sp)
library(reshape2)

###base map#########
library(sf)
world<-ne_countries(scale="medium", returnclass = "sf")
load("A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Hawaii.RData")
oahu_df <- fortify(Hawaii)
not_sea<-oahu_df%>%
  filter(z>0)%>%
  dplyr::select(x,y)
island_map <- ggplot(data = world) +
  geom_raster(data = not_sea, aes(x = x, y = y)) +
  geom_sf() +
  coord_sf(xlim = c(206, 177), ylim = c(18, 32), label_graticule = "SW") + 
  theme_bw() +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")

PMNM<-read_sf("A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/PMNM_Shp_File/hi_noaa_nwhi_papahanaumokuakea.shp")
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
load("A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Hawaii.RData")
Hawaii_df <-fortify(Hawaii)

#August SST###############

#data extraction####
junk <- GET('https://oceanwatch.pifsc.noaa.gov/erddap/griddap/CRW_sst_v3_1_monthly.nc?sea_surface_temperature%5B(2023-08-01T12:00:00Z):1:(2023-08-31T12:00:00Z)%5D%5B(18):1:(32)%5D%5B(177):1:(206)%5D',
            write_disk("junk.nc", overwrite=TRUE))
nc <- nc_open('junk.nc')
v1 <- nc$var[[1]]
sst<- ncvar_get(nc,v1) 
lon <- v1$dim[[1]]$vals 
lat <- v1$dim[[2]]$vals 
nc_close(nc) 
rm(junk,v1)
file.remove('junk.nc')

df<-sst[,,1] 
rownames(df, do.NULL = TRUE, prefix = "row")
rownames(df) <- lon
colnames(df, do.NULL = TRUE, prefix = "col")
colnames(df) <-lat
dfreal<-as.data.frame(as.table(df))
dfreal<-dfreal%>%rename(x=Var1,y=Var2)
dfreal$x <- as.character(dfreal$x)
dfreal$x <- as.numeric(dfreal$x)
dfreal$y <- as.numeric(as.character(dfreal$y))
dfaug<-dfreal

#universal theme####
universal_theme <- theme_bw() +
  theme(text = element_text(size=10),
    axis.title = element_text(size=10),
    axis.text = element_text(size=9),
    legend.title = element_text(size=9), 
    legend.text = element_text(size=7), 
    strip.text.x = element_text(size=10),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.01, 0.01, 0.01, 0.1), "cm"),
    legend.position = "right",
    legend.justification = c(0, 0.5),  
    legend.margin = margin(l = -5),   
    legend.box.spacing = unit(0.2, "cm"))

###aug sst# map#####
saug <- ggplot() + 
  geom_tile(data=dfaug, aes(x=x, y=y, fill=Freq)) +
  scale_fill_gradientn(colors=c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", 
                                "yellow", "#FF7F00", "red", "#7F0000"),
                       na.value="gray90", limits = c(24, 29)) +
  guides(fill = "none") + # Hides legend
  universal_theme +       # Applies all formatting
  coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE) +
  geom_path(data = re_ordered_EEZ, aes(x = Lon, y = X2), color = "white", linewidth = 1) +
  scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1),
                     labels = function(x) { x - 360 }) +
  scale_y_continuous(breaks = round(seq(min(Hawaii_df$y), max(Hawaii_df$y), by = 2),1))+
  ylab(expression("Latitude (" * degree *"N)")) +
  xlab(expression("Longitude (" * degree * "W)"))

saug <- saug + geom_raster(data = not_sea, aes(x = x, y = y))

#October SST###############
#data extraction####
junk <- GET('https://oceanwatch.pifsc.noaa.gov/erddap/griddap/CRW_sst_v3_1_monthly.nc?sea_surface_temperature%5B(2023-10-01T12:00:00Z):1:(2023-11-01T12:00:00Z)%5D%5B(18):1:(32)%5D%5B(177):1:(206)%5D',
            write_disk("junk.nc", overwrite=TRUE))
nc <- nc_open('junk.nc')
v1 <- nc$var[[1]]
sst<- ncvar_get(nc,v1) 
lon <- v1$dim[[1]]$vals 
lat <- v1$dim[[2]]$vals 
nc_close(nc) 
rm(junk,v1)
file.remove('junk.nc')

df<-sst[,,1] 
rownames(df, do.NULL = TRUE, prefix = "row")
rownames(df) <- lon
colnames(df, do.NULL = TRUE, prefix = "col")
colnames(df) <-lat
dfreal<-as.data.frame(as.table(df))
dfreal<-dfreal%>%rename(x=Var1,y=Var2)
dfreal$x <- as.character(dfreal$x)
dfreal$x <- as.numeric(dfreal$x)
dfreal$y <- as.numeric(as.character(dfreal$y))

###oct sst# map#####
legend_title<-"SST (°C)"
soct <- ggplot() + 
  geom_tile(data=dfreal, aes(x=x, y=y, fill=Freq)) +
  scale_fill_gradientn(colors=c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F",
                                "yellow", "#FF7F00", "red", "#7F0000"),
                       na.value="gray90", limits = c(24, 29)) +
  guides(fill = guide_colorbar(title = legend_title,
                               barwidth = unit(0.3, "cm"), 
                               barheight = unit(1.2, "cm"))) +
  universal_theme + 
  coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE) +
  geom_path(data = re_ordered_EEZ, aes(x = Lon, y = X2), color = "white", linewidth = 1) +
  scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1),
                     labels = function(x) { x - 360 }) +
  scale_y_continuous(breaks = round(seq(min(Hawaii_df$y), max(Hawaii_df$y), by = 2),1))+
  ylab(expression("Latitude (" * degree *"N)")) +
  xlab(expression("Longitude (" * degree * "W)")) +
  labs(fill=legend_title)

soct <- soct + geom_raster(data = not_sea, aes(x = x, y = y))


#chla read in monthly 4km ocean color data (.nc file)#########
data = "A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/20230801-20230831_cmems_obs-oc_glo_bgc-plankton_myint_l4-olci-4km_P1M.nc"
chlor <-stack(data, varname = "CHL")
load("A:/My Drive/billfish_not_github/bathy_HI_30m_buffered_poly_1km.RData")
chlor_masked <- raster::mask(chlor, bathy_buffered_poly_1km, inverse = TRUE)
dfreal <- terra::as.data.frame(chlor_masked, xy = TRUE, na.rm = TRUE)
names(dfreal) <- c("x", "y", "Freq")
dfreal<-dfreal%>%
  filter(y>17)%>%
  filter(y<33)%>%
  mutate(x=ifelse(x<0, x+360, x))
ca<-dfreal

#october chla#
path <- ("A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/20231001-20231031_cmems_obs-oc_glo_bgc-plankton_myint_l4-olci-4km_P1M.nc")
chlor <-stack(path, varname = "CHL")
load("A:/My Drive/billfish_not_github/bathy_HI_30m_buffered_poly_1km.RData")
chlor_masked <- raster::mask(chlor, bathy_buffered_poly_1km, inverse = TRUE)
dfreal <- terra::as.data.frame(chlor_masked, xy = TRUE, na.rm = TRUE)
names(dfreal) <- c("x", "y", "Freq")
dfreal<-dfreal%>%
  filter(y>17)%>%
  filter(y<33)%>%
  mutate(x=ifelse(x<0, x+360, x))

#map AUGUST CHLOROPHYLL (No Legend)###############
chlaug <- ggplot() + 
  geom_tile(data=ca, aes(x=x, y=y, fill=log(Freq))) +
  scale_fill_distiller(palette = "Greens", direction = 1, na.value = "gray70") +
  guides(fill = "none") + # Hides legend
  universal_theme + 
  coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE) +
  geom_path(data = re_ordered_EEZ, aes(x = Lon, y = X2), color = "white", linewidth = 1) +
  scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1),
                     labels = function(x) { x - 360 }) +
  ylab(expression("Latitude (" * degree *"N)")) +
  xlab(expression("Longitude (" * degree * "W)")) +
  labs(fill=legend_title)

chlaug <- chlaug + geom_raster(data = not_sea, aes(x = x, y = y))

#####map OCTOBER CHLOROPHYLL (With Legend)############
legend_title <- expression("Log\nChlorophyll"*" (mg" ~ m^{-3} ~ ")")

chloct <- ggplot() + 
  geom_tile(data=dfreal, aes(x=x, y=y, fill=log(Freq))) +
  scale_fill_distiller(palette = "Greens", direction = 1, na.value = "gray70") +
  guides(fill = guide_colorbar(title = legend_title,
                               barwidth = unit(0.3, "cm"), 
                               barheight = unit(1.2, "cm"))) +
  universal_theme +
  coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE) +
  geom_path(data = re_ordered_EEZ, aes(x = Lon, y = X2), color = "white", linewidth = 1) +
  scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1),
                     labels = function(x) { x - 360 }) +
  ylab(expression("Latitude (" * degree *"N)")) +
  xlab(expression("Longitude (" * degree * "W)"))

chloct <- chloct + geom_raster(data = not_sea, aes(x = x, y = y))


#all 4 layout##########
ggsave("sst_aug_uniform.png", plot = saug , width = 5.5, height = 3.5)
ggsave("sst_oct_uniform.png", plot = soct, width = 6, height = 4)
ggsave("chla_aug_uniform.png", plot = chlaug, width = 5.5, height = 3.5)
ggsave("chla_oct_uniform.png", plot = chloct, width = 6, height = 4)

# --- FINAL COMBINED PLOT ---
combined_plots <- (saug + soct) / (chlaug + chloct) + 
  plot_layout(heights = c(1, 1)) & 
  theme(plot.margin = unit(c(0.1, 0.1, -0.2, 0.5), "cm"), 
        plot.tag = element_text(size = 12, face = "bold", 
                                margin = margin(r = 10, b = 5))) # Added margin(r=5, b=5) to nudge letters away from axis


final_plot <- combined_plots +
  plot_annotation(tag_levels = 'A')

ggsave(filename = "fig2_final_combined_plot.png",
       path="A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Map Figures/",
       plot = final_plot,
       width = 8,
       height = 3.8, 
       dpi = 300)



#CRP vs PRP's coordinates############
load("A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Hawaii.RData")
Hawaii_df <-marmap::fortify.bathy(Hawaii)
library(sf)
PMNM<-read_sf("A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/PMNM_Shp_File/hi_noaa_nwhi_papahanaumokuakea.shp")
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
jsky<-read.csv("A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/HICEAS_IKMT_Locations_Times_KY_JS.xlsx - Sheet1 (2).csv")
cruz.list.save$ndas.data[[2]]
lat_dd<-cruz.list.save$ndas.data[[2]]$y
lon_dd<-cruz.list.save$ndas.data[[2]]$x
coord2<-as_tibble(cbind(lat_dd,lon_dd))
ggplot() +geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4",limits=c(-10000,0),
                      na.value=(color="#003300"),guide="none")+
  #geom_point(data=ik,mapping=aes(y=lat_dd, x=ifelse(lon_dd<0, lon_dd+360, lon_dd)), color="blue", alpha=0.5)+
  geom_label(data=ik,mapping=aes(y=lat_dd, x=ifelse(lon_dd<0, lon_dd+360, lon_dd),
                                 label=Station), color="blue4", alpha=0.5)+
  #geom_point(data=jsky,mapping=aes(y=Lat_DD, x=ifelse(Lon_DD<0, Lon_DD+360,Lon_DD)), color="red", alpha=0.5)+
  #geom_point(data=coord2,mapping=aes(y=lat_dd, x=lon_dd),color="red", alpha=0.2, size=3) +
  #coord_sf(xlim=c(191.5,192.5), ylim=c(21.5, 22.5)) + for zooming in
  coord_sf(xlim=c(178.4,204.9), ylim=c(18.6, 31.4)) +
  scale_x_continuous(breaks = round(seq(min(Hawaii_df$x, na.rm=T), max(Hawaii_df$x,na.rm=T), by = 5),1),labels = function(x) { x - 360 })+
  theme_bw()+theme(axis.title = element_text(size=11),axis.text=element_text(size=10),
                   axis.ticks.length = unit(0.25, "cm"),
                   plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                   panel.grid = element_blank(),text = element_text(size=10))+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  geom_label(data=jsky,mapping=aes(y=Lat_DD,
                                   x=ifelse(Lon_DD<0, Lon_DD+360,Lon_DD), 
                                   label=Station), color="red2", alpha=0.5)+
  geom_label(data=jsky,mapping=aes(y=corrected_latdd_for_plotting,
                                   x=ifelse(corrected_londd_for_plotting<0,
                                            corrected_londd_for_plotting+360,
                                            corrected_londd_for_plotting), 
                                   label=Station), color="green4", alpha=0.5)+
  scale_color_manual(
    name = "Data Source",
    values = c(
      "original" = "blue4",     # Bluish Green
      "CRPs" = "red2",         # Orange
      "Transcribed" = "green4"   # Reddish Purple
    )
  ) 
 

  
  
  
####with less horrible colors######  
ggplot() +geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4",limits=c(-10000,0),
                      na.value=(color="#003300"),guide="none")+
  #geom_point(data=ik,mapping=aes(y=lat_dd, x=ifelse(lon_dd<0, lon_dd+360, lon_dd)), color="blue", alpha=0.5)+
  geom_label(data=ik,mapping=aes(y=lat_dd, x=ifelse(lon_dd<0, lon_dd+360, lon_dd),
                                 label=Station), color="#009E73", alpha=0.5)+
  #geom_point(data=jsky,mapping=aes(y=Lat_DD, x=ifelse(Lon_DD<0, Lon_DD+360,Lon_DD)), color="red", alpha=0.5)+
  #geom_point(data=coord2,mapping=aes(y=lat_dd, x=lon_dd),color="red", alpha=0.2, size=3) +
  #coord_sf(xlim=c(191.5,192.5), ylim=c(21.5, 22.5)) + for zooming in
  coord_sf(xlim=c(178.4,204.9), ylim=c(18.6, 31.4)) +
  scale_x_continuous(breaks = round(seq(min(Hawaii_df$x, na.rm=T), max(Hawaii_df$x,na.rm=T), by = 5),1),labels = function(x) { x - 360 })+
  theme_bw()+theme(axis.title = element_text(size=11),axis.text=element_text(size=10),
                   axis.ticks.length = unit(0.25, "cm"),
                   plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                   panel.grid = element_blank(),text = element_text(size=10),
                   legend.position = "inside",legend.position.inside=c(0.22,0.03),
                   legend.justification = c(1, 0),legend.background = element_rect(fill = "transparent"),
                   legend.key = element_rect(fill = "transparent"))+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  geom_label(data=jsky,mapping=aes(y=Lat_DD,
                                   x=ifelse(Lon_DD<0, Lon_DD+360,Lon_DD), 
                                   label=Station), color="#E69F00", alpha=0.5)+
  geom_label(data=jsky,mapping=aes(y=corrected_latdd_for_plotting,
                                   x=ifelse(corrected_londd_for_plotting<0,
                                            corrected_londd_for_plotting+360,
                                            corrected_londd_for_plotting), 
                                   label=Station), color="#CC79A7", alpha=0.5)+

  scale_color_manual(
    name = "Data Source",
    values = c(
      "original" = "#009E73",     # Bluish Green
      "CRPs" = "#E69F00",         # Orange
      "Transcribed" = "#CC79A7"   # Reddish Purple
    )
  )
####gemini fix#################
library(ggplot2)
library(sf)

# Defined to make the code cleaner
# Adjust xlim upper bound from 204.9 to 208 to catch points at -154 longitude (206)
map_xlim <- c(178.4, 208) 
map_ylim <- c(18.6, 31.4)

final_map <- ggplot() +
  # 1. The Base Map
  geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradient(
    high = "lightskyblue1", low = "lightsteelblue4",
    limits = c(-10000, 0), na.value = "#003300", guide = "none"
  ) +
  # 5. Coordinate System & Theme
  coord_sf(xlim = map_xlim, ylim = map_ylim) +
  scale_x_continuous(
    breaks = round(seq(min(Hawaii_df$x, na.rm = T), max(Hawaii_df$x, na.rm = T), by = 5), 1),
    labels = function(x) { x - 360 }
  ) +
  geom_path(data = re_ordered_EEZ, aes(x = Lon, y = X2), color = "white", linewidth = 1) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    axis.ticks.length = unit(0.25, "cm"),
    plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
    panel.grid = element_blank(),
    text = element_text(size = 10),
    legend.position = "inside",
    legend.position.inside = c(0.22, 0.03),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )+
  # 3. RED Labels (from 'jsky' dataset - standard lat/lon)
  # Note: str(jsky) shows uppercase 'Lat_DD' and 'Lon_DD'
  geom_label(
    data = jsky, 
    mapping = aes(y = Lat_DD, 
                  x = ifelse(Lon_DD < 0, Lon_DD + 360, Lon_DD), 
                  label = Station), 
    color = "red", alpha = 0.1
  ) +
  # 2. BLUE Labels (from 'ik' dataset)
  # Note: str(ik) shows lowercase 'lat_dd' and 'lon_dd'
  geom_label(
    data = ik, 
    mapping = aes(y = lat_dd, 
                  x = ifelse(lon_dd < 0, lon_dd + 360, lon_dd), 
                  label = Station), 
    color = "blue", alpha = 0.1
  ) +
  
  # 4. PURPLE Labels (from 'jsky' dataset - corrected lat/lon)
  # We subset data!=0 inline here so we don't plot points at the equator
  geom_label(
    data = subset(jsky, corrected_latdd_for_plotting != 0), 
    mapping = aes(y = corrected_latdd_for_plotting, 
                  x = ifelse(corrected_londd_for_plotting < 0, 
                             corrected_londd_for_plotting + 360, 
                             corrected_londd_for_plotting), 
                  label = Station), 
    color = "green", alpha = 0.1
  ) 
  


# Explicitly print the map to see it
print(final_map)


kym<- ggplot() +
  # 1. The Base Map
  geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradient(
    high = "lightskyblue1", low = "lightsteelblue4",
    limits = c(-10000, 0), na.value = "#003300", guide = "none"
  ) +
  # 5. Coordinate System & Theme
  coord_sf(xlim = map_xlim, ylim = map_ylim) +
  scale_x_continuous(
    breaks = round(seq(min(Hawaii_df$x, na.rm = T), max(Hawaii_df$x, na.rm = T), by = 5), 1),
    labels = function(x) { x - 360 }
  ) +  labs(title="CRP/Justin Points")+
  geom_path(data = re_ordered_EEZ, aes(x = Lon, y = X2), color = "white", linewidth = 1) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    axis.ticks.length = unit(0.25, "cm"),
    plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
    panel.grid = element_blank(),
    text = element_text(size = 10),
    legend.position = "inside",
    legend.position.inside = c(0.22, 0.03),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )+geom_label(
    data = jsky, 
    mapping = aes(y = Lat_DD, 
                  x = ifelse(Lon_DD < 0, Lon_DD + 360, Lon_DD), 
                  label = Station), 
    color = "red", alpha = 0.1
  ) 


dre<- ggplot() +
  # 1. The Base Map
  geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradient(
    high = "lightskyblue1", low = "lightsteelblue4",
    limits = c(-10000, 0), na.value = "#003300", guide = "none"
  ) +
  # 5. Coordinate System & Theme
  coord_sf(xlim = map_xlim, ylim = map_ylim) +
  scale_x_continuous(
    breaks = round(seq(min(Hawaii_df$x, na.rm = T), max(Hawaii_df$x, na.rm = T), by = 5), 1),
    labels = function(x) { x - 360 }
  ) +
  geom_path(data = re_ordered_EEZ, aes(x = Lon, y = X2), color = "white", linewidth = 1) +
  theme_bw() +  labs(title="'original' Points with - for joining lon vlaues")+
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    axis.ticks.length = unit(0.25, "cm"),
    plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
    panel.grid = element_blank(),
    text = element_text(size = 10),
    legend.position = "inside",
    legend.position.inside = c(0.22, 0.03),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"))+
    geom_label(
      data = ik, 
      mapping = aes(y = lat_dd, 
                    x = ifelse(lon_dd < 0, lon_dd + 360, lon_dd), 
                    label = Station), 
      color = "blue", alpha = 0.1
    ) 
corrected<- ggplot() +
  # 1. The Base Map
  geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradient(
    high = "lightskyblue1", low = "lightsteelblue4",
    limits = c(-10000, 0), na.value = "#003300", guide = "none"
  ) +
  # 5. Coordinate System & Theme
  coord_sf(xlim = map_xlim, ylim = map_ylim) +
  scale_x_continuous(
    breaks = round(seq(min(Hawaii_df$x, na.rm = T), max(Hawaii_df$x, na.rm = T), by = 5), 1),
    labels = function(x) { x - 360 }
  ) +
  geom_path(data = re_ordered_EEZ, aes(x = Lon, y = X2), color = "white", linewidth = 1) +
  theme_bw() +
  labs(title=" not actually Corrected Points")+
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    axis.ticks.length = unit(0.25, "cm"),
    plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
    panel.grid = element_blank(),
    text = element_text(size = 10),
    legend.position = "inside",
    legend.position.inside = c(0.22, 0.03),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"))+
  geom_label(
    data = subset(jsky, corrected_latdd_for_plotting != 0), 
    mapping = aes(y = corrected_latdd_for_plotting, 
                  x = ifelse(corrected_londd_for_plotting < 0, 
                             corrected_londd_for_plotting + 360, 
                             corrected_londd_for_plotting), 
                  label = Station), 
    color = "green", alpha = 0.1
  ) 
library(patchwork)

kym/dre/corrected

