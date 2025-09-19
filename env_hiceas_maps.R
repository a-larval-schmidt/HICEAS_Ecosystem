#temp, chla, SSH means over each leg
# libraries####
library(ncdf4)
library(httr)
library(tidyverse)
library(ggOceanMaps)
library(marmap) #masks as.raster() from package:grDevices
library(ggmap)#masks qmap from ggOceanMaps
library(shape)
library(sf)#Linking to GEOS 3.13.1, GDAL 3.11.0, PROJ 9.6.0; sf_use_s2() is TRUE
###base map#########
library(sf)
world<-ne_countries(scale="medium", returnclass = "sf")
load("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Hawaii.RData")
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
island_map
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
island_map+geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)

#August SST###############
#data extraction####
junk <- GET('https://oceanwatch.pifsc.noaa.gov/erddap/griddap/CRW_sst_v3_1_monthly.nc?sea_surface_temperature%5B(2023-08-01T12:00:00Z):1:(2023-08-31T12:00:00Z)%5D%5B(18):1:(32)%5D%5B(177):1:(206)%5D',
            write_disk("junk.nc", overwrite=TRUE))
nc <- nc_open('junk.nc')
names(nc$var)
v1 <- nc$var[[1]]
sst<- ncvar_get(nc,v1) #Extract analysed_sst, reads data from the netCDF file, only works if you have already opened the file, shows as a multi-dimensional array
dim(sst)
dates <- as.POSIXlt(v1$dim[[3]]$vals,origin='1970-01-01',tz='GMT') #get the dates for each time step
lon <- v1$dim[[1]]$vals #gives vector of longitude
#lon=lon-360
lat <- v1$dim[[2]]$vals #gives vector of latitude
nc_close(nc) #this step is important, otherwise you risk data loss
rm(junk,v1)
file.remove('junk.nc')

df<-sst[,,1] #because it is already month mean for august no need to take mean of the whole month
rownames(df, do.NULL = TRUE, prefix = "row")
rownames(df) <- lon
colnames(df, do.NULL = TRUE, prefix = "col")
colnames(df) <-lat
dfreal<-as.data.frame(as.table(df))
dfreal<-dfreal%>%rename(y=Var1,x=Var2)


###sst# map#####
sst_plot<-ggplot()+ geom_tile(data=dfreal,aes(x=x,y=y, fill=Freq))+
  scale_fill_gradientn(colors=c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", 
                                "yellow", "#FF7F00", "red", "#7F0000"),
                       na.value="gray90",limits = c(24, 32))+
  theme(legend.title=element_text(size=20),legend.text=element_text(size=20),
        axis.text = element_text(size=10), axis.title = element_text(size=20),
        legend.direction = "vertical", legend.box = "vertical",strip.text.x=element_text(size=20))+
  labs(fill="SST (°C)")+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")+labs(title ="Sampling Stations and Mean SST August 2023")
sst_plot+geom_raster(data = not_sea, aes(x = x, y = y))
#map august points##
#load points form hiceas_sample_map.R
august<-se23%>%
  filter(month==8)
saug<-sst_plot+geom_point(data=august, mapping=aes(x=lon_dd, y=lat_dd))+geom_raster(data = not_sea, aes(x = x, y = y))
#ggsave(filename="Sampling Stations and Mean SST August 2023.png",plot=saug,height=5, width=8, units="in", dpi=300)
#October SST###############
#data extraction####
junk <- GET('https://oceanwatch.pifsc.noaa.gov/erddap/griddap/CRW_sst_v3_1_monthly.nc?sea_surface_temperature%5B(2023-10-01T12:00:00Z):1:(2023-11-01T12:00:00Z)%5D%5B(18):1:(32)%5D%5B(177):1:(206)%5D',
            write_disk("junk.nc", overwrite=TRUE))
nc <- nc_open('junk.nc')
names(nc$var)
v1 <- nc$var[[1]]
sst<- ncvar_get(nc,v1) #Extract analysed_sst, reads data from the netCDF file, only works if you have already opened the file, shows as a multi-dimensional array
dim(sst)
dates <- as.POSIXlt(v1$dim[[3]]$vals,origin='1970-01-01',tz='GMT') #get the dates for each time step
lon <- v1$dim[[1]]$vals #gives vector of longitude
#lon=lon-360
lat <- v1$dim[[2]]$vals #gives vector of latitude
nc_close(nc) #this step is important, otherwise you risk data loss
rm(junk,v1)
file.remove('junk.nc')

df<-sst[,,1] #because it is already month mean for august no need to take mean of the whole month
rownames(df, do.NULL = TRUE, prefix = "row")
rownames(df) <- lon
colnames(df, do.NULL = TRUE, prefix = "col")
colnames(df) <-lat
dfreal<-as.data.frame(as.table(df))
dfreal<-dfreal%>%rename(y=Var1,x=Var2)

###sst# map#####
sst_plot<-ggplot()+ geom_tile(data=dfreal,aes(x=x,y=y, fill=Freq))+
  scale_fill_gradientn(colors=c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F",
                                "yellow", "#FF7F00", "red", "#7F0000"),
                       na.value="gray90",limits = c(24, 32))+
  theme(legend.title=element_text(size=20),legend.text=element_text(size=20),
        axis.text = element_text(size=10), axis.title = element_text(size=20),
        legend.direction = "vertical", legend.box = "vertical",strip.text.x=element_text(size=20))+
  labs(fill="SST (°C)")+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")+labs(title ="Sampling Stations and Mean SST October 2023")
sst_plot+geom_raster(data = not_sea, aes(x = x, y = y))
october<-se23%>%
  filter(month==10)
soct<-sst_plot+geom_point(data=october, mapping=aes(x=lon_dd, y=lat_dd))+geom_raster(data = not_sea, aes(x = x, y = y))
ggsave(filename="Sampling Stations and Mean SST October 2023.png",plot=soct,height=5, width=8, units="in", dpi=300)
#bind sst######
library(patchwork)
saug+soct+plot_annotation(tag_levels=c("A","B"))


#August chla###########
#data extraction####
junk <- GET('https://oceanwatch.pifsc.noaa.gov/erddap/griddap/noaa_snpp_chla_monthly.nc?chlor_a%5B(2023-08-01T12:00:00Z):1:(2023-08-31T12:00:00Z)%5D%5B(18):1:(32)%5D%5B(177):1:(206)%5D',
            write_disk("junk.nc", overwrite=TRUE))
nc <- nc_open('junk.nc')
names(nc$var)
v1 <- nc$var[[1]]
sst<- ncvar_get(nc,v1) #Extract analysed_sst, reads data from the netCDF file, only works if you have already opened the file, shows as a multi-dimensional array
dim(sst)
dates <- as.POSIXlt(v1$dim[[3]]$vals,origin='1970-01-01',tz='GMT') #get the dates for each time step
lon <- v1$dim[[1]]$vals #gives vector of longitude
#lon=lon-360
lat <- v1$dim[[2]]$vals #gives vector of latitude
nc_close(nc) #this step is important, otherwise you risk data loss
rm(junk,v1)
file.remove('junk.nc')

df<-sst[,,1] #because it is already month mean for august no need to take mean of the whole month
rownames(df, do.NULL = TRUE, prefix = "row")
rownames(df) <- lon
colnames(df, do.NULL = TRUE, prefix = "col")
colnames(df) <-lat
dfreal<-as.data.frame(as.table(df))
dfreal<-dfreal%>%rename(y=Var1,x=Var2)
#mapping####
chl_map<-ggplot()+ geom_tile(data=dfreal,aes(x=x,y=y, fill=Freq))+
  scale_fill_gradientn(colors=c("darkseagreen1","aquamarine4"),na.value="gray 90",limits = c(0, 1))+
  theme(legend.title=element_text(size=20),legend.text=element_text(size=20),
        axis.text = element_text(size=10), axis.title = element_text(size=20),
        legend.direction = "vertical", legend.box = "vertical",strip.text.x=element_text(size=20))+
  labs(fill="Chlorophyll Concentration, mg/m-3")+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")+labs(title ="Sampling Stations and Mean Chla August 2023")
  
chlaug<-chl_map+geom_point(data=august, mapping=aes(x=lon_dd, y=lat_dd))+geom_raster(data = not_sea, aes(x = x, y = y))
#october chla###########
#data extraction####
junk <- GET('https://oceanwatch.pifsc.noaa.gov/erddap/griddap/noaa_snpp_chla_monthly.nc?chlor_a%5B(2023-10-01T12:00:00Z):1:(2023-10-31T12:00:00Z)%5D%5B(18):1:(32)%5D%5B(177):1:(206)%5D',
            write_disk("junk.nc", overwrite=TRUE))
nc <- nc_open('junk.nc')
names(nc$var)
v1 <- nc$var[[1]]
sst<- ncvar_get(nc,v1) #Extract analysed_sst, reads data from the netCDF file, only works if you have already opened the file, shows as a multi-dimensional array
dim(sst)
dates <- as.POSIXlt(v1$dim[[3]]$vals,origin='1970-01-01',tz='GMT') #get the dates for each time step
lon <- v1$dim[[1]]$vals #gives vector of longitude
#lon=lon-360
lat <- v1$dim[[2]]$vals #gives vector of latitude
nc_close(nc) #this step is important, otherwise you risk data loss
rm(junk,v1)
file.remove('junk.nc')

df<-sst[,,1] #because it is already month mean for august no need to take mean of the whole month
rownames(df, do.NULL = TRUE, prefix = "row")
rownames(df) <- lon
colnames(df, do.NULL = TRUE, prefix = "col")
colnames(df) <-lat
dfreal<-as.data.frame(as.table(df))
dfreal<-dfreal%>%rename(y=Var1,x=Var2)
#mapping####
chl_map<-ggplot()+ geom_tile(data=dfreal,aes(x=x,y=y, fill=Freq))+
  scale_fill_gradientn(colors=c("darkseagreen1","aquamarine4"),na.value="gray 90", limits = c(0, 1))+
  theme(legend.title=element_text(size=20),legend.text=element_text(size=20),
        axis.text = element_text(size=10), axis.title = element_text(size=20),
        legend.direction = "vertical", legend.box = "vertical",strip.text.x=element_text(size=20))+
  labs(fill="Chlorophyll Concentration, mg/m-3")+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")+labs(title ="Sampling Stations and Mean Chla October 2023")

chloct<-chl_map+geom_point(data=october, mapping=aes(x=lon_dd, y=lat_dd))+geom_raster(data = not_sea, aes(x = x, y = y))


chlaug+chloct+plot_annotation()
#to fit all plots together#####
library(patchwork)
(saug+soct)/(chlaug+chloct)+plot_annotation(tag_levels=c("A","B","C","D"))

#ssh#####
#data source: https://www.nnvl.noaa.gov/StoryMaps/DITC/ENSO/getdata_SSH_Anomaly.html
ssh<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/SSHA.monthly.202308.bbox=32,182,18,-155.csv")
#convert to 360
ssh<-ssh%>%mutate(x=ifelse(Longitude< 0, Longitude + 360, Longitude))
ssh_plot<-ggplot()+geom_tile(data=ssh,aes(x=x,y=Latitude, fill=Data.Value))+
  scale_fill_gradientn(colors=c("#00007F", "white", "#7F0000"),na.value="gray90",limits = c(-200,200))+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm"))+
  labs(fill="SSH Anomaly", title = "August")+
  theme_bw()+geom_sf()
ssha<-ssh_plot+geom_raster(data = not_sea, aes(x = x, y = y))

ssh<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/SSHA.monthly.202310.bbox=32,182,18,-155.csv")
ssh<-ssh%>%mutate(x=ifelse(Longitude< 0, Longitude + 360, Longitude))
ssh_plot<-ggplot()+geom_tile(data=ssh,aes(x=x,y=Latitude, fill=Data.Value))+
  scale_fill_gradientn(colors=c("#00007F", "white", "#7F0000"),na.value="gray90",limits = c(-200,200))+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm"))+
  labs(fill="SSH Anomaly", title = "October")+
  theme_bw()+geom_sf()
ssho<-ssh_plot+geom_point(data=october, mapping=aes(x=lon_dd, y=lat_dd))+geom_raster(data = not_sea, aes(x = x, y = y))

library(patchwork)
# Your original layout
layout <- (saug + soct) / (chlaug + chloct) / (ssha + ssho)

# Add this line to shrink the margins on ALL plots
layout_tight <- layout & theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

# Add your final annotation
final_plot <- layout_tight + plot_annotation(tag_levels = c("A", "B", "C", "D", "E", "F"))

print(final_plot)
