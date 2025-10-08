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
PMNM<-read_sf("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/PMNM_Shp_File/hi_noaa_nwhi_papahanaumokuakea.shp")
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
legend_title<-"SST (°C)"
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
saug<-ggplot()+ geom_tile(data=dfreal,aes(x=x,y=y, fill=Freq))+
  scale_fill_gradientn(colors=c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", 
                                "yellow", "#FF7F00", "red", "#7F0000"),
                       na.value="gray90",limits = c(26, 29))+
  guides(fill = guide_colorbar(title = NULL,barwidth = unit(0.5, "cm"), barheight = unit(2.5, "cm"))) +
  theme(legend.title=element_text(size=10),legend.text=element_text(size=10),
        axis.text = element_text(size=5), axis.title = element_text(size=10),
        legend.direction = "vertical", legend.box = "vertical",strip.text.x=element_text(size=10))+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")+
  labs(fill=legend_title)
  
saug<-saug+geom_raster(data = not_sea, aes(x = x, y = y))
#map august points##
#load points form hiceas_sample_map.R
#august<-se23%>%filter(month==8)
#saug<-sst_plot+geom_point(data=august, mapping=aes(x=lon_dd, y=lat_dd))+geom_raster(data = not_sea, aes(x = x, y = y))
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
soct<-ggplot()+ geom_tile(data=dfreal,aes(x=x,y=y, fill=Freq))+
  scale_fill_gradientn(colors=c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F",
                                "yellow", "#FF7F00", "red", "#7F0000"),
                       na.value="gray90",limits = c(26, 29))+
  guides(fill = guide_colorbar(title = NULL,barwidth = unit(0.5, "cm"), barheight = unit(2.5, "cm"))) +
  theme(legend.title=element_text(size=10),legend.text=element_text(size=10),
        axis.text = element_text(size=5), axis.title = element_text(size=10),
        legend.direction = "vertical", legend.box = "vertical",strip.text.x=element_text(size=10))+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")
soct+geom_raster(data = not_sea, aes(x = x, y = y))
#october<-se23%>%filter(month==10)
#geom_point(data=october, mapping=aes(x=lon_dd, y=lat_dd))
soct<-sst_plot+geom_raster(data = not_sea, aes(x = x, y = y))
#ggsave(filename="Sampling Stations and Mean SST October 2023.png",plot=soct,height=5, width=8, units="in", dpi=300)
#bind sst######
library(patchwork)
#saug+soct+plot_annotation(tag_levels=c("A","B"))


#August chla###########
#data extraction####
#junk <- GET('https://oceanwatch.pifsc.noaa.gov/erddap/griddap/noaa_snpp_chla_monthly.nc?chlor_a%5B(2023-08-01T12:00:00Z):1:(2023-08-31T12:00:00Z)%5D%5B(18):1:(32)%5D%5B(177):1:(206)%5D',write_disk("junk.nc", overwrite=TRUE))
junk <- GET('https://oceanwatch.pifsc.noaa.gov/erddap/griddap/noaa_snpp_chla_weekly_750m_v2023.nc?chlor_a%5B(2023-08-01T12:00:00Z):1:(2023-08-31T12:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(32):1:(18)%5D%5B(177):1:(206)%5D',write_disk("junk.nc", overwrite=TRUE))

nc <- nc_open('junk.nc')
names(nc$var)
v1 <- nc$var[[1]]
chla<- ncvar_get(nc,v1) #Extract analysed_sst, reads data from the netCDF file, only works if you have already opened the file, shows as a multi-dimensional array
dim(chla)
dates <- as.POSIXlt(v1$dim[[3]]$vals,origin='1970-01-01',tz='GMT') #get the dates for each time step
lon <- v1$dim[[1]]$vals #gives vector of longitude
#lon=lon-360
lat <- v1$dim[[2]]$vals #gives vector of latitude
nc_close(nc) #this step is important, otherwise you risk data loss
rm(junk,v1)
file.remove('junk.nc')

df=apply(chla[,,1:(dim(chla)[3])],c(1,2),mean,na.rm=TRUE) #calculates mean over all dates of the 3rd element of chla. In this case it is the avgerage value over X days
rownames(df, do.NULL = TRUE, prefix = "row")
rownames(df) <- lon
colnames(df, do.NULL = TRUE, prefix = "col")
colnames(df) <-lat
dfreal<-as.data.frame(as.table(df))
dfreal<-dfreal%>%rename(y=Var1,x=Var2)
#mapping####
legend_title <-expression("Chlorophyll Concentration" * "\n" * "(mg" ~ m^{-3} ~ ")")
chla_min<-(quantile(dfreal$Freq, na.rm=T)[1])
chla_max<-(quantile(dfreal$Freq, na.rm=T)[5])

chlaug<-ggplot()+ geom_tile(data=dfreal,aes(x=x,y=y, fill=Freq))+
  scale_fill_gradientn(colors=c("gray90","darkseagreen1","aquamarine4"),na.value="gray 90",limits=c(chla_min, chla_max))+
  guides(
    fill = guide_colorbar(
      title = NULL, 
      barwidth = unit(0.5, "cm"), 
      barheight = unit(2.5, "cm")
    )
  ) +
  theme(legend.title=element_text(size=10),legend.text=element_text(size=10),
        axis.text = element_text(size=5), axis.title = element_text(size=10),
        legend.direction = "vertical", legend.box = "vertical",strip.text.x=element_text(size=10))+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")
chlaug  
  
chlaug<-chlaug+geom_raster(data = not_sea, aes(x = x, y = y))
#october chla###########
#data extraction####
#junk <- GET('https://oceanwatch.pifsc.noaa.gov/erddap/griddap/noaa_snpp_chla_monthly.nc?chlor_a%5B(2023-10-01T12:00:00Z):1:(2023-10-31T12:00:00Z)%5D%5B(18):1:(32)%5D%5B(177):1:(206)%5D',write_disk("junk.nc", overwrite=TRUE))
junk <- GET('https://oceanwatch.pifsc.noaa.gov/erddap/griddap/noaa_snpp_chla_monthly.nc?chlor_a%5B(2023-10-01T12:00:00Z):1:(2023-10-31T12:00:00Z)%5D%5B(18):1:(32)%5D%5B(0.01875000000001137177):1:(206)%5D',write_disk("junk.nc", overwrite=TRUE))

nc <- nc_open('junk.nc')
names(nc$var)
v1 <- nc$var[[1]]
chla<- ncvar_get(nc,v1) #Extract analysed_sst, reads data from the netCDF file, only works if you have already opened the file, shows as a multi-dimensional array
dim(chla)
dates <- as.POSIXlt(v1$dim[[3]]$vals,origin='1970-01-01',tz='GMT') #get the dates for each time step
lon <- v1$dim[[1]]$vals #gives vector of longitude
#lon=lon-360
lat <- v1$dim[[2]]$vals #gives vector of latitude
nc_close(nc) #this step is important, otherwise you risk data loss
rm(junk,v1)
file.remove('junk.nc')

df=apply(chla[,,1:31],c(1,2),mean,na.rm=TRUE) #calculates mean ssh over the 31 vlaues of the 3rd element of ssh. In this case it is the avgerage sst over 31 days
rownames(df, do.NULL = TRUE, prefix = "row")
rownames(df) <- lon
colnames(df, do.NULL = TRUE, prefix = "col")
colnames(df) <-lat
dfreal<-as.data.frame(as.table(df))
dfreal<-dfreal%>%rename(y=Var1,x=Var2)
#mapping####
legend_title <-expression("Chlorophyll Concentration" * "\n" * "(mg" ~ m^{-3} ~ ")")
wrapped_legend_title <- str_wrap(legend_title, width = 15) # Adjust width as needed

chloct<-ggplot()+ geom_tile(data=dfreal,aes(x=x,y=y, fill=Freq))+
  scale_fill_gradientn(colors=c("gray90","darkseagreen1","aquamarine4"),na.value="gray 90")+
  guides(
    fill = guide_colorbar(
      title = NULL, # <--- THIS IS THE KEY CHANGE
      barwidth = unit(0.5, "cm"), 
      barheight = unit(2.5, "cm")
    )
  ) +
  theme(legend.title=element_text(size=10),legend.text=element_text(size=10),
        axis.text = element_text(size=5), axis.title = element_text(size=10),
        legend.direction = "vertical", legend.box = "vertical",strip.text.x=element_text(size=10))+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")

chloct<-chloct+geom_raster(data = not_sea, aes(x = x, y = y))


#to fit all plots together#####

#(saug+soct)/(chlaug+chloct)+plot_annotation(tag_levels=c("A","B","C","D"))

#ssh#####

junk <- GET("https://oceanwatch.pfeg.noaa.gov/noaa_wide/erddap/griddap/noaacwBLENDEDsshDaily.nc?sla%5B(2023-08-01T00:00:00Z):1:(2023-08-31T00:00:00Z)%5D%5B(18):1:(32)%5D%5B(-179.875):1:(-155)%5D",
            write_disk("junk.nc", overwrite=TRUE))
nc <- nc_open('junk.nc')
names(nc$var)
v1 <- nc$var[[1]]
sst<- ncvar_get(nc,v1) #Extract analysed_sst, reads data from the netCDF file, only works if you have already opened the file, shows as a multi-dimensional array
dim(sst)
dates <- as.POSIXlt(v1$dim[[3]]$vals,origin='1970-01-01',tz='GMT') #get the dates for each time step
lon <- v1$dim[[1]]$vals #gives vector of longitude
lon=lon+360
lat <- v1$dim[[2]]$vals #gives vector of latitude
nc_close(nc) #this step is important, otherwise you risk data loss
rm(junk,v1)
file.remove('junk.nc')

df=apply(sst[,,1:31],c(1,2),mean,na.rm=TRUE) #calculates mean ssh over the 31 vlaues of the 3rd element of ssh. In this case it is the avgerage sst over 31 days
rownames(df, do.NULL = TRUE, prefix = "row")
rownames(df) <- lon
colnames(df, do.NULL = TRUE, prefix = "col")
colnames(df) <-lat
# Your data frame dfreal is created here
dfreal <- as.data.frame(as.table(df))
dfreal <- dfreal %>% rename(x = Var1, y = Var2)

# THE FIX: Convert x and y columns to numeric
# The `as.numeric()` function will ensure they are treated as continuous values
dfreal$x <- as.numeric(as.character(dfreal$x))
dfreal$y <- as.numeric(as.character(dfreal$y))

# Now, create the corrected plot
legend_title <- "SSH Anomaly (m)"
ssha<- ggplot() +
  # Add landmass layer first using geom_tile()
  geom_tile(data = not_sea, aes(x = x, y = y), fill = "grey50") +
  
  # Add the SSH data layer next
  geom_tile(data = dfreal, aes(x = x, y = y, fill = Freq)) +
  
  # Add the border path
  geom_path(data = re_ordered_EEZ, aes(x = Lon, y = X2), color = "white", linewidth = 1) +
  
  # Set the fill color scale
  scale_fill_gradientn(
    colors = c("#00007F", "white", "#7F0000"),
    na.value = "gray90",
    limits = c(-0.3, 0.3) 
  ) +
  guides(fill = guide_colorbar(title = NULL,barwidth = unit(0.5, "cm"), barheight = unit(2.5, "cm"))) +
  # Set the plot coordinates
  coord_sf(
    xlim = c(177, 206),
    ylim = c(18, 32),
    expand = FALSE,
    default_crs = NULL 
  ) +
  
  # Add other plot elements
  theme_bw() +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude") +
  labs(fill = legend_title)

ssha<-ssha+geom_raster(data = not_sea, aes(x = x, y = y))#+geom_point(data=august, mapping=aes(x=lon_dd, y=lat_dd))

#ggsave(filename="Sampling Stations and Mean SSH Aug 2023.png",plot=ssha,height=5, width=8, units="in", dpi=300)


junk <- GET('https://oceanwatch.pfeg.noaa.gov/noaa_wide/erddap/griddap/noaacwBLENDEDsshDaily.nc?sla%5B(2023-10-01T00:00:00Z):1:(2023-10-31T00:00:00Z)%5D%5B(18):1:(32)%5D%5B(-179.875):1:(-155)%5D',
            write_disk("junk.nc", overwrite=TRUE))
nc <- nc_open('junk.nc')
names(nc$var)
v1 <- nc$var[[1]]
sst<- ncvar_get(nc,v1) #Extract analysed_sst, reads data from the netCDF file, only works if you have already opened the file, shows as a multi-dimensional array
dim(sst)
dates <- as.POSIXlt(v1$dim[[3]]$vals,origin='1970-01-01',tz='GMT') #get the dates for each time step
lon <- v1$dim[[1]]$vals #gives vector of longitude
lon=lon+360
lat <- v1$dim[[2]]$vals #gives vector of latitude
nc_close(nc) #this step is important, otherwise you risk data loss
rm(junk,v1)
file.remove('junk.nc')

df=apply(sst[,,1:31],c(1,2),mean,na.rm=TRUE) #calculates mean ssh over the 31 values of the 3rd element of ssh. In this case it is the average ssh over 31 days
rownames(df, do.NULL = TRUE, prefix = "row")
rownames(df) <- lon
colnames(df, do.NULL = TRUE, prefix = "col")
colnames(df) <-lat
# Your data frame dfreal is created here
dfreal <- as.data.frame(as.table(df))
dfreal <- dfreal %>% rename(x = Var1, y = Var2)

# THE FIX: Convert x and y columns to numeric
# The `as.numeric()` function will ensure they are treated as continuous values
dfreal$x <- as.numeric(as.character(dfreal$x))
dfreal$y <- as.numeric(as.character(dfreal$y))

# Now, create the corrected plot
ssho<- ggplot() +
  # Add landmass layer first using geom_tile()
  geom_tile(data = not_sea, aes(x = x, y = y), fill = "grey50") +
  
  # Add the SSH data layer next
  geom_tile(data = dfreal, aes(x = x, y = y, fill = Freq)) +
  
  # Add the border path
  geom_path(data = re_ordered_EEZ, aes(x = Lon, y = X2), color = "white", linewidth = 1) +
  
  # Set the fill color scale
  scale_fill_gradientn(
    colors = c("#00007F", "white", "#7F0000"),
    na.value = "gray90",
    limits = c(-0.3, 0.3) 
  ) +
  guides(
    fill = guide_colorbar(
      title = NULL, # Use the wrapped title here
      title.hjust = 0.5,
      barwidth = unit(0.5, "cm"),
      barheight = unit(2.5, "cm")
    )
  ) +  # Set the plot coordinates
  coord_sf(
    xlim = c(177, 206),
    ylim = c(18, 32),
    expand = FALSE,
    default_crs = NULL 
  ) +
  
  # Add other plot elements
  theme_bw() +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude") 
ssho<-ssho+geom_raster(data = not_sea, aes(x = x, y = y))

library(patchwork)
#all 6 layout##########
ggsave("sst_aug_uniform.png", plot = saug , width = 6, height = 4)
ggsave("sst_oct_uniform.png", plot = soct, width = 6, height = 4)
ggsave("chla_aug_uniform.png", plot = chlaug, width = 6, height = 4)
ggsave("chla_oct_uniform.png", plot = chloct, width = 6, height = 4)
ggsave("ssh_august_uniform.png", plot = ssha, width = 6, height = 4)
ggsave("ssh_october_uniform.png", plot = ssho, width = 6, height = 4)

combined_plots <-(saug + soct) / (chlaug + chloct) / (ssha + ssho)

# Add this line to shrink the margins on ALL plots
layout_tight <-  combined_plots&
  theme(
    plot.margin = unit(c(0.001, 0.1, 0.1, 0.001), "cm"),
    plot.tag.position = c(0, 1),
    plot.tag = element_text(size = 10, face = "bold")
  )
# Add your final annotation
final_plot <- layout_tight +
  plot_annotation(tag_levels = 'A')
final_plot
# Save the complete, combined plot
#ggsave(filename = "final_combined_plot.png",path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Map Figures/",plot = final_plot,
 #width = 8,  # Set the overall width for the entire figure
#  height = 10,  # Set the overall height
 # dpi = 300)


#MODELLLED ENV DATA###########
#add chla from Jessie!!!!!####
library(raster)
library(tidyverse)
library(sp)
library(sf)
library(reshape2)
#August####
#read in monthly 4km ocean color data (.nc file)
data = "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/20230801-20230831_cmems_obs-oc_glo_bgc-plankton_myint_l4-olci-4km_P1M.nc"
chlor <-stack(data, varname = "CHL")
# # use buffered bathymetry file to mask shallow waters around MHI
load("C:/Users/Andrea.Schmidt/Documents/billfish_not_github/bathy_HI_30m_buffered_poly_1km.RData")
chlor_masked <- raster::mask(chlor, bathy_buffered_poly_1km, inverse = TRUE)

dfreal <- terra::as.data.frame(chlor_masked, xy = TRUE, na.rm = TRUE)
names(dfreal) <- c("x", "y", "Freq")
dfreal<-dfreal%>%
  mutate(x=ifelse(x<0, x+360, x))
print(paste("Data frame dimensions:", dim(dfreal)[1], "rows and", dim(dfreal)[2], "columns"))

#map
legend_title <-expression("Chlorophyll Concentration" * "\n" * "(mg" ~ m^{-3} ~ ")")
wrapped_legend_title <- str_wrap(legend_title, width = 15) # Adjust width as needed
chlaug<-ggplot()+ geom_tile(data=dfreal,aes(x=x,y=y, fill=Freq))+
  scale_fill_gradientn(colors=c("gray90","darkseagreen1","aquamarine4"),na.value="gray 90",limits=c(0,0.1))+
  guides(
    fill = guide_colorbar(
      title = NULL, # <--- THIS IS THE KEY CHANGE
      barwidth = unit(0.5, "cm"), 
      barheight = unit(2.5, "cm")
    )
  ) +
  theme(legend.title=element_text(size=10),legend.text=element_text(size=10),
        axis.text = element_text(size=5), axis.title = element_text(size=10),
        legend.direction = "vertical", legend.box = "vertical",strip.text.x=element_text(size=10))+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")

chlaug

##October
data = "C:/Users/Andrea.Schmidt/Desktop/Cruises/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/20231001-20231031_cmems_obs-oc_glo_bgc-plankton_myint_l4-olci-4km_P1M.nc"
chlor <-stack(data, varname = "CHL")
# # use buffered bathymetry file to mask shallow waters around MHI
load("C:/Users/Andrea.Schmidt/Documents/billfish_not_github/bathy_HI_30m_buffered_poly_1km.RData")
chlor_masked <- raster::mask(chlor, bathy_buffered_poly_1km, inverse = TRUE)

dfreal <- terra::as.data.frame(chlor_masked, xy = TRUE, na.rm = TRUE)
names(dfreal) <- c("x", "y", "Freq")

print(paste("Data frame dimensions:", dim(dfreal)[1], "rows and", dim(dfreal)[2], "columns"))

#map
legend_title <-expression("Chlorophyll Concentration" * "\n" * "(mg" ~ m^{-3} ~ ")")
wrapped_legend_title <- str_wrap(legend_title, width = 15) # Adjust width as needed
chlaug<-ggplot()+ geom_tile(data=dfreal,aes(x=x,y=y, fill=Freq))+
  scale_fill_gradientn(colors=c("gray90","darkseagreen1","aquamarine4"),na.value="gray 90")+
  guides(
    fill = guide_colorbar(
      title = NULL, # <--- THIS IS THE KEY CHANGE
      barwidth = unit(0.5, "cm"), 
      barheight = unit(2.5, "cm")
    )
  ) +
  theme(legend.title=element_text(size=10),legend.text=element_text(size=10),
        axis.text = element_text(size=5), axis.title = element_text(size=10),
        legend.direction = "vertical", legend.box = "vertical",strip.text.x=element_text(size=10))+
  theme_bw()+coord_sf(xlim = c(206, 177), ylim = c(18, 32), expand=FALSE)+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  theme(axis.ticks.length = unit(0.25, "cm")) +
  ylab("Latitude") +
  xlab("Longitude")
chloct