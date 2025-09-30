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
library(sf)
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
library(patchwork)
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
#swchlInfo <- as.info("noaacwNPPN20S3ASCIDINEOF2kmDaily", url="https://coastwatch.noaa.gov/erddap/")
#"Chlorophyll (Gap-filled DINEOF), NOAA S-NPP NOAA-20 VIIRS and Copernicus S-3AOLCI, Science Quality, Global 2km, 2018-recent, Daily"
chla_Match <- rxtracto(swchlInfo, parameter = 'chlor_a', 
                      xcoord = ik$lon_dd, ycoord = ik$lat_dd, tcoord =ik$date, 
                      zcoord = rep(0, length(ik$lon_dd)),
                      zName = "altitude",
                      xlen = .2, ylen = .2, progress_bar = TRUE)
ik$chla_daily<-chla_Match$`mean chlor_a`
#monument boundary####
library(sf)
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
  mutate(month=month(date))%>%
  mutate(month = factor(month, levels = c(7, 8, 10), labels = c("July", "August", "October"))) %>%
  mutate(day_of_year=yday(date))%>%
  mutate(toyos=as.numeric(Tow.yo.Count))%>%
  mutate(total=as.numeric(total))%>%
  mutate(density=(total/mean_flow))%>%
  dplyr::select(Leg, Station,lat_dd,lon_dd,mean_flow, month, depth_m,total, toyos,tow_duration,
                sst,chla_daily,
                FLAG, family, subfamily,genus, species,density,day_of_year)
#doy and month labeling for plots#####

# Define the start of each survey month by Day of Year (DOY)
# This assumes the survey starts on July 23rd for the sake of an example.
july_doy <- yday(as.Date("2023-07-01"))
august_doy <- yday(as.Date("2023-08-01"))
october_doy <- yday(as.Date("2023-10-01"))

# Define the full range of the survey in DOY to ensure consistency
min_doy <- yday(as.Date("2023-07-01"))
max_doy <- yday(as.Date("2023-11-01"))
doy_minor_breaks <- seq(from = min_doy, to = max_doy, by = 5)
# Define the range for the major breaks (e.g., every 10 days)
doy_major_breaks <- seq(from = min_doy, to = max_doy, by = 10)
# Define the range for the minor breaks (e.g., every 5 days)
doy_minor_breaks <- seq(from = min_doy, to = max_doy, by = 5)
#combo plot#####
# Assuming min_doy, max_doy, july_doy, august_doy, october_doy are defined

# 1. Define the Major Break points using a regular 10-day interval
doy_combined_breaks <- seq(from = min_doy, to = max_doy, by = 10)

# 2. Define the exact Day of Year for month starts and their names
month_start_doy <- c(july_doy, august_doy, october_doy)
month_names <- c("July", "August", "October")

# Keep track of which months have already been labeled
months_labeled <- logical(length(month_start_doy)) 

# 3. Create the Two-Line Combined Label Vector
combined_labels_two_line <- character(length(doy_combined_breaks))

for (i in 1:length(doy_combined_breaks)) {
  doy <- doy_combined_breaks[i]
  
  # Day of Year number (always available)
  doy_line <- as.character(doy)
  
  # Month name (default to empty)
  month_line <- "" 
  
  # Check each month start to see if this break is the FIRST one to pass it
  for (m in 1:length(month_start_doy)) {
    if (!months_labeled[m] && doy >= month_start_doy[m]) {
      # This is the first 10-day break ON or AFTER the month start date
      month_line <- month_names[m]
      months_labeled[m] <- TRUE # Mark this month as labeled
      break # Exit the inner loop once a month is labeled for this break
    }
  }
  
  # Combine the lines: DOY on top, Month Name on the bottom (only if a month name exists)
  if (month_line != "") {
    combined_labels_two_line[i] <- paste(doy_line, month_line, sep = "\n")
  } else {
    # Only display the DOY number
    combined_labels_two_line[i] <- doy_line
  }
}
#load maps and taxa name list####
load("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Hawaii.RData")
Hawaii_df <- fortify(Hawaii)

etelinae_expression <- list("Etelinae" = expression("'opakapaka or" * "\n" * "ehu or gindai or" * "\n" * "lehi or onaga or" * "\n" * "kalekale"))

hawaiian_names <- c(
  "Thunnus" = "'ahi",
  "Aprion" = "uku",
  "Acanthocybium" = "ono",
  "Coryphaena" = "mahi mahi",
  "Istiophoridae" = "a'u",
  "Katsuwonus"="aku",
  "Scombridae"= "Tunas & ono",
  "Etelinae"= "'opakapaka/onaga/gindai/ehu/lehi/kalekale")
taxa_to_plot <- list(
  family = c("Scombridae", "Istiophoridae"),
  subfamily = c("Etelinae"),
  genus = c("Coryphaena", "Aprion", "Katsuwonus", "Thunnus", "Acanthocybium")
)
#map density values####
tuna <- se23 %>% filter(family == "Scombridae")
marlin<- se23 %>% filter(family == "Istiophiridae")
jobfish<- se23 %>% filter(subfamily == "Etelinae")
mahi<- se23 %>% filter(genus == "Coryphaena")
uku<- se23 %>% filter(genus == "Aprion")
aku<- se23 %>% filter(genus == "Kaysuwonus")
ahi<- se23 %>% filter(genus == "Thunnus")
ono<- se23 %>% filter(genus == "Acanthocybium")
sword<- se23 %>% filter(genus == "Xiphias")

all_density_data <- list(tuna, marlin,jobfish,mahi,uku, aku,ahi, ono, sword)

all_densities <- unlist(lapply(all_density_data, function(df) df$density))
GLOBAL_MAX_DENSITY <- max(all_densities, na.rm = TRUE)
GLOBAL_MIN_DENSITY <- min(all_densities, na.rm = TRUE)
MIN_POINT_SIZE <- 2
MAX_POINT_SIZE <- 6
#size histograms##########
str(se23)
se23<-se23%>%unite("taxa",sep="_", family:genus, remove=F)

str(larv)
mas_long<-larv%>%
  unite("taxa",sep="_", family:genus)%>%
  dplyr::select(c(station, taxa,count_1_mm:count_40_mm))%>%
  pivot_longer(count_1_mm:count_40_mm,names_to="length",
               values_to="length_occurence",values_drop_na = T) #pivot longer &coalesce both lengths and frequencies columns or #melt, shape, reshape

mas_clean<-mas_long%>%
  mutate("length"=gsub("mm","",length), .keep="unused")%>%
  mutate("length"=gsub("count","",length))%>%
  mutate("length"=gsub("_","",length))%>%
  mutate("length"=as.numeric(length))
mas_long_clean<-mas_clean%>%filter(length_occurence!=0)#distinct function eliminated vials where multiple size classes were counted
se23<-left_join(mas_long_clean, se23,by=c(station="Station", taxa="taxa"),relationship = "many-to-many")

#for loop####
for (level in names(taxa_to_plot)) {
  for (fish_taxa in taxa_to_plot[[level]]) {
    
    # --- Step 3: Filter the data for the current taxon ---
    # The filtering column changes based on the level.
    if (level == "family") {
      df <- se23 %>% filter(family == fish_taxa)
      anti_df <- se23 %>% filter(family != fish_taxa)
    } else if (level == "subfamily") {
      df <- se23 %>% filter(subfamily == fish_taxa)
      anti_df <- se23 %>% filter(subfamily != fish_taxa)
    } else if (level == "genus") {
      df <- se23 %>% filter(genus == fish_taxa)
      anti_df <- se23 %>% filter(genus != fish_taxa)
    }
    
    # Check if the filtered data frame is empty to avoid errors
    if (nrow(df) == 0) {
      next # Skip to the next iteration if no data is found
    }
    
    #find quantile min and max and 1,2,3,4 then make a column so each value that falls into a given quantile is in its own column
    station_quantiles <- df %>%
      reframe(quantile_values = quantile(density, na.rm = TRUE),
              quantile_level = names(quantile(density, na.rm = TRUE))) %>%
      rename(total_count = quantile_values)
    fish_min_quant<-station_quantiles$total_count[1]
    fish_max_quant<-station_quantiles$total_count[5]
    quant_breaks<-c(station_quantiles$total_count[1],station_quantiles$total_count[2], station_quantiles$total_count[3],station_quantiles$total_count[4],station_quantiles$total_count[5])
    DENSITY_BREAKS <- unique(quant_breaks)
    
    #progress bar
    hawaii_name <- hawaiian_names[fish_taxa]
    if (!is.na(hawaii_name)) {
      cat("Processing", fish_taxa, " (", hawaii_name, ")\n")
    } else {
      cat("Processing", fish_taxa, " (Hawaiian name not found)\n")
    }
    
    # --- Step 4: Create and save the plots ---
    hist<-ggplot(data=df, aes(x=length))+geom_histogram(stat="bin", binwidth = 1)+
      labs(x= "Standard Length (mm)",y=paste("frequency of ", fish_taxa))+theme_bw()
    ggsave(paste0(fish_taxa, "_length_histogram.png"), plot = hist, width = 6, height = 5,
           path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")

     p16 <- ggplot(data = df, aes(x = day_of_year, y =density)) +
      geom_point(size = 3)+
      geom_boxplot(outliers = F, alpha = 0.8,aes(group=Leg), fill="pink", color="red")+
      labs(x = "Day of Year",
           y = bquote(paste(.(fish_taxa)," Density (fish" ~ m^{-3} ~ ")")))+
      #adds break every 5 days to x axis while keeping months
      #scale_x_continuous(breaks = c(july_doy, august_doy, october_doy),labels = c("July", "August", "October"),minor_breaks = doy_minor_breaks, limits = c(min_doy, max_doy)) +
      #version 3 combo labels      
      scale_x_continuous(breaks = doy_combined_breaks,labels = combined_labels_two_line,limits = c(min_doy, max_doy)) +
      theme_bw()
    ggsave(paste0(fish_taxa, "_density_vs_doy_boxplot.png"), plot = p16, width = 6, height = 5,
           path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
   
    # Plot 2: Density vs SST
    p2 <- ggplot(data = df, aes(x = sst, y = density)) +geom_point() +labs(x =expression("SST (" * degree * "C)"),y = bquote(paste(.(fish_taxa)," Density (fish" ~ m^{-3} ~ ")")))+theme_bw()
    ggsave(paste0(fish_taxa, "_density_vs_sst.png"), plot = p2, width = 6, height = 5,path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
    
    # Plot 3: Density vs Chla
    p3 <- ggplot(data = df, aes(x = chla_daily, y = density)) +geom_point() +labs(x = expression("Chlorophyll Concentration" * "\n" * "(mg" ~ m^{-3} ~ ")"),y = bquote(paste(.(fish_taxa),"  Density (fish" ~ m^{-3} ~ ")")))+theme_bw()
    ggsave(paste0(fish_taxa, "_density_vs_chla.png"), plot = p3, width = 6, height = 5,path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
   
    # Turn off scientific notation globally
    options(scipen = 5)
    options(digits = 6) 

     map_plot<-ggplot() +
      geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +labs(fill = "Depth (m)")+
      scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4",limits=c(-10000,0),
                          na.value="black",guide="none")+
      geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density),shape=19) +
       geom_point(data=anti_df,mapping=aes(y=lat_dd, x=lon_dd), shape=1, size=2)+
       coord_sf(xlim=c(178.3,204.6), ylim=c(18.6, 31.4)) +
      scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1)) +
       scale_size_continuous(limits = c(fish_min_quant,fish_max_quant),range = c(MIN_POINT_SIZE, MAX_POINT_SIZE),breaks = DENSITY_BREAKS,name=bquote(paste(.(fish_taxa),"  Density (fish" ~ m^{-3} ~ ")")))+
       theme_bw()+theme(axis.ticks.length = unit(0.25, "cm"),panel.grid = element_blank(),text = element_text(size=16))+
      ylab("Latitude")+xlab("Longitude")
    ggsave(paste0(fish_taxa, "_map.png"), plot = map_plot, path = "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/", width = 8, height = 6)
 
    combined_plots<-map_plot/(p16+hist)/(p2+p3)
    layout_tight <-  combined_plots&
      theme(
        plot.margin = unit(c(0.001, 0.01, 0.01, 0.001), "cm"),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 10, face = "bold")
      )
    # Add your final annotation
    final_plot <- layout_tight +
      plot_annotation(tag_levels = 'A')
    final_plot
    ggsave(plot = final_plot,width =8 ,height = 8,dpi = 300,
    filename =paste0(fish_taxa,"_final_combined_plot.png"),path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
     }
}



# To revert to R's default behavior:
options(scipen = 0)
