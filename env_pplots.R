###chla, sst per station###
# libraries####
library(tidyverse)      # Includes dplyr, tidyr, ggplot2, stringr (and is sufficient for most of your code)
library(lubridate)
library(rerddap)
library(rerddapXtracto)
library(sf)
library(TSP)
library(mgcv)
library(mgcViz)
library(visreg)
library(patchwork)
library(ncdf4)
library(httr)
#station data#####
#coordinate conversion##
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
#sst data extraction####
#DONT USE THIS ONE, IT DROPPED 6 STATIONS: ikd<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/env_pplots_data_w_sse.csv")#ik[,1]<-NULL
#legend_title<-"SST (°C)"
swchlInfo <- as.info("noaacrwsstDaily", url="https://coastwatch.noaa.gov/erddap/")
SST_Match <- rxtracto(swchlInfo, parameter = 'analysed_sst', xcoord = ik$lon_dd, ycoord = ik$lat_dd, tcoord =ik$date , xlen = .2, ylen = .2, progress_bar = TRUE)
ik$sst<-SST_Match$`mean analysed_sst`
#chla####
legend_title <-expression("Chlorophyll Concentration" * "\n" * "(mg" ~ m^{-3} ~ ")")
chla_df<-read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/chla_extraction_glorys_4km.csv")
chla_df<-chla_df%>%mutate(date=ymd(date))
ik2<-ik%>%dplyr::mutate(date=lubridate::ymd(date))
library(data.table)
ik2_dt<-ik2
chla_dt<-chla_df
setDT(ik2_dt)
setDT(chla_dt)
TOLERANCE_DEGREE <- 0.03# This is ~3.3 km, ensuring you cover the 4km cell
ik2_dt[, `:=`(ik_date = as.Date(date),
  lat_min = lat_dd - TOLERANCE_DEGREE,
  lat_max = lat_dd + TOLERANCE_DEGREE,
  lon_min = lon_dd - TOLERANCE_DEGREE,
  lon_max = lon_dd + TOLERANCE_DEGREE)][, date := NULL]

chla_dt[, `:=`(chla_date = as.Date(date))][, date := NULL]
setkey(chla_dt, chla_date)
spatial_matches_dt<- chla_dt[ik2_dt,
                  on = c("lat > lat_min","lat < lat_max", # 1. LATITUDE RANGE MATCH: chla_dt$lat must be between ik2_dt$lat_min and ik2_dt$lat_max
                         "lon > lon_min","lon < lon_max"),# 2. LONGITUDE RANGE MATCH: chla_dt$lon must be between ik2_dt$lon_min and ik2_dt$lon_max
                  allow.cartesian = TRUE,nomatch=0]
final_ik3_dt <- spatial_matches_dt %>%as_tibble() %>% 
  mutate(date_diff = abs(chla_date - ik_date)) %>%
  filter(date_diff <= days(2)) %>%
  group_by(ik_date, lat_dd, lon_dd) %>%
  filter(date_diff == min(date_diff)) %>%
  slice(1) %>% 
  ungroup() %>%
  select(Station,chla_value, date_diff)
ik3<-left_join(ik2,final_ik3_dt, by=join_by("Station"=="Station"))
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
se23<-left_join(ik3, larv,by=join_by("Station"=="station"))
se23<-se23%>%
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
  mutate(density=(total/(mean_flow/100)))%>%
  dplyr::select(Leg, Station,lat_dd,lon_dd,mean_flow, month, depth_m,total, toyos,tow_duration,
                sst,
                chla_value,
                FLAG, family, subfamily,genus, species,density,day_of_year)

min_global_sst=min(se23$sst,na.rm=T)
max_global_sst=max(se23$sst,na.rm=T)
global_chla_min=min(se23$chla_value,na.rm=T)
global_chla_max=max(se23$chla_value,na.rm=T)
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
Hawaii_df <-fortify(Hawaii)
print("despite the above error message Hawaii_df still loads as needed")

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
species_names<-c("Aprion"="virescens","Acanthocybium" = "solandri",
                 "Coryphaena"="spp.", "Katsuwonus"="pelamis","Thunnus"="spp.")

taxa_to_plot <- list(
  family = c("Scombridae", "Istiophoridae"),
  subfamily = c("Etelinae"),
  genus = c("Coryphaena", "Aprion", "Katsuwonus", "Thunnus", "Acanthocybium"))
#map density values####
tuna <- se23 %>% filter(family == "Scombridae")
marlin<- se23 %>% filter(family == "Istiophoridae")
jobfish<- se23 %>% filter(subfamily == "Etelinae")
mahi<- se23 %>% filter(genus == "Coryphaena")
uku<- se23 %>% filter(genus == "Aprion")
aku<- se23 %>% filter(genus == "Katsuwonus")
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
tuna <- larv %>% filter(family == "Scombridae")
marlin<- larv %>% filter(family == "Istiophoridae")
jobfish<- larv %>% filter(subfamily == "Etelinae")
mahi<- larv %>% filter(genus == "Coryphaena")
uku<- larv %>% filter(genus == "Aprion")
aku<- larv %>% filter(genus == "Katsuwonus")
ahi<- larv %>% filter(genus == "Thunnus")
ono<- larv %>% filter(genus == "Acanthocybium")
larv2 <-rbind(tuna, marlin,jobfish,mahi,uku, aku,ahi, ono)

mas_long<-larv2%>%
  unite("taxa",sep="_", family:genus)%>%
  dplyr::select(c(station, taxa,count_1_mm:count_40_mm))%>%
  pivot_longer(count_1_mm:count_40_mm,names_to="length",
               values_to="length_occurence",values_drop_na = T) #pivot longer &coalesce both lengths and frequencies columns or #melt, shape, reshape

mas_clean<-mas_long%>%
  mutate("length"=gsub("mm","",length), .keep="unused")%>%
  mutate("length"=gsub("count","",length))%>%
  mutate("length"=gsub("_","",length))%>%
  mutate("length"=as.numeric(length))
mas_long_clean<-mas_clean%>%filter(length_occurence!=0)%>%#distinct function eliminated vials where multiple size classes were counted
  mutate(station=as.factor(station))
se23<-se23%>%mutate(station=as.factor(Station))
se23<-left_join(mas_long_clean, se23,by=join_by("station"=="station"),relationship == "many-to-many")#missing values here, but also df is WAYYY too long
global_length_min=min(se23$length,na.rm=T)
global_length_max=max(se23$length,na.rm=T)
#plot theme####
plot_theme<-c(axis.title = element_text(size=10),axis.text=element_text(size=10))
map_theme<-c(axis.title = element_text(size=11),axis.text=element_text(size=10),
             axis.ticks.length = unit(0.25, "cm"),
             plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
             panel.grid = element_blank(),text = element_text(size=10),
             legend.position = "inside",legend.position.inside=c(0.22,0.03),
             legend.justification = c(1, 0),legend.background = element_rect(fill = "transparent"),
             legend.key = element_rect(fill = "transparent"))
#for loop####
for (level in names(taxa_to_plot)) {
  for (fish_taxa in taxa_to_plot[[level]]) {
    
    # --- Step 3: Filter the data for the current taxon ---
    #allow for access to full scientific name
    current_species <- species_names[fish_taxa]
    current_species_name <- ifelse(is.na(current_species),"", paste0(current_species))
    full_scientific_name <- paste(fish_taxa, current_species_name)
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
    summary_filename <- paste0(fish_taxa, "_summary_stats.txt")
    full_path <- file.path("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/", summary_filename)
      capture.output(
      cat("--- Summary Statistics for: ", full_scientific_name, " ---\n\n"),
      print(summary(df)),
      file = full_path,
      append = FALSE)
    # --- Step 4: Create and save the plots ---
    hist<-ggplot(data=df, aes(x=length))+
      geom_histogram(stat="bin", binwidth = 1)+
      labs(x= "Standard Length (mm)",y=paste("frequency of ", fish_taxa))+theme_bw()+
      scale_x_continuous(breaks = seq(global_length_min,global_length_max, by=1),
                         limits = c(global_length_min, global_length_max)) +
      theme(plot_theme)
    ggsave(paste0(fish_taxa, "_length_histogram.png"), plot = hist, width = 4, height = 4,
           path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
#day of year with month box on top
     p16 <- ggplot(data = df, aes(x = day_of_year, y =density)) +
      geom_point(size = 3)+
      geom_boxplot(outliers = F, alpha = 0.8,aes(group=Leg), fill="pink", color="red")+
      labs(x = "Day of Year",
           y = bquote(paste(.(fish_taxa)," Density (fish 100" ~ m^{-3} ~ ")")))+
      scale_x_continuous(breaks = doy_combined_breaks,labels = combined_labels_two_line,limits = c(min_doy, max_doy)) +
      theme_bw()+theme(plot_theme)


    ggsave(paste0(fish_taxa, "_density_vs_doy_boxplot.png"), plot = p16, width = 4, height = 4,
           path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
 #day of year only
    p17 <- ggplot(data = df, aes(x = day_of_year, y =density)) +
      geom_point(size = 3)+
      labs(x = "Day of Year",y = bquote(paste(.(fish_taxa)," Density (fish 100" ~ m^{-3} ~ ")")))+
      scale_x_continuous(breaks = doy_combined_breaks,
                         labels = combined_labels_two_line,
                         limits=c(min(se23$day_of_year, na.rm=T),max(se23$day_of_year, na.rm=T))) +
      theme_bw()+theme(plot_theme)
    ggsave(paste0(fish_taxa, "_density_vs_doy.png"), plot = p17, width = 4, height = 4,
           path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
    
    # Plot 2: Density vs SST
    p2 <- ggplot(data = df, aes(x = sst, y = density)) +geom_point(size=3)+
      scale_x_continuous(breaks = seq(min_global_sst, max_global_sst, by=0.5),
                         labels = function(x) { return(as.character(round(x, digits = 1)))})+
      #scale_y_continuous(breaks=seq(GLOBAL_MIN_DENSITY,GLOBAL_MAX_DENSITY, by=0.),limits=c(GLOBAL_MIN_DENSITY,GLOBAL_MAX_DENSITY),labels = function(y) { return(as.character(round(y, digits = 3)))})+
      labs(x =expression("SST (" * degree * "C)"),
           y = bquote(paste(.(fish_taxa)," Density (fish 100" ~ m^{-3} ~ ")")))+theme_bw()+
      theme(plot_theme)
    
    ggsave(paste0(fish_taxa, "_density_vs_sst.png"), plot = p2, width = 4, height = 4,path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
    
    # Plot 3: Density vs Chla
    p3 <- ggplot(data = df, aes(x =chla_value, y = density)) +
      geom_point(size=3)+
      scale_x_continuous(breaks = seq(global_chla_min,global_chla_max,by=0.01), 
                         limits=c(global_chla_min,global_chla_max),
                         labels = function(x) { return(as.character(round(x, digits = 3)))})+
      #scale_y_continuous(breaks=seq(GLOBAL_MIN_DENSITY,GLOBAL_MAX_DENSITY, by=0.),limits=c(GLOBAL_MIN_DENSITY,GLOBAL_MAX_DENSITY),labels = function(y) { return(as.character(round(y, digits = 3)))})+
    labs(x = expression("Chlorophyll Concentration" * "\n" * "(mg" ~m^{-3}~")"),
         y = bquote(paste(.(fish_taxa),"  Density (fish 100" ~ m^{-3} ~")")))+
      theme_bw()+theme(plot_theme)
    ggsave(paste0(fish_taxa, "_density_vs_chla.png"), plot = p3, width = 4, height = 4,path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
   
    
     map_plot<-ggplot() +
      geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +labs(fill = "Depth (m)")+
      scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4",limits=c(-10000,0),
                          na.value=(color="#003300"),guide="none")+
      geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density, color=day_of_year),shape=19,alpha=0.8) +
       #geom_point(data=anti_df,mapping=aes(y=lat_dd, x=lon_dd), shape=1, size=2)+
       coord_sf(xlim=c(178.4,204.9), ylim=c(18.6, 31.4)) +
      scale_x_continuous(breaks = round(seq(min(Hawaii_df$x, na.rm=T), max(Hawaii_df$x,na.rm=T), by = 5),1),
                         labels = function(x) { x - 360 })+
       scale_color_gradientn(colors=c("yellow","orange","orangered","red3","red4"),
                             na.value="gray 90",limits=c(min(se23$day_of_year, na.rm=T),max(se23$day_of_year, na.rm=T)))+
       scale_size_continuous(limits = c(fish_min_quant,fish_max_quant),
                             range = c(MIN_POINT_SIZE, MAX_POINT_SIZE),
                             breaks = DENSITY_BREAKS,
                             name=bquote("Density(fish 100" ~ m^{-3}~")"))+
       theme_bw()+theme(axis.title = element_text(size=11),axis.text=element_text(size=10),
                        axis.ticks.length = unit(0.25, "cm"),
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        panel.grid = element_blank(),text = element_text(size=10),
                        legend.position = "inside",legend.position.inside=c(0.22,0.03),
                        legend.justification = c(1, 0),legend.background = element_rect(fill = "transparent"),
                        legend.key = element_rect(fill = "transparent"))+
       geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
      #antidf label
        #annotate("text", x=178.4, y=18.35,size=3, color="black", hjust=0,label=paste("Stations where no",(fish_taxa),"were found"))+
       #taxa name/label
       annotate("text", x=190, y=31.5,size=4, color="black", hjust=0,
                fontface = ifelse(level == "genus", "bold.italic", "bold"),
                label=paste(full_scientific_name))+
       annotate("point", x=177.8, y=18.35, color="black",shape=1, size=2)+
      ylab(expression("Latitude (" * degree *"N)"))+xlab(expression("Longitude (" * degree * "W)"))
     map_plot
     ggsave(paste0(fish_taxa, "_map.png"),
            plot = map_plot,
            path = "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/",
            width = 8, height = 6)
 
#legend.position c(0.95, 0.1),This moves the legend anchor point to the bottom-right area. 0.95 is near the right edge of the plot (1.0 is the far right), and 0.1 is just above the bottom edge (0.0 is the far bottom).
#legend.justification	c(1, 0), 	This tells ggplot to use the bottom-right corner of the legend box as the anchor point. 1 is right-justified, and 0 is bottom-justified. This ensures the legend box expands up and to the left from the specified position
#legend.background	element_rect(fill = "transparent"),	This is crucial for placing the legend over the map background without blocking it with a white or grey box, making it visually integrate better.    
    
  #Combining all 5 plots together  
     # 1. Define the rows, setting widths inside the row definition
     # The middle row contains the Boxplot (p17) and Histogram (hist)
     p_middle_row <- p17 + hist + plot_layout(widths = c(1, 1)) # Use equal widths (1:1)
     
     # The bottom row contains the SST (p2) and Chla (p3) scatterplots
     p_bottom_row <- p2 + p3 + plot_layout(widths = c(1, 1)) # Use equal widths (1:1)
     
     # 2. Combine all three rows vertically using the / operator
     final_plot_combined <- map_plot / p_middle_row / p_bottom_row
     
     # 3. Apply the layout and styling to the combined plot object
     final_plot_layout <- final_plot_combined +
       # Set the relative heights for the three rows (Map, Middle Row, Bottom Row)
       # Map is usually tallest (e.g., 3 units), the other rows are equal (e.g., 1 unit each)
       plot_layout(heights = c(3, 1.5, 1.5)) +
       
       # Apply all consistent styling to all sub-plots using the & operator
       # We will use your tight margin and tag settings here
       theme(
         plot.margin = unit(c(0.1, 0.01, 0.01, 0.001), "cm"),
         plot.tag.position = c(0, 1),
         plot.tag = element_text(size = 10, face = "bold")
       )
     
     # 4. Add plot tags (A, B, C, etc.)
     final_plot_tagged <- final_plot_layout +
       plot_annotation(tag_levels = 'A')
     
     # 5. Save the final plot
     ggsave(
       plot = final_plot_tagged, # Use the object with tags and layout
       width = 8,                 # Use a reasonable, often slightly smaller width
       height = 10,               # Use a larger height to accommodate the three rows
       dpi = 300,
       filename = paste0(fish_taxa, "_final_combined_plot.png"),
       path = "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/combined figures/")
     }
}


#hiceas sampling map#####
fig1<-ggplot() +
  geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +labs(fill = "Depth (m)")+
  scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4",limits=c(-10000,0),
                      na.value=(color="#003300"),guide="none")+
  geom_point(data=se23,mapping=aes(y=lat_dd, x=lon_dd,fill=month), color="black",shape=19,size=3) +
  scale_fill_manual(values = c("chartreuse4","#800074","#c99b38","#996600"),labels = c("July","Aug","Oct","Nov"),name="Collection Month") +
  coord_sf(xlim=c(178.4,204.9), ylim=c(18.6, 31.4)) +
  scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1),
                     labels = function(x) { x - 360 }) +
  theme_bw()+theme(axis.title = element_text(size=11),axis.text=element_text(size=10),
                   axis.ticks.length = unit(0.25, "cm"),
                   plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                   panel.grid = element_blank(),text = element_text(size=10),
                   legend.position = "inside",legend.position.inside=c(1,0.59),
                   legend.justification = c(1, 0),legend.background = element_rect(fill = "white"),
                   legend.key = element_rect(fill = "white"))+
  geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
  labs(y=expression("Latitude (" * degree *"N)"),
       x=(expression("Longitude (" * degree * "W)")))
fig1
ggsave("fig1.png",
       plot = fig1,
       path = "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Map Figures/",
       width = 8, height = 6)

#doy/ density##########
for (level in names(taxa_to_plot)) {
  for (fish_taxa in taxa_to_plot[[level]]) {
    
    # --- Step 3: Filter the data for the current taxon ---
    #allow for access to full scientific name
    current_species <- species_names[fish_taxa]
    current_species_name <- ifelse(is.na(current_species),"", paste0(current_species))
    full_scientific_name <- paste(fish_taxa, current_species_name)
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
    
    pp<-ggplot() +
      geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +labs(fill = "Depth (m)")+
      scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4",limits=c(-10000,0),
                          na.value=(color="#003300"),guide="none")+
      geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density, color=day_of_year),shape=19,alpha=0.8) +
      #geom_point(data=anti_df,mapping=aes(y=lat_dd, x=lon_dd), shape=1, size=2)+
      scale_color_gradientn(colors=c("yellow","orange","orangered","red3","red4"),
                            na.value="gray 90",limits=c(min(se23$day_of_year),max(se23$day_of_year)))+
      coord_sf(xlim=c(178.4,204.9), ylim=c(18.6, 31.4)) +
      scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1),
                         labels = function(x) { x - 360 }) +
      scale_size_continuous(limits = c(fish_min_quant,fish_max_quant),
                            range = c(MIN_POINT_SIZE, MAX_POINT_SIZE),
                            breaks = DENSITY_BREAKS,
                            name=bquote("Density(fish 100" ~ m^{-3}~")"))+
      theme_bw()+theme(axis.title = element_text(size=11),axis.text=element_text(size=10),
                       axis.ticks.length = unit(0.25, "cm"),
                       plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                       panel.grid = element_blank(),text = element_text(size=10),
                       legend.position = "inside",legend.position.inside=c(0.22,0.03),
                       legend.justification = c(1, 0),legend.background = element_rect(fill = "transparent"),
                       legend.key = element_rect(fill = "transparent"))+
      geom_path(data = re_ordered_EEZ,aes(x = Lon, y = X2), color = "white",linewidth = 1)+
      #antidf label
      #annotate("text", x=178.4, y=18.35,size=3, color="black", hjust=0,label=paste("Stations where no",(fish_taxa),"were found"))+
      #taxa name/label
      annotate("text", x=190, y=31.5,size=4, color="black", hjust=0,
               fontface = ifelse(level == "genus", "bold.italic", "bold"),
               label=paste(full_scientific_name))+
      annotate("point", x=177.8, y=18.35, color="black",shape=1, size=2)+
      ylab(expression("Latitude (" * degree *"N)"))+xlab(expression("Longitude (" * degree * "W)"))
    ggsave(
      plot =pp, # Use the object with tags and layout
      width = 8,                 # Use a reasonable, often slightly smaller width
      height = 10,               # Use a larger height to accommodate the three rows
      dpi = 300,
      filename = paste0(fish_taxa, "_doy_color_map.png"),
      path = "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Map Figures/")
    
  }
}
#stats#####
#univariate pearson corellations as further diagnostics####
summary_file <- file.path("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/cor_tests.txt")
sink(summary_file) # Redirects R output to the specified file
ggplot()+geom_point(data=se23, aes(x=log(chla_value), y=sst))
cor.test(se23$sst, log(se23$chla_value),method="pearson")
ggplot()+geom_point(data=se23, aes(x=chla_value, y=sst))
cor.test(se23$sst, se23$chla_value,method="pearson") 
cor.test(se23$sst, se23$density,method="pearson") 
cor.test(as.numeric(se23$station), se23$density,method="pearson")
cor.test(se23$day_of_year, se23$density,method="pearson") 
cor.test(se23$depth_m, se23$density,method="pearson") 
cor.test(se23$family, se23$density,method="pearson") 
cor.test(log(se23$chla_value), se23$density,method="pearson") 
cor.test(se23$chla_value,se23$density,method="pearson")
cor.test(se23$lat_dd,se23$density,method="pearson")
cor.test(se23$lon_dd,se23$density,method="pearson")
sink()
library(coefplot)
# 1. Load Necessary Libraries
#JNP:for the models, assuming zero-inflation, the family should be tweedie(link="log").
# I would run new models with that family, and can you also print the text from the gam.check 
#diagnostics so we can see the effective degrees of freedom for each? 
#concern is that with such low sample sizes,we may be running up against too many predictors 
library(mgcv)
library(mgcViz)
library(visreg)
library(coefplot)
output_dir <- "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Stat Figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
for (level in names(taxa_to_plot)) {
  for (fish_taxa in taxa_to_plot[[level]]) {
    
    # 3. Create Clean Names for Files
    # Replace spaces and special characters for a safe file name
    clean_taxa_name <- gsub(" ", "_", fish_taxa) 
    
    current_species <- species_names[fish_taxa]
    current_species_name <- ifelse(is.na(current_species), "", paste0(current_species))
    full_scientific_name <- paste(fish_taxa, current_species_name)
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
    df$pa<-"present"
    anti_df$pa<-"absent"
    anti_df$density<-0
    df<-rbind(df,anti_df)
    # Check if the filtered data frame is empty to avoid errors
    if (nrow(df) == 0) {
      next # Skip to the next iteration if no data is found
    }
    
    # Check if the filtered data frame is empty to avoid errors
    if (nrow(df) == 0) {
      cat(paste("Skipping", full_scientific_name, "- No data found.\n"))
      next # Skip to the next iteration if no data is found
    }
    #Count non-zero density observations
    non_zero_density <- sum(df$density > 0, na.rm = TRUE)
    min_observations <- 5 # Set a threshold for minimum non-zero observations
    
    if (non_zero_density < min_observations) {
      cat(paste("Skipping", full_scientific_name, 
                "- Only", non_zero_density, 
                "non-zero density observations. Need at least", min_observations, ".\n"))
      # Save a brief note in the summary file to document the skip
      summary_file <- file.path(output_dir, paste0(clean_taxa_name, "_Model_Summary.txt"))
      sink(summary_file)
      cat(paste("Skipped analysis: Insufficient non-zero density observations (n=", non_zero_density, ") for GAM.\n"))
      sink()
      next # Skip to the next taxa
    }
    
    
    cat(paste("--- Analyzing:", full_scientific_name, "---\n"))
    
    # 7. Open the PDF graphics device BEFORE tryCatch
    # If this fails (rare), the tryCatch won't even run.
    pdf_file <- file.path(output_dir, paste0(clean_taxa_name, "_Diagnostic_Plots.pdf"))
    pdf(pdf_file, width = 10, height = 7)
    tryCatch({
      
      # 4. Fit the Gaussian GAM Model
      model <- gam(density ~ s(sst, k=4) + s(log(chla_value), k=4) +
                     factor(month) + te(lon_dd, lat_dd, k=4), 
                   family = gaussian(link="identity"), data = df)
    
    # 5. Save Model Summary and ANOVA to a Text File using sink()
    summary_file <- file.path(output_dir, paste0(clean_taxa_name, "_Model_Summary.txt"))
    sink(summary_file) # Redirects R output to the specified file
    cat(paste("GAM Analysis for:", full_scientific_name, "\n\n"))
    cat("---------- MODEL SUMMARY ----------\n")
    print(summary(model))
    cat("\n\n---------- ANOVA TABLE ----------\n")
    print(anova(model))
    cat("\n\n---------- CORRELATION SST vs DENSITY ----------\n")
    print(cor.test(df$sst, df$density, method="pearson"))
    cat("n\n-------------Pearson Cortests vs density----------\n")
    cor.test(df$sst, log(df$chla_value),method="pearson")
    cor.test(df$sst, df$chla_value,method="pearson") 
    cor.test(df$sst, df$density,method="pearson") 
    cor.test(as.numeric(df$station), df$density,method="pearson")
    cor.test(df$day_of_year, df$density,method="pearson") 
    cor.test(df$depth_m, df$density,method="pearson") 
    cor.test(df$family, df$density,method="pearson") 
    cor.test(log(df$chla_value), df$density,method="pearson") 
    cor.test(df$chla_value,df$density,method="pearson")
    cor.test(df$lat_dd,df$density,method="pearson")
    cor.test(df$lon_dd,df$density,method="pearson")
    sink() # Closes the redirection, sending output back to the console
    cat(paste("Saved statistical summary to:", summary_file, "\n"))
    
    # 6. Create mgcViz object and assign dynamically
    viz_object_name <- paste0(clean_taxa_name, "_fig")
    assign(viz_object_name, getViz(model)) 
    mahifig <- get(viz_object_name)
    dummy_plot_output <- print(plot(mahifig, allTerms=TRUE), pages = 1)
  
    # mgcViz plots (using the dynamically created object)
    print(plot(mahifig, allTerms=TRUE), pages = 1)
    
    # Standard GAM Checks
    gam.check(model)
    
    # Standard R Diagnostic Plots (Model 1-4)
    plot(model) 
    
    # Residual Plots
    qqnorm(resid(model), main = paste("Q-Q Plot of Residuals -", full_scientific_name)) 
    qqline(resid(model))
    plot(resid(model) ~ fitted(model), main = paste("Residuals vs. Fitted -", full_scientific_name), 
         xlab = "Fitted Values", ylab = "Residuals")
    hist(resid(model), main = paste("Histogram of Residuals -", full_scientific_name))
    
    # Influence/Outlier Plots
    plot(cooks.distance(model), type="h", main = paste("Cook's Distance -", full_scientific_name)) 
    
    # Visualization of SST effect
    visreg(model, "sst", scale = "response", ylab = "Density (Response Scale)", 
           xlab = "SST (deg.C)", alpha=0.01, rug=1, 
           main = paste("Partial Effect of SST -", full_scientific_name))
    }, error = function(e) {
      # This runs if an error occurs. 
      cat(paste("--- ERROR DURING ANALYSIS for", full_scientific_name, ":", conditionMessage(e), "---\n"))
      
      # Append the error message to the summary file
      summary_file <- file.path(output_dir, paste0(clean_taxa_name, "_Model_Summary.txt"))
      sink(summary_file, append = TRUE)
      cat(paste("\n\n*** ERROR STOPPED MODELING ***\n", conditionMessage(e), "\n"))
      sink()
      
      # We do NOT use dev.off() here. The 'finally' block handles it.
      
      # We must explicitly return to ensure the loop moves on gracefully
      return(NULL) 
      
    }, finally = {
      if (dev.cur() != 1) { # Close the active device, but only if one is open (dev.cur() != 1)
        Sys.sleep(0.1) # Add a tiny sleep to allow I/O operations to complete
        dev.off() # Close the currently active device (should be the PDF)
        cat(paste("Successfully closed graphics device.\n"))
      } else {
        cat("No active device to close.\n")
      }
      
      cat(paste("Saved output files to:", pdf_file, "\n"))
    })}}



output_dir <- "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Stat Figures"
##WIP cor tests###############
str(se23)
se23<-se23%>%unite("taxa",sep="_", family:genus, remove=F)
str(larv)
df <- larv %>% filter(family == "Scombridae")
anti_df <- larv %>% filter(family != "Scombridae")
df<- larv %>% filter(family == "Istiophoridae")
anti_df <- larv %>% filter(family != "Istiophoridae")
df<- larv %>% filter(subfamily == "Etelinae")
anti_df <- larv %>% filter(subfamily != "Etelinae")
df<- larv %>% filter(genus == "Coryphaena")
anti_df <- larv %>% filter(genus != "Coryphaena")
df<- larv %>% filter(genus == "Aprion")
anti_df <- larv %>% filter(genus != "Aprion")
df<- larv %>% filter(genus == "Katsuwonus")
anti_df <- larv %>% filter(genus != "Katsuwonus")
df<- larv %>% filter(genus == "Thunnus")
anti_df <- larv %>% filter(genus != "Thunnus")
df<- larv %>% filter(genus == "Acanthocybium")
anti_df <- larv %>% filter(genus != "Acanthocybium")

    df$pa <- 1
    anti_df$pa <- 0
    anti_df$density <- 0
    df <- rbind(df, anti_df)

    model <- gam(pa ~ s(sst, k=4) + s(log(chla_value), k=4) +
                   factor(month) + te(lon_dd, lat_dd, k=4), 
                 family = binomial(link="logit"), data = df)
    
    coefplot(model,  title = paste0(clean_taxa_name," Coefficient Plot"))

        print(cor.test(df$sst, log(df$chla_value), method = "pearson"))
    print(cor.test(df$sst, df$chla_value, method = "pearson")) 
    print(cor.test(df$sst, df$density, method = "pearson")) 
    print(cor.test(as.numeric(df$station), df$density, method = "pearson"))
    print(cor.test(df$day_of_year, df$density, method = "pearson")) 
    print(cor.test(df$depth_m, df$density, method = "pearson")) 
    # Note: cor.test on two non-numeric/factor columns will likely fail.
    # print(cor.test(df$family, df$density, method = "pearson")) 
    print(cor.test(log(df$chla_value), df$density, method = "pearson")) 
    print(cor.test(df$chla_value, df$density, method = "pearson"))
    print(cor.test(df$lat_dd, df$density, method = "pearson"))
    print(cor.test(df$lon_dd, df$density, method = "pearson"))
    


#binomial#####
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
for (level in names(taxa_to_plot)) {
  for (fish_taxa in taxa_to_plot[[level]]) {
    
    # 3. Create Clean Names for Files
    # Replace spaces and special characters for a safe file name
    clean_taxa_name <- gsub(" ", "_", fish_taxa) 
    
    current_species <- species_names[fish_taxa]
    current_species_name <- ifelse(is.na(current_species), "", paste0(current_species))
    full_scientific_name <- paste(fish_taxa, current_species_name)
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
    df$pa<-1#present
    anti_df$pa<-0#absent
    df<-rbind(df,anti_df)
    # Check if the filtered data frame is empty to avoid errors
    if (nrow(df) == 0) {
      next # Skip to the next iteration if no data is found
    }
    
    # Check if the filtered data frame is empty to avoid errors
    if (nrow(df) == 0) {
      cat(paste("Skipping", full_scientific_name, "- No data found.\n"))
      next # Skip to the next iteration if no data is found
    }
    cat(paste("--- Analyzing:", full_scientific_name, "---\n"))
    pdf_file <- file.path(output_dir, paste0(clean_taxa_name, "_Binomial_Diagnostic_Plots.pdf"))
    pdf(pdf_file, width = 10, height = 7)
    tryCatch({
#4.5 Fit binomial model
model <- gam(pa ~ s(sst, k=4) + s(log(chla_value), k=4) +
               factor(month) + te(lon_dd, lat_dd, k=4), 
             family = binomial(link="logit"), data = df)

# 5. Save Model Summary and ANOVA to a Text File using sink()
summary_file <- file.path(output_dir, paste0(clean_taxa_name, "_Binomial Model_Summary.txt"))
sink(summary_file) # Redirects R output to the specified file
cat(paste("binomial GAM Analysis for:", full_scientific_name, "\n\n"))
cat("---------- MODEL SUMMARY ----------\n")
print(summary(model))
cat("\n\n---------- ANOVA TABLE ----------\n")
print(anova(model))
cat("\n\n---------- CORRELATION SST vs pa ----------\n")
print(cor.test(df$sst, df$pa, method="pearson"))
cat("\n\n---------- CORRELATION CHLA vs pa ----------\n")
print(cor.test(df$chla_value, df$pa, method="pearson"))
cat("\n\n---------- GAM CHECK ----------\n")
print(gam.check(model))
sink() # Closes the redirection, sending output back to the console
cat(paste("Saved statistical summary to:", summary_file, "\n"))

# 6. Create mgcViz object and assign dynamically
viz_object_name <- paste0(clean_taxa_name, "_fig")
assign(viz_object_name, getViz(model)) 
mahifig <- get(viz_object_name)
dummy_plot_output <- print(plot(mahifig, allTerms=TRUE), pages = 1)

# mgcViz plots (using the dynamically created object)
print(plot(mahifig, allTerms=TRUE), pages = 1)

# Standard GAM Checks
gam.check(model)

# Standard R Diagnostic Plots (Model 1-4)
plot(model) 

# Residual Plots
qqnorm(resid(model), main = paste("Q-Q Plot of Residuals -", full_scientific_name)) 
qqline(resid(model))
plot(resid(model) ~ fitted(model), main = paste("Residuals vs. Fitted -", full_scientific_name), 
     xlab = "Fitted Values", ylab = "Residuals")
hist(resid(model), main = paste("Histogram of Residuals -", full_scientific_name))

# Influence/Outlier Plots
plot(cooks.distance(model), type="h", main = paste("Cook's Distance -", full_scientific_name)) 

# Visualization of SST effect
visreg(model, "sst", scale = "response", ylab = "presence/absence (Response Scale)", 
       xlab = "SST (deg.C)", alpha=0.01, rug=1, 
       main = paste("Partial Effect of SST -", full_scientific_name))

# Visualization of chla effect
visreg(model, "daily chla", scale = "response", ylab = "presence/absence (Response Scale)", 
       xlab = "chla mg/-m3)", alpha=0.01, rug=1, 
       main = paste("Partial Effect of chla -", full_scientific_name))

  }, error = function(e) {
    # This runs if an error occurs. 
    cat(paste("--- ERROR DURING ANALYSIS for", full_scientific_name, ":", conditionMessage(e), "---\n"))
    
    # Append the error message to the summary file
    summary_file <- file.path(output_dir, paste0(clean_taxa_name, "_Model_Summary.txt"))
    sink(summary_file, append = TRUE)
    cat(paste("\n\n*** ERROR STOPPED MODELING ***\n", conditionMessage(e), "\n"))
    sink()
    
    # We do NOT use dev.off() here. The 'finally' block handles it.
    
    # We must explicitly return to ensure the loop moves on gracefully
    return(NULL) 
    
  }, finally = {
    if (dev.cur() != 1) { # Close the active device, but only if one is open (dev.cur() != 1)
      Sys.sleep(0.1) # Add a tiny sleep to allow I/O operations to complete
      dev.off() # Close the currently active device (should be the PDF)
      cat(paste("Successfully closed graphics device.\n"))
    } else {
      cat("No active device to close.\n")
    }
    
    cat(paste("Saved output files to:", pdf_file, "\n"))
  })}}



#tweedie####
#Diagnostic Plots for the Tweedie Model (tw_model)
#tw model useful for examining pa and abundance
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
for (level in names(taxa_to_plot)) {
  for (fish_taxa in taxa_to_plot[[level]]) {
    
    # 3. Create Clean Names for Files
    # Replace spaces and special characters for a safe file name
    clean_taxa_name <- gsub(" ", "_", fish_taxa) 
    
    current_species <- species_names[fish_taxa]
    current_species_name <- ifelse(is.na(current_species), "", paste0(current_species))
    full_scientific_name <- paste(fish_taxa, current_species_name)
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
    df$pa<-"present"
    anti_df$pa<-"absent"
    anti_df$density<-0
    df<-rbind(df,anti_df)
    # Check if the filtered data frame is empty to avoid errors
    if (nrow(df) == 0) {
      next # Skip to the next iteration if no data is found
    }
    
    # Check if the filtered data frame is empty to avoid errors
    if (nrow(df) == 0) {
      cat(paste("Skipping", full_scientific_name, "- No data found.\n"))
      next # Skip to the next iteration if no data is found
    }
    #Count non-zero density observations
    non_zero_density <- sum(df$density > 0, na.rm = TRUE)
    min_observations <- 5 # Set a threshold for minimum non-zero observations
    
    if (non_zero_density < min_observations) {
      cat(paste("Skipping", full_scientific_name, 
                "- Only", non_zero_density, 
                "non-zero density observations. Need at least", min_observations, ".\n"))
      # Save a brief note in the summary file to document the skip
      summary_file <- file.path(output_dir, paste0(clean_taxa_name, "_Model_Summary.txt"))
      sink(summary_file)
      cat(paste("Skipped analysis: Insufficient non-zero density observations (n=", non_zero_density, ") for GAM.\n"))
      sink()
      next # Skip to the next taxa
    }
    
    
cat(paste("--- Running gam.check for Tweedie Model ---\n"))
# 7. Open the PDF graphics device BEFORE tryCatch
# If this fails (rare), the tryCatch won't even run.
pdf_file <- file.path(output_dir, paste0(clean_taxa_name, "_Tweedie_Diagnostic_Plots.pdf"))
pdf(pdf_file, width = 10, height = 7)
tryCatch({
  
model <- gam(density ~ s(sst, k=4) + s(log(chla_value), k=4) +
                  factor(month) + te(lon_dd, lat_dd, k=4), 
                family = tw(link="log"), data = df)

# 5. Save Model Summary and ANOVA to a Text File using sink()
summary_file <- file.path(output_dir, paste0(clean_taxa_name, "_Tweedie_Model_Summary.txt"))
sink(summary_file) # Redirects R output to the specified file
cat(paste("GAM Analysis for:", full_scientific_name, "\n\n"))
cat("---------- MODEL SUMMARY ----------\n")
print(summary(model))
cat("\n\n---------- ANOVA TABLE ----------\n")
print(anova(model))
cat("\n\n---------- CORRELATION SST vs DENSITY ----------\n")
print(cor.test(df$sst, df$density, method="pearson"))
cat("\n\n---------- GAM CHECK ----------\n")
print(gam.check(model))
sink() # Closes the redirection, sending output back to the console
cat(paste("Saved statistical summary to:", summary_file, "\n"))



# Standard GAM Diagnostic Plots (Model 1-4 for Deviance Residuals)
plot(model, main=paste("Tweedie GAM Diagnostics -", full_scientific_name))

#gam.check
gam.check(model)

# Visualization of SST effect for the Tweedie Model
visreg(model, "sst", scale = "response", ylab = "Density (Tweedie Response Scale)",
       xlab = "SST (deg.C)", alpha=0.01, rug=1,
       main = paste("Tweedie Partial Effect of SST -", full_scientific_name))
}, error = function(e) {
  # This runs if an error occurs. 
  cat(paste("--- ERROR DURING ANALYSIS for", full_scientific_name, ":", conditionMessage(e), "---\n"))
  
  # Append the error message to the summary file
  summary_file <- file.path(output_dir, paste0(clean_taxa_name, "_Model_Summary.txt"))
  sink(summary_file, append = TRUE)
  cat(paste("\n\n*** ERROR STOPPED MODELING ***\n", conditionMessage(e), "\n"))
  sink()
  
  # We do NOT use dev.off() here. The 'finally' block handles it.
  
  # We must explicitly return to ensure the loop moves on gracefully
  return(NULL) 
  
}, finally = {
  if (dev.cur() != 1) { # Close the active device, but only if one is open (dev.cur() != 1)
    Sys.sleep(0.1) # Add a tiny sleep to allow I/O operations to complete
    dev.off() # Close the currently active device (should be the PDF)
    cat(paste("Successfully closed graphics device.\n"))
  } else {
    cat("No active device to close.\n")
  }
  
  cat(paste("Saved output files to:", pdf_file, "\n"))
})}}


#NON-fish###########
load("A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Hawaii.RData")
Hawaii_df <-marmap::fortify.bathy(Hawaii)

not_fish<-read.csv("A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/SE2303 ID Log - Non-fish counts.csv")
library(tidyverse)
library(ggrepel)
ik<-read.csv("A:/My Drive/crusies/HICEAS_23/IKMT Log - Sheet1.csv")
#convert ship coordinates into DMS then to DD
ik<-ik %>% mutate(lat_deg= as.numeric(substr(Lat.In,1,2)))#http://127.0.0.1:40641/graphics/4ee7ccc4-d91e-43e0-8b31-d427019a6d2c.png
ik<-ik %>% mutate(lon_deg= as.numeric(substr(Lon.In,1,4)))
ik<-ik %>% mutate(lat_minsec= as.numeric(substr(Lat.In,3,10)))
ik<-ik %>% mutate(long_minsec= as.numeric(substr(Lon.In,5,10)))
ik<-ik %>% mutate(lat_dd= lat_deg+(lat_minsec/60))
ik<-ik %>% mutate(lon_dd= lon_deg+(long_minsec/60))
ik<-ik%>%filter(lon_dd<1000)
ik<-ik%>%mutate(lon_dd=ifelse(lon_dd<0, lon_dd+360,lon_dd))

ik2<-ik%>%dplyr::select(Station, lat_dd, lon_dd)%>%
  rename(station="Station")
ik3<-ik2%>%left_join(not_fish, ik2, by="station")

ik3 <- ik3 %>%mutate(halo = ifelse(is.na(halobates_count) | as.character(halobates_count) %in% c("0", ""),
      "absent","present"),
    phyl = ifelse(is.na(phyllosome_count) | as.character(phyllosome_count) %in% c("0", ""),
      "absent","present"),
    phys = ifelse(is.na(physallia_count) | as.character(physallia_count) %in% c("0", ""),
      "absent","present"),
    squid = ifelse(is.na(cephalopods_YN) | cephalopods_YN %in% c("N", ""),
      "absent","present"))


ik3_clean_coords <- ik3 %>%
  select(station, lat_dd, lon_dd, halo, phyl, phys, squid) %>%
  distinct(station, .keep_all = TRUE)

ik4 <- ik3_clean_coords %>%
  pivot_longer(cols = c(halo, phyl, phys, squid),
               names_to = "taxa",values_to = "presence_absence")

larv2<-ik3_clean_coords%>%
  group_by(station)%>%
  summarize(halo, phyl, phys, squid)

kable(larv2,col.names =c("Station Number", "Halobates sp.", "Phyllosomes (Palinuridae or Scyllaridae)","Physalia sp.","Cephalopoda"),
      caption = "Presence or absence of non-fish from all stations")%>%
  kable_styling(latex_options = "striped",table.env='table*')#%>%
  save_kable("nonfish.png")

taxa_names<- c("halo"= "Halobates sp.","phyl"="Phyllosomes (Palinuridae or Scyllaridae)",
                      "phys"="Physalia sp.","squid"="Cephalopoda")

map_theme <- theme(axis.title = element_text(size = 11), 
  axis.text = element_text(size = 10),axis.ticks.length = unit(0.25, "cm"),
  plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),strip.text.x = element_text(size = 10),
  panel.grid = element_blank(), text = element_text(size = 10),#legend.position = "inside", 
  #legend.position.inside = c(0.22, 0.03),legend.justification = c(1, 0), 
  legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent"))

map_plot <- ggplot() +
  geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +
  labs(fill = "Depth (m)") +
  scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4", limits = c(-10000, 0),na.value = (color = "#003300"), guide = "none") +
  geom_point(data = ik4, mapping = aes(y = lat_dd, x = lon_dd, shape = presence_absence)) +
  facet_wrap(~taxa, labeller = ggplot2::as_labeller(taxa_names))+
  scale_shape_manual(values = c("absent" = 1, "present" = 19), name = "Observation")+
  coord_sf(xlim = c(178.4, 204.9), ylim = c(18.6, 31.4)) +
    scale_x_continuous(breaks = round(seq(min(Hawaii_df$x, na.rm=T), max(Hawaii_df$x,na.rm=T), by = 5),1),
                       labels = function(x) { x - 360 })+
  theme_bw() + map_theme +ylab(expression("Latitude (" * degree * "N)")) +
  xlab(expression("Longitude (" * degree * "W)"))
ggsave(
  plot =map_plot, # Use the object with tags and layout
  width = 8,                 # Use a reasonable, often slightly smaller width
  height = 10,               # Use a larger height to accommodate the three rows
  dpi = 300,
  filename = "nonfish.png",
  path = "A:/My Drive/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Map Figures/")
