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
#read.csv(C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/env_pplots_data_w_sse.csv")
legend_title<-"SST (°C)"
swchlInfo <- as.info("noaacrwsstDaily", url="https://coastwatch.noaa.gov/erddap/")

SST_Match <- rxtracto(swchlInfo, parameter = 'analysed_sst', 
                      xcoord = ik$lon_dd, ycoord = ik$lat_dd, tcoord =ik$date , 
                      xlen = .2, ylen = .2, progress_bar = TRUE)
ik$sst<-SST_Match$`mean analysed_sst`
#chla####
legend_title <-expression("Chlorophyll Concentration" * "\n" * "(mg" ~ m^{-3} ~ ")")
#read.csv("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/env_pplots_data_w_chla.csv")
swchlInfo <- as.info("noaacwN20VIIRSchlaDaily", url="https://coastwatch.noaa.gov/erddap/")
#swchlInfo <- as.info("noaacwNPPN20S3ASCIDINEOF2kmDaily", url="https://coastwatch.noaa.gov/erddap/")
#"Chlorophyll (Gap-filled DINEOF), NOAA S-NPP NOAA-20 VIIRS and Copernicus S-3AOLCI, Science Quality, Global 2km, 2018-recent, Daily"
chla_Match <- rxtracto(swchlInfo, parameter = 'chlor_a', 
                      xcoord = ik$lon_dd, ycoord = ik$lat_dd, tcoord =ik$date, 
                      zcoord = rep(0, length(ik$lon_dd)),
                      zName = "altitude",
                      xlen = .2, ylen = .2, progress_bar = TRUE)
ik$chla_daily<-chla_Match$`mean chlor_a`
#write.csv(ik, file="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Data/env_pplots_data_w_chla.csv")

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
  mutate(density=(total/(mean_flow/100)))%>%
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
species_names<-c("Aprion"="virescens","Acanthocybium" = "solandri",
                 "Coryphaena"="spp.", "Katsuwonus"="pelamis","Thunnus"="spp.")

taxa_to_plot <- list(
  family = c("Scombridae", "Istiophoridae"),
  subfamily = c("Etelinae"),
  genus = c("Coryphaena", "Aprion", "Katsuwonus", "Thunnus", "Acanthocybium")
)
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
    
    # --- Step 4: Create and save the plots ---
    hist<-ggplot(data=df, aes(x=length))+geom_histogram(stat="bin", binwidth = 1)+
      labs(x= "Standard Length (mm)",y=paste("frequency of ", fish_taxa))+theme_bw()+
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
                         labels = combined_labels_two_line,limits = c(min_doy, max_doy)) +
      theme_bw()+theme(plot_theme)
    ggsave(paste0(fish_taxa, "_density_vs_doy.png"), plot = p17, width = 4, height = 4,
           path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
    
    # Plot 2: Density vs SST
    p2 <- ggplot(data = df, aes(x = sst, y = density)) +geom_point(size=3) +
      labs(x =expression("SST (" * degree * "C)"),
           y = bquote(paste(.(fish_taxa)," Density (fish 100" ~ m^{-3} ~ ")")))+theme_bw()+
      theme(plot_theme)
    
    ggsave(paste0(fish_taxa, "_density_vs_sst.png"), plot = p2, width = 4, height = 4,path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
    
    # Plot 3: Density vs Chla
    p3 <- ggplot(data = df, aes(x = log(chla_daily), y = density)) +geom_point(size=3)+
      labs(x = expression("Log (e) Chlorophyll Concentration" * "\n" * "(mg" ~m^{-3}~")"),
           y = bquote(paste(.(fish_taxa),"  Density (fish 100" ~ m^{-3} ~")")))+theme_bw()+
      theme(plot_theme)
    ggsave(paste0(fish_taxa, "_density_vs_chla.png"), plot = p3, width = 4, height = 4,path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
   
    
     map_plot<-ggplot() +
      geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +labs(fill = "Depth (m)")+
      scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4",limits=c(-10000,0),
                          na.value=(color="#003300"),guide="none")+
      geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density),shape=19,alpha=0.8) +
       geom_point(data=anti_df,mapping=aes(y=lat_dd, x=lon_dd), shape=1, size=2)+
       coord_sf(xlim=c(178.4,204.9), ylim=c(18.6, 31.4)) +
      scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1)) +
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
        annotate("text", x=178.4, y=18.35,size=3, color="black", hjust=0,
                label=paste("Stations where no",(fish_taxa),"were found"))+
       #taxa name/label
       annotate("text", x=200, y=31.5,size=4, color="black", hjust=0,
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
       path = "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/combined figures/"
     )
     }
}


#hiceas sampling map#####
fig1<-ggplot() +
  geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +labs(fill = "Depth (m)")+
  scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4",limits=c(-10000,0),
                      na.value=(color="#003300"),guide="none")+
  geom_point(data=se23,mapping=aes(y=lat_dd, x=lon_dd,color=month),shape=19,size=3) +
  scale_color_manual(values = c("chartreuse4","#800074","#c99b38","#996600"),labels = c("July","Aug","Oct","Nov"),name="Collection Month") +
  coord_sf(xlim=c(178.4,204.9), ylim=c(18.6, 31.4)) +
  scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1)) +
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

#env in points###########
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
    
p1<-ggplot() +
  geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +labs(fill = "Depth (m)")+
  scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4",limits=c(-10000,0),
                      na.value=(color="#003300"),guide="none")+
  geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density, color=chla_daily),shape=19,alpha=0.8) +
  scale_color_gradientn(colors=c("gray90","darkseagreen1","aquamarine4"),na.value="gray 90",limits=c(0.01,0.2))+
    coord_sf(xlim=c(178.4,204.9), ylim=c(18.6, 31.4)) +
  scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1)) +
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
  annotate("text", x=178.4, y=18.35,size=3, color="black", hjust=0,
           label=paste("Stations where no",(fish_taxa),"were found"))+
  #taxa name/label
  annotate("text", x=200, y=31.5,size=4, color="black", hjust=0,
           fontface = ifelse(level == "genus", "bold.italic", "bold"),
           label=paste(full_scientific_name))+
  annotate("point", x=177.8, y=18.35, color="black",shape=1, size=2)+
  ylab(expression("Latitude (" * degree *"N)"))+xlab(expression("Longitude (" * degree * "W)"))

p2<-ggplot() +
  geom_raster(data = Hawaii_df, aes(x = x, y = y, fill = z)) +labs(fill = "Depth (m)")+
  scale_fill_gradient(high = "lightskyblue1", low = "lightsteelblue4",limits=c(-10000,0),
                      na.value=(color="#003300"),guide="none")+
  geom_point(data=df,mapping=aes(y=lat_dd, x=lon_dd, size=density, color=sst),shape=19,alpha=0.8) +
  scale_color_gradientn(colors=c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", 
                                "yellow", "#FF7F00", "red", "#7F0000"),
                       na.value="gray90",limits = c(26, 29))+
    coord_sf(xlim=c(178.4,204.9), ylim=c(18.6, 31.4)) +
  scale_x_continuous(breaks = round(seq(min(Hawaii_df$x), max(Hawaii_df$x), by = 5),1)) +
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
  annotate("text", x=178.4, y=18.35,size=3, color="black", hjust=0,
           label=paste("Stations where no",(fish_taxa),"were found"))+
  #taxa name/label
  annotate("text", x=200, y=31.5,size=4, color="black", hjust=0,
           fontface = ifelse(level == "genus", "bold.italic", "bold"),
           label=paste(full_scientific_name))+
  annotate("point", x=177.8, y=18.35, color="black",shape=1, size=2)+
  ylab(expression("Latitude (" * degree *"N)"))+xlab(expression("Longitude (" * degree * "W)"))
  pp<-p1/p2
  # 5. Save the final plot
  ggsave(
    plot =pp, # Use the object with tags and layout
    width = 8,                 # Use a reasonable, often slightly smaller width
    height = 10,               # Use a larger height to accommodate the three rows
    dpi = 300,
    filename = paste0(fish_taxa, "_sst_chla_map.png"),
    path = "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Map Figures/"
  )

  }
}

