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
  dplyr::select(Station,lat_dd,lon_dd,mean_flow, month, depth_m,total, toyos,tow_duration,
                sst,chla_daily,FLAG, family, subfamily,genus, species,density,day_of_year)

# Define the start of each survey month by Day of Year (DOY)
# This assumes the survey starts on July 23rd for the sake of an example.
july_doy <- yday(as.Date("2023-07-01"))
august_doy <- yday(as.Date("2023-08-01"))
october_doy <- yday(as.Date("2023-10-01"))

# Define the full range of the survey in DOY to ensure consistency
min_doy <- yday(as.Date("2023-07-01"))
max_doy <- yday(as.Date("2023-11-01"))

# Load necessary data and define plot variables outside the loop
load("C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Hawaii.RData")
Hawaii_df <- fortify(Hawaii)

hawaiian_names <- c(
  "Thunnus" = "ahi",
  "Aprion" = "uku",
  "Acanthocybium" = "ono",
  "Coryphaena" = "mahi mahi",
  "Istiophoridae" = "a'u",
  "Xiphiidae" = "a'u kū"
)
taxa_to_plot <- list(
  family = c("Scombridae", "Istiophoridae","Xiphiidae"),
  subfamily = c("Etelinae"),
  genus = c("Coryphaena", "Aprion", "Katsuwonus", "Thunnus", "Acanthocybium")
)
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
    
    hawaii_name <- hawaiian_names[fish_taxa]
    if (!is.na(hawaii_name)) {
      cat("Processing", fish_taxa, " (", hawaii_name, ")\n")
    } else {
      cat("Processing", fish_taxa, " (Hawaiian name not found)\n")
    }
    
    # --- Step 4: Create and save the plots ---
    
    # Plot 1: Boxplot of density by month
    p1 <- ggplot(data = se23 %>% filter(!!sym(level) == fish_taxa), aes(x = month, y = density)) +
      geom_boxplot(outliers = F) +
      labs(title = paste(fish_taxa, "Density by Month"),
           subtitle = "Outliers removed", y = "Density", x = "Month") +
      theme_bw()
    
    ggsave(paste0(fish_taxa, "_monthly_density.png"), plot = p1, width = 6, height = 5,
           path = "C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
    
    # Plot 1.5: Density vs Day of Year (no changes needed)
    p15 <- ggplot(data = df, aes(x = day_of_year, y = density)) +
      geom_point(size = 3) +
      labs(title = paste(fish_taxa, "Density by Day of Year"),
           y = "Density", x = "Month") +
      scale_x_continuous(breaks = c(yday(as.Date("2023-07-01")), yday(as.Date("2023-08-01")), yday(as.Date("2023-10-01"))),
                         labels = c("July", "August", "October")) +
      theme_bw()
    
    ggsave(paste0(fish_taxa, "_density_vs_doy.png"), plot = p15, width = 6, height = 5,
           path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
    
    p16 <- ggplot(data = df, aes(x = day_of_year, y = density)) +
      geom_point(size = 3) +geom_boxplot(outliers = F, alpha = 0.8,aes(group=month))+
      labs(
        title = paste(fish_taxa, "Density by Day of Year"),
        y = "Density", x = "Month"
      ) +
      # Set the x-axis to be continuous and define the breaks and labels
      scale_x_continuous(
        breaks = c(july_doy, august_doy, october_doy),
        labels = c("July", "August", "October"),
        limits = c(min_doy, max_doy)
      ) +
      theme_bw()
    ggsave(paste0(fish_taxa, "_density_vs_doy_boxplot.png"), plot = p16, width = 6, height = 5,
           path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
    
    # Plot 2: Density vs SST
    p2 <- ggplot(data = df, aes(x = sst, y = density)) +
      geom_point() +
      labs(title = paste(fish_taxa, "Density vs SST"),
           y = "Density", x = "SST (°C)") +
      theme_bw()
    ggsave(paste0(fish_taxa, "_density_vs_sst.png"), plot = p2, width = 6, height = 5,
           path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
    
    # Plot 3: Density vs Chla
    p3 <- ggplot(data = df, aes(x = chla_daily, y = density)) +
      geom_point() +
      labs(title = paste(fish_taxa, "Density vs Chla"),
           y = "Density", x = expression("Chlorophyll Concentration" * "\n" * "(mg" ~ m^{-3} ~ ")")) +
      theme_bw()
    ggsave(paste0(fish_taxa, "_density_vs_chla.png"), plot = p3, width = 6, height = 5,
           path="C:/Users/Andrea.Schmidt/Desktop/crusies/HICEAS_23/Ichthyoplankton Projects/Ichthyoplankton Projects/Plot Figures/")
    # --- Plot 5: Map of density locations with quantile sizes ---
    
    # 1. Inside your main for loop, get the data for the current taxon.
    # We will use the 'se23' data frame and filter it within the loop.
    
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
    
  }
}
