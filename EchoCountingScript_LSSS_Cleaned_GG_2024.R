##########################################################################################
####################                                                  ####################
####################    ESTIMATING FISH DENSITY FROM ECHO-COUNTING    ####################
####################              (Nuuk fjord system)                 ####################
##########################################################################################

## Copyright 2024 Technical University of Denmark
## Author: Giacomo Gardella <giaga@aqua.dtu.dk>


##########################################################################################
## This R script can be used to estimate the fish density in boxes of 0.05 deg latitude ##
## 0.1 deg longitude and 10 m depth using the method of echo-counting. Data needed      ##
## is a report 16 file (Sv data) from LSSS exported with 0.005 nmi (or 0.1 nmi)         ##
## horizontal grid size and 10 m of vertical grid size and with a minimum Sv of -99 dB. ##
## In addition a "TS distribution" file should be exported from LSSS, containing        ##
## the results of the SED algorithm. Finall, a file containg all pings shuld be         ##
## generated in LSSS. This can be done by generating a report 16 with the horizontal    ##
## grid size of 1 ping and the vertical grid size equal to teh max depth encountered.   ##
##########################################################################################

### Getting started

# Load relevant packages           
library("ggplot2")
library('ggpubr')
library(dplyr)
library(jpeg)
library(png)
library(fields)
library(devtools)
library(sp)
library(mgcv)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(VGAM)
library(raster)
library(rnaturalearthhires)

# Define parameters - EBA (from calibration of echo sounder)
alpha_beam <- 10^(-20.7/10)   # -20.7 dB, convert to steradians

# Set working directory
setwd("C:/Giacomo Gardella/Comparisons LSSS")

### Pings - load, rename, aggregate into rectangles
df_pings_2024 = read.table("pings_2024.txt", quote = "", header=T, sep="", dec=".", blank.lines.skip = FALSE, stringsAsFactors=F)
names(df_pings_2024) = c("Year", "Month", "Day", "Time_UTC", "Log1_ping", "Log2_ping", "Lat", "Lon", "BDMIN","BDMAX","OBJECT","CH", "DepthMin_m", "DepthMax_m", "DepthMean_m", "UPINLM", "NASC_UHC","TOTAL") #"OTHER", "NASC_MallotusVilosus","KRILL","TOTAL")
dim(df_pings_2024)
df_pings_2024$Rect_Lat = round(df_pings_2024$Lat*20)/20
df_pings_2024$Rect_Lon = round(df_pings_2024$Lon*10)/10
df_pings_2024$n_pings = 1
df_pings_2024 = aggregate(n_pings~Rect_Lat+Rect_Lon, df_pings_2024, sum)

### Load report 16 of entire fjord - -99 db threshold; Nuuk_Fjord_99db_2 ; Nuuk_Fjord_2024FMNight_01nmi10m_99dB; for higher resolution Nuuk_Fjord_2024FMNight_0005nmi10m_99dB.txt
dfAcousticData_99dB = read.table("Nuuk_Fjord_2024FMNight_0005nmi10m_99dB_Copy.txt", quote = "", header=T, sep="", dec=".", blank.lines.skip = FALSE, stringsAsFactors=F)
names(dfAcousticData_99dB) = c("Year", "Month", "Day", "Time_UTC", "Log1_ping", "Log2_ping", "Lat", "Lon", "BDMIN","BDMAX","OBJECT","CH", "DepthMin_m", "DepthMax_m", "DepthMean_m", "UPINLM", "NASC_UHC","TOTAL") #"OTHER", "NASC_MallotusVilosus","KRILL","TOTAL")
dim(dfAcousticData_99dB)
head(dfAcousticData_99dB, 10)

# Aggregate into boxes
dfAcousticData_99dB$Rect_Lat = round(dfAcousticData_99dB$Lat*20)/20
dfAcousticData_99dB$Rect_Lon = round(dfAcousticData_99dB$Lon*10)/10
dfAcousticData_99dB = unique(dfAcousticData_99dB[,c("Rect_Lat","Rect_Lon","DepthMin_m","DepthMax_m","DepthMean_m")])

# Calculate ping volume
dfAcousticData_99dB$Volume_m3 = (alpha_beam/3)*(dfAcousticData_99dB$DepthMax_m^3-dfAcousticData_99dB$DepthMin_m^3)

# Merge to add number of pings
dfAcousticData_99dB = merge(dfAcousticData_99dB, df_pings_2024, by = c("Rect_Lat", "Rect_Lon"), all.x=TRUE) # Always make sure that necessary variables are included in the bz argument and put all=TRUE, as we want to keep the information of the cells that have no fish
table(dfAcousticData_99dB$n_pings) #check no NA
head(dfAcousticData_99dB)

### Load file containing all detected targets
dfSED <- read.table("Nuuk_Fjord_20cmCodSettings_Phase8_3dB_copy.txt",sep=",",header=T,skip=17)
names(dfSED) = c("Date", "Time", "Lat", "Lon", "Depth_m", "TSC", "TSU", "AlongshipAngle", "AthwartshipAngle","sV_of_peak")
head(dfSED)
dim(dfSED)

dfSED$DepthMin_m = floor(dfSED$Depth_m/10)*10
dfSED$DepthMax_m = dfSED$DepthMin_m+10
dfSED$DepthMean_m = (dfSED$DepthMin_m+dfSED$DepthMax_m)/2

dfSED$Rect_Lat = round(dfSED$Lat*20)/20
dfSED$Rect_Lon = round(dfSED$Lon*10)/10

# Add dummy variable that represents the number of fish (n detected targets)
dfSED$N = 1

# Aggregate N based on cells
dfSED = aggregate(N~Rect_Lat+Rect_Lon+DepthMin_m+DepthMax_m+DepthMean_m, dfSED, sum)

# Merge
dfAcousticData_Counting_merged = merge(dfAcousticData_99dB, dfSED, by = c("Rect_Lat", "Rect_Lon", "DepthMin_m", "DepthMax_m", "DepthMean_m"), all.x=TRUE)
# Always make sure that necessary variables are included in the bz argument and put all=TRUE, as we want to keep the information of the cells that have no fish
dim(dfAcousticData_Counting_merged)
# check the dimensions of the data frame!

# Substitute NA values for the varible N generated when merging with 0
dfAcousticData_Counting_merged$N[is.na(dfAcousticData_Counting_merged$N)] <- 0
head(dfAcousticData_Counting_merged)
sum(dfAcousticData_Counting_merged$N)
# Note: 100 detection lost

# Total N per ping
dfAcousticData_Counting_merged$N_perping = dfAcousticData_Counting_merged$N/dfAcousticData_Counting_merged$n_pings

# Split cod - model derived from separate script
dfAcousticData_Counting_merged$N_cod_perping <- dfAcousticData_Counting_merged$N_perping*(1+tanh(-0.00235*pi*(dfAcousticData_Counting_merged$DepthMean_m-227)))/2
sum(dfAcousticData_Counting_merged$N_cod_perping)

### Derive cod density

# Volume - m3 = m2x10m
dfAcousticData_Counting_merged$Density_N_per_m3 = dfAcousticData_Counting_merged$N_cod_perping/dfAcousticData_Counting_merged$Volume_m3
head(dfAcousticData_Counting_merged)

# Area - m2
dfAcousticData_Counting_merged_Area = dfAcousticData_Counting_merged
dfAcousticData_Counting_merged_Area$DensityA_N_per_m2 = dfAcousticData_Counting_merged_Area$Density_N_per_m3*10
dfAcousticData_Counting_merged_Area = aggregate(DensityA_N_per_m2~Rect_Lat+Rect_Lon, dfAcousticData_Counting_merged_Area, sum)

# Area density in nmi
dfAcousticData_Counting_merged_Area$DensityA_N_per_nmi2 = dfAcousticData_Counting_merged_Area$DensityA_N_per_m2*(1852^2)
head(dfAcousticData_Counting_merged_Area)


##########################################################################################
#####                                 CREATE MAPS                                    ##### 
##########################################################################################

# Define box around Nuuk Fjord
#nuuk_bbox <- st_bbox(c(xmin = min(dfAcousticData_Counting_merged_Area$Rect_Lon), ymin = min(dfAcousticData_Counting_merged_Area$Rect_Lat), xmax = max(dfAcousticData_Counting_merged_Area$Rect_Lon), ymax = max(dfAcousticData_Counting_merged_Area$Rect_Lat)),crs = st_crs(4326))
#nuuk_bbox <- as.array(nuuk_bbox)

# Get the countries from the `rnaturalearth` package
#nuuk_countries <- ne_countries(scale = "large", returnclass = "sf") %>% st_make_valid() %>%
#  st_intersection(st_as_sfc(nuuk_bbox))

# Get the landmasses from the `rnaturalearthdata` package
#landmasses <- ne_download(type = "land", category = "physical", returnclass = "sf")

# Load "Prop_Area_Rect_Sea" file
setwd("C:/Giacomo Gardella/Comparisons LSSS")
true_rect_area = read.csv("Prop_Area_Rect_Sea.csv")

### Area density map
AD_2024_map_PeakPdev8 <- ggplot() +
  geom_tile(data=true_rect_area, aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_Counting_merged_Area, aes(x = Rect_Lon, y = Rect_Lat, fill = DensityA_N_per_nmi2), color="black") +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 210000),  
                      breaks = seq(0, 200000, by = 50000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude",title = "Echo-counting - 2024", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=14016702", size = 4, hjust = 0) + 
  theme_classic()
AD_2024_map_PeakPdev8
#ggsave(filename="new Area D map 2024 Peak Phase8.png", plot=AD_2024_map_PeakPdev8, bg="white", width=5.3, height=4, dpi=500)


##########################################################################################
####                   ESTIMATE THE TOTAL NUMBER OF COD IN THE FJORD                  ####
##########################################################################################

### Simple mean
#mean_AreaDensity_cod_count <- mean(dfAcousticData_Counting_merged_Area$DensityA_N_per_nmi2) # fish per nmi^2


### Weighted average - weights are the proportion of the area of each rectangle that is at sea

# Load "Prop_Area_Rect_Sea" file
setwd("C:/Giacomo Gardella/Comparisons LSSS")
true_rect_area = read.csv("Prop_Area_Rect_Sea.csv")

# Merge
dfAcousticData_Counting_merged_Area = merge(dfAcousticData_Counting_merged_Area, true_rect_area, by=c("Rect_Lat","Rect_Lon"), all.x=TRUE)

# Multiply desnities by corresponding weights
dfAcousticData_Counting_merged_Area_2 = dfAcousticData_Counting_merged_Area[!is.na(dfAcousticData_Counting_merged_Area$Prop_Area_Sea_Norm),]
dfAcousticData_Counting_merged_Area_2$Prop_Area_Sea_Norm = dfAcousticData_Counting_merged_Area_2$Prop_Area_Sea/sum(dfAcousticData_Counting_merged_Area_2$Prop_Area_Sea)
dfAcousticData_Counting_merged_Area_2$Density_N_nmi2_cod_weighted = dfAcousticData_Counting_merged_Area_2$DensityA_N_per_nmi2 *dfAcousticData_Counting_merged_Area_2$Prop_Area_Sea_Norm

# Add to get weighted mean
mean_AreaDensity_cod_count <- sum(dfAcousticData_Counting_merged_Area_2$Density_N_nmi2_cod_weighted) # fish per nmi^2


# Scale to entire fjord
dfBathymetryData = read.csv("Godthåbsfjord_Ameralik_strata_20200414.csv", sep=";")
numArea_nm2 = (sum(dfBathymetryData$areakm2)*(1000^2))/(1852^2)

Cod_total_count <- mean_AreaDensity_cod_count*numArea_nm2
Cod_total_count # weighted 14016702 - if recalculate true_area_rect this value will always be slightly different because the proportions are based on random points
#not weighted 14895393


# Aggregated by layer
dfAcousticData_Counting_merged_Layer = aggregate(Density_N_per_m3~DepthMean_m, dfAcousticData_Counting_merged, mean)

# Compute smoothers
smoother_int_2024 <-sreg(I(rev(dfAcousticData_integration_layer$DepthMean_m)),dfAcousticData_integration_layer$Density_N_per_m3, lambda = 0.2)
predVD_int_2024 <- predict(smoother_int_2024)
smoother_count_2024 <-sreg(I(rev(dfAcousticData_Counting_merged_Layer$DepthMean_m)),dfAcousticData_Counting_merged_Layer$Density_N_per_m3, lambda = 0.2)
predVD_count_2024 <- predict(smoother_count_2024)

##########################################################################################
####                        INVESTIGATE VERTICAL PROFILE                              ####
##########################################################################################

# Vertical profile
plot(dfAcousticData_integration_layer$Density_N_per_m3,dfAcousticData_integration_layer$DepthMean_m, col="white",
     ylim=rev(range(dfAcousticData_integration_layer$DepthMean_m)),
     xlab=expression("Volume density (cod"~m^3~")"), ylab="Depth (m)", xlim=c(0,0.00006))

lines(predVD_int_2024,dfAcousticData_integration_layer$DepthMean_m, lwd=6, col="purple")
lines(predVD_count_2024,dfAcousticData_Counting_merged_Layer$DepthMean_m, lwd=6,col="orange")

legend(x = "bottomright",bty="n",legend = c("Echo-integration - 2024", "Echo-counting - 2024"),lty = c(1, 1), col = c("purple","orange"),lwd = 6, cex=0.8)


