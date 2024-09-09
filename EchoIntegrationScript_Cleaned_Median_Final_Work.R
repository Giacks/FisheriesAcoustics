##########################################################################################
####################                                                  ####################
####################   ESTIMATING FISH DENSITY FROM ECHO-INTEGRATION  ####################
####################              (Nuuk fjord system)                 ####################
##########################################################################################

## Copyright 2024 Technical University of Denmark
## Author: Giacomo Gardella <giaga@aqua.dtu.dk>


##########################################################################################
## This R script can be used to estimate the fish density in boxes of 0.05 deg latitude ##
## 0.1 deg longitude and 10 m depth using the method of echo-integration. Data needed   ##
## is a report 16 file (Sv data) from LSSS exported with 0.1 nmi horizontal grid size   ##
## and 10 m of vertical grid size and with a minimum Sv of -99 dB. Then, report 16 files##
## should be exported from LSSS with a different Sv threshold for each 10 m depth layer ##
## (this is because LSSS does not support a depth-dependent Sv threshold), as explained ##
## by Rudstam et al. 2009. In addition a "TS distribution" file should be exported from ##
##LSSS, containing the results of the SED algorithm.                                    ##
##########################################################################################


##### Getting started

# Load relevant packages           
library("ggplot2")
library('ggpubr')
library(dplyr)
library(jpeg)
library(png)
library(fields)
library(devtools)
library(readxl)
library(sp)
library(mgcv)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(VGAM)
library(raster)
library(rnaturalearthhires)

#Set working directory
setwd("C:/Giacomo Gardella/Comparisons LSSS")

### Load in data file with Sv - TSU values, i.e. corresponding Sv threshold of a -45 dB target for each depth layer
sv_TSumin_lsss <- as.data.frame(read_excel("C:/Giacomo Gardella/Comparisons LSSS/sv_thresholds_2024_FMnight.xlsx", sheet=1))
head(sv_TSumin_lsss)
dim(sv_TSumin_lsss)

# Add depth layer information
sv_TSumin_lsss$DepthMin_m = floor(sv_TSumin_lsss$depth/10)*10
sv_TSumin_lsss$DepthMax_m = sv_TSumin_lsss$DepthMin_m+10
sv_TSumin_lsss$DepthMean_m = (sv_TSumin_lsss$DepthMin_m+sv_TSumin_lsss$DepthMax_m)/2

# Aggregate based on layer
sv_summarized_layer <- sv_TSumin_lsss %>%
  dplyr::group_by(DepthMean_m) %>%
  dplyr::summarize(Sv_mean_round = round(mean(Sv)),
            Sv_mean = mean(Sv),
            Sv_sd = sd(Sv),
            TSU_mean_round = round(mean(TSU)),
            TSU_mean = mean(TSU),
            TSU_sd = sd(TSU),
            depth_average_round = round(mean(depth)),
            depth_average = mean(depth),
            depth_sd = sd(depth),
            n = length(TSU))
sv_summarized_layer <- as.data.frame(sv_summarized_layer)
head(sv_summarized_layer)

### Load report 16 of entire fjord - -99 db threshold: Nuuk_Fjord_2024FMNight_01nmi10m_99dB; Nuuk_Fjord_99db_2.txt
dfAcousticData_99dB = read.table("Nuuk_Fjord_2024FMNight_01nmi10m_99dB_2.txt", quote = "", header=T, sep="", dec=".", blank.lines.skip = FALSE, stringsAsFactors=F)
names(dfAcousticData_99dB) = c("Year", "Month", "Day", "Time_UTC", "Log1_ping", "Log2_ping", "Lat", "Lon", "BDMIN","BDMAX","OBJECT","CH", "DepthMin_m", "DepthMax_m", "DepthMean_m", "UPINLM", "NASC_UHC","TOTAL") #"OTHER", "NASC_MallotusVilosus","KRILL","TOTAL")
dim(dfAcousticData_99dB)

# Add depth and rectangle information
dfAcousticData_99dB$Rect_Lat = round(dfAcousticData_99dB$Lat*20)/20
dfAcousticData_99dB$Rect_Lon = round(dfAcousticData_99dB$Lon*10)/10
dfAcousticData_99dB$Samples = 1

# Aggregate to get samples for each box, i.e. how many small cells are available and thus how much data
dfAcousticData_99dB = aggregate(Samples~Rect_Lat+Rect_Lon+DepthMin_m+DepthMax_m+DepthMean_m, dfAcousticData_99dB, sum)
dim(dfAcousticData_99dB)
head(dfAcousticData_99dB)

##########################################################################################
####                                LOAD SV FILES                                     ####
##########################################################################################

# Change working directory. List all CSV files in the directory
setwd("C:/Giacomo Gardella/Comparisons LSSS/LSSS Sv thresholds Nuuk 2024")
file_list_Nuuk2024_sv_lsss <- list.files(pattern = "*.txt")

# Create empty data frame
dfAcousticData_Empty <- data.frame(Rect_Lat=0, Rect_Lon=0, DepthMin_m=0, DepthMax_m=0, DepthMean_m=0, NASC_UHC=0, Sv_threshold=0)

# Read each .txt file and assign it to a variable with the same name as the file if it needs to be stored
for (i in seq_along(file_list_Nuuk2024_sv_lsss)) {
  
  # Extract name of file
  file_name <- gsub(".txt", "", file_list_Nuuk2024_sv_lsss[i])
  
  # Read file
  intdata = read.table(file_list_Nuuk2024_sv_lsss[i], quote = "", header=T, sep="", dec=".", blank.lines.skip = FALSE, stringsAsFactors=F)
  names(intdata) = c("Year", "Month", "Day", "Time_UTC", "Log1_ping", "Log2_ping", "Lat", "Lon", "BDMIN","BDMAX","OBJECT","CH", "DepthMin_m", "DepthMax_m", "DepthMean_m", "UPINLM", "NASC_UHC","TOTAL")
  
  intdata$Rect_Lat = round(intdata$Lat*20)/20
  intdata$Rect_Lon = round(intdata$Lon*10)/10
  intdata = aggregate(NASC_UHC~Rect_Lat+Rect_Lon+DepthMin_m+DepthMax_m+DepthMean_m, intdata, mean)
  dim(intdata)
  head(intdata)
  
  # Add Sv threshold variable to intdata
  intdata$Sv_threshold = 0-as.numeric(substr(file_name,21,22))
  
  # Find depth mean values for given threshold
  sv_values_layers <- subset(sv_summarized_layer, Sv_mean_round==intdata$Sv_threshold[1])
  intdata <- subset(intdata, DepthMean_m %in% unique(sv_values_layers$DepthMean_m))
  
  #assign(file_name, intdata)
  
  dfAcousticData_Empty <- rbind(dfAcousticData_Empty, intdata)
}
dim(dfAcousticData_Empty)
head(dfAcousticData_Empty)

# Clean and rename data frame
dfAcousticData_Integration_Combined = dfAcousticData_Empty
dfAcousticData_Integration_Combined = subset(dfAcousticData_Integration_Combined, Rect_Lat!=0)
dim(dfAcousticData_Integration_Combined)
head(dfAcousticData_Integration_Combined)


# Merge
dfAcousticData_Integration_merged = merge(dfAcousticData_99dB, dfAcousticData_Integration_Combined, by = c("Rect_Lat","Rect_Lon","DepthMin_m","DepthMax_m","DepthMean_m"), all=TRUE) # Always make sure that necessary variables are included in the by argument and put all=TRUE, as we want to keep the information of the cells that have no fish
dfAcousticData_Integration_merged$NASC_UHC[is.na(dfAcousticData_Integration_merged$NASC_UHC)] <- 0
dim(dfAcousticData_Integration_merged)
head(dfAcousticData_Integration_merged)

#Set working directory
setwd("C:/Giacomo Gardella/Comparisons LSSS")

# Load TS distribution file
dfSED <- read.table("Nuuk_Fjord_20cmCodSettings_Phase8_MaxGain3_copy.txt",sep=",",header=T,skip=17)
names(dfSED) = c("Date", "Time", "Lat", "Lon", "Depth_m", "TSC", "TSU", "AlongshipAngle", "AthwartshipAngle","sV_of_peak")
head(dfSED)

dfSED$DepthMin_m = floor(dfSED$Depth_m/10)*10
dfSED$DepthMax_m = dfSED$DepthMin_m+10
dfSED$DepthMean_m = (dfSED$DepthMin_m+dfSED$DepthMax_m)/2

# Calculate the backscattering cross-section corresponding to each TS
dfSED$Sigma_bs <- 10^(dfSED$TSC/10)

dfSED$Rect_Lat = round(dfSED$Lat*20)/20
dfSED$Rect_Lon = round(dfSED$Lon*10)/10
dfSED$Samples = 1

# Aggregate to take mean of sigma_bs and the total number of samples
dfSED_Mean = aggregate(Sigma_bs~Rect_Lat+Rect_Lon+DepthMin_m+DepthMax_m+DepthMean_m, dfSED, mean) # Mean sigma_bs
dfSED_Samples = aggregate(Samples~Rect_Lat+Rect_Lon+DepthMin_m+DepthMax_m+DepthMean_m, dfSED, sum) # Number of samples

# Merge
dfSED = merge(dfSED_Mean, dfSED_Samples, by=c("Rect_Lat","Rect_Lon","DepthMin_m","DepthMax_m","DepthMean_m"), all=TRUE)

# Discard mean sigma_bs derived from less than 20 targets
dfSED$Sigma_bs[dfSED$Samples<20] <- NA
head(dfSED)
dim(dfSED)

# Merge TS distribution and NASC data frames
dfAcousticData_Integration_merged = merge(dfAcousticData_Integration_merged, dfSED, by = c("Rect_Lat","Rect_Lon","DepthMin_m","DepthMax_m","DepthMean_m"), all.x=TRUE)
dfAcousticData_Integration_merged$NASC_UHC[dfAcousticData_Integration_merged$Samples.x<5] <- NA
dim(dfAcousticData_Integration_merged)

dfAcousticData_Integration_merged <- subset(dfAcousticData_Integration_merged, Samples.x>=5)


##########################################################################################
####                                INTERPOLATION                                     ####
##########################################################################################


# Rename file
dfAcousticData_Integration_int <- dfAcousticData_Integration_merged
dim(dfAcousticData_Integration_int)


# Holes are filled by taking the mean sigma_bs of the same fjord area and 10 m depth interval.
# The reason for this is that TS changes more with depth than it does with fjord area (see spatial patterns of TS script).

### Add fjord Area to dfAcousticData_int

# Polygons were made in Google Earth and need to be read in this specific order:
# IceFjord_Inner -> Kapisillit -> IceFjord_outer -> Qoqqut -> ParallelIslands -> Ameralik_Inner -> Ameralik_outer -> Entrance

dfAcousticData_Integration_int$Area <- "NA"

dfAcousticData_Integration_int$Area[dfAcousticData_Integration_int$Area=="NA" & point.in.polygon(point.x=dfAcousticData_Integration_int$Rect_Lon,
                             point.y=dfAcousticData_Integration_int$Rect_Lat,
                             pol.x=c(-50.95070228685925, -50.75698716725122, -50.37229263961623, -50.07940324745356, -49.78160492154712, -49.32447220698425, -49.45884964929406, -50.95070228685925),
                             pol.y=c(64.8758824100782, 64.60580369850079, 64.58796572761753, 64.37812440557245, 64.26825365591989, 64.25687926550468, 64.90458172862124, 64.8758824100782))==1]="IceFjord_Inner"

dfAcousticData_Integration_int$Area[dfAcousticData_Integration_int$Area=="NA" & point.in.polygon(point.x= dfAcousticData_Integration_int$Rect_Lon,
                                                 point.y=dfAcousticData_Integration_int$Rect_Lat,
                                                 pol.x=c(-50.55914325503232, -50.66125709309799, -50.72315988150282, -50.6692232043405, -49.88819352279002, -49.32447220698425, -49.45884964929406, -50.55914325503232),
                                                 pol.y=c(64.46513783208208, 64.47965135451753, 64.45352158499692, 64.40743071416496, 64.26275972317113, 64.25687926550468, 64.90458172862124, 64.46513783208208))==1]="Kapisillit"

dfAcousticData_Integration_int$Area[dfAcousticData_Integration_int$Area=="NA" & point.in.polygon(point.x= dfAcousticData_Integration_int$Rect_Lon,
                                                 point.y=dfAcousticData_Integration_int$Rect_Lat,
                                                 pol.x=c(-51.65007092266644, -51.49362974864551, -50.87570887016447, -50.78365687565665, -49.67668456823876, -49.32447220698425, -49.45884964929406, -51.65007092266644),
                                                 pol.y=c(64.81296920414044, 64.54638372114822, 64.52638946360194, 64.60866138190788, 64.65817728213058, 64.25687926550468, 64.90458172862124, 64.81296920414044))==1]="IceFjord_Outer"

dfAcousticData_Integration_int$Area[dfAcousticData_Integration_int$Area=="NA" & point.in.polygon(point.x= dfAcousticData_Integration_int$Rect_Lon,
                                                 point.y=dfAcousticData_Integration_int$Rect_Lat,
                                                 pol.x=c(-51.23886996016934, -51.4625304806928, -51.2759694914723, -50.79784254136037, -50.39396150867955, -49.32447220698425, -49.45884964929406, -50.81474878035305, -51.10821857357152, -51.23886996016934),
                                                 pol.y=c(64.22800787807211, 64.20277095384839, 64.12849690504177, 64.21270406861619, 64.38183856832522, 64.25687926550468, 64.90458172862124, 64.64188000405915, 64.33054310023292, 64.22800787807211))==1]="Qoqqut"

dfAcousticData_Integration_int$Area[dfAcousticData_Integration_int$Area=="NA" & point.in.polygon(point.x= dfAcousticData_Integration_int$Rect_Lon,
                                                 point.y=dfAcousticData_Integration_int$Rect_Lat,
                                                 pol.x=c(-51.64039568325947, -51.77102097768694, -51.24679553435306, -51.18772381473647, -50.93365267578322, -49.32447220698425, -49.45884964929406, -51.64039568325947),
                                                 pol.y=c(64.80961490236811, 64.26757475463523, 64.21838473823684, 64.24945505298705, 64.49926095524991, 64.25687926550468, 64.90458172862124, 64.80961490236811))==1]="ParallelIslands"

dfAcousticData_Integration_int$Area[dfAcousticData_Integration_int$Area=="NA" & point.in.polygon(point.x= dfAcousticData_Integration_int$Rect_Lon,
                                                 point.y=dfAcousticData_Integration_int$Rect_Lat,
                                                 pol.x=c(-50.91053247317214, -51.10045418910667, -51.071982638489, -51.02161413014755, -50.22579051605558, -49.32447220698425, -49.45884964929406, -50.91053247317214),
                                                 pol.y=c(64.20458832721056, 64.17456308319618, 64.13971865911087, 64.05011239038087, 64.06137993089069, 64.25687926550468, 64.90458172862124, 64.20458832721056))==1]="Ameralik_Inner"

dfAcousticData_Integration_int$Area[dfAcousticData_Integration_int$Area=="NA" & point.in.polygon(point.x= dfAcousticData_Integration_int$Rect_Lon,
                                                 point.y=dfAcousticData_Integration_int$Rect_Lat,
                                                 pol.x=c(-50.91053247317214, -51.11549459087487, -51.57199974485408, -51.60064320824927, -51.22988782260224, -49.32447220698425, -49.45884964929406, -50.91053247317214),
                                                 pol.y=c(64.20458832721056, 64.15144395544114, 64.08536575921345, 64.02725986986758, 63.97785036825515, 64.25687926550468, 64.90458172862124, 64.20458832721056))==1]="Ameralik_Outer"

dfAcousticData_Integration_int$Area[dfAcousticData_Integration_int$Area=="NA" & point.in.polygon(point.x= dfAcousticData_Integration_int$Rect_Lon,
                                                 point.y=dfAcousticData_Integration_int$Rect_Lat,
                                                 pol.x=c(-52.13099306060508, -52.04364934352697, -51.80002886474912, -51.34445078832292, -51.22988782260224, -49.32447220698425, -49.45884964929406, -52.13099306060508),
                                                 pol.y=c(64.27509388766522, 63.95076361403334, 63.87495103972306, 63.91977487038407, 63.97785036825515, 64.25687926550468, 64.90458172862124, 64.27509388766522))==1]="Entrance"
dfAcousticData_Integration_int$Area_F <- as.factor(dfAcousticData_Integration_int$Area)

dim(dfAcousticData_Integration_int)

### Load TS distribution file and repeat previous step
dfSED_int <- read.table("Nuuk_Fjord_20cmCodSettings_Phase8_MaxGain3_copy.txt",sep=",",header=T,skip=17)
names(dfSED_int) = c("Date", "Time", "Lat", "Lon", "Depth_m", "TSC", "TSU", "AlongshipAngle", "AthwartshipAngle","sV_of_peak")
head(dfSED_int)

dfSED_int$DepthMin_m = floor(dfSED_int$Depth_m/10)*10
dfSED_int$DepthMax_m = dfSED_int$DepthMin_m+10
dfSED_int$DepthMean_m = (dfSED_int$DepthMin_m+dfSED_int$DepthMax_m)/2
dfSED_int$Sigma_bs <- 10^(dfSED_int$TSC/10)

dfSED_int$Area <- "NA"

dfSED_int$Area[dfSED_int$Area=="NA" & point.in.polygon(point.x=dfSED_int$Lon,
                             point.y=dfSED_int$Lat,
                             pol.x=c(-50.95070228685925, -50.75698716725122, -50.37229263961623, -50.07940324745356, -49.78160492154712, -49.32447220698425, -49.45884964929406, -50.95070228685925),
                             pol.y=c(64.8758824100782, 64.60580369850079, 64.58796572761753, 64.37812440557245, 64.26825365591989, 64.25687926550468, 64.90458172862124, 64.8758824100782))==1]="IceFjord_Inner"

dfSED_int$Area[dfSED_int$Area=="NA" & point.in.polygon(point.x= dfSED_int$Lon,
                                                 point.y=dfSED_int$Lat,
                                                 pol.x=c(-50.55914325503232, -50.66125709309799, -50.72315988150282, -50.6692232043405, -49.88819352279002, -49.32447220698425, -49.45884964929406, -50.55914325503232),
                                                 pol.y=c(64.46513783208208, 64.47965135451753, 64.45352158499692, 64.40743071416496, 64.26275972317113, 64.25687926550468, 64.90458172862124, 64.46513783208208))==1]="Kapisillit"

dfSED_int$Area[dfSED_int$Area=="NA" & point.in.polygon(point.x= dfSED_int$Lon,
                                                 point.y=dfSED_int$Lat,
                                                 pol.x=c(-51.65007092266644, -51.49362974864551, -50.87570887016447, -50.78365687565665, -49.67668456823876, -49.32447220698425, -49.45884964929406, -51.65007092266644),
                                                 pol.y=c(64.81296920414044, 64.54638372114822, 64.52638946360194, 64.60866138190788, 64.65817728213058, 64.25687926550468, 64.90458172862124, 64.81296920414044))==1]="IceFjord_Outer"

dfSED_int$Area[dfSED_int$Area=="NA" & point.in.polygon(point.x= dfSED_int$Lon,
                                                 point.y=dfSED_int$Lat,
                                                 pol.x=c(-51.23886996016934, -51.4625304806928, -51.2759694914723, -50.79784254136037, -50.39396150867955, -49.32447220698425, -49.45884964929406, -50.81474878035305, -51.10821857357152, -51.23886996016934),
                                                 pol.y=c(64.22800787807211, 64.20277095384839, 64.12849690504177, 64.21270406861619, 64.38183856832522, 64.25687926550468, 64.90458172862124, 64.64188000405915, 64.33054310023292, 64.22800787807211))==1]="Qoqqut"

dfSED_int$Area[dfSED_int$Area=="NA" & point.in.polygon(point.x= dfSED_int$Lon,
                                                 point.y=dfSED_int$Lat,
                                                 pol.x=c(-51.64039568325947, -51.77102097768694, -51.24679553435306, -51.18772381473647, -50.93365267578322, -49.32447220698425, -49.45884964929406, -51.64039568325947),
                                                 pol.y=c(64.80961490236811, 64.26757475463523, 64.21838473823684, 64.24945505298705, 64.49926095524991, 64.25687926550468, 64.90458172862124, 64.80961490236811))==1]="ParallelIslands"

dfSED_int$Area[dfSED_int$Area=="NA" & point.in.polygon(point.x= dfSED_int$Lon,
                                                 point.y=dfSED_int$Lat,
                                                 pol.x=c(-50.91053247317214, -51.10045418910667, -51.071982638489, -51.02161413014755, -50.22579051605558, -49.32447220698425, -49.45884964929406, -50.91053247317214),
                                                 pol.y=c(64.20458832721056, 64.17456308319618, 64.13971865911087, 64.05011239038087, 64.06137993089069, 64.25687926550468, 64.90458172862124, 64.20458832721056))==1]="Ameralik_Inner"

dfSED_int$Area[dfSED_int$Area=="NA" & point.in.polygon(point.x= dfSED_int$Lon,
                                                 point.y=dfSED_int$Lat,
                                                 pol.x=c(-50.91053247317214, -51.11549459087487, -51.57199974485408, -51.60064320824927, -51.22988782260224, -49.32447220698425, -49.45884964929406, -50.91053247317214),
                                                 pol.y=c(64.20458832721056, 64.15144395544114, 64.08536575921345, 64.02725986986758, 63.97785036825515, 64.25687926550468, 64.90458172862124, 64.20458832721056))==1]="Ameralik_Outer"

dfSED_int$Area[dfSED_int$Area=="NA" & point.in.polygon(point.x= dfSED_int$Lon,
                                                 point.y=dfSED_int$Lat,
                                                 pol.x=c(-52.13099306060508, -52.04364934352697, -51.80002886474912, -51.34445078832292, -51.22988782260224, -49.32447220698425, -49.45884964929406, -52.13099306060508),
                                                 pol.y=c(64.27509388766522, 63.95076361403334, 63.87495103972306, 63.91977487038407, 63.97785036825515, 64.25687926550468, 64.90458172862124, 64.27509388766522))==1]="Entrance"
dfSED_int$Area_F <- as.factor(dfSED_int$Area)
dfSED_int$n <- 1
dim(dfSED_int)
Sigma_bs_Area_Depth10m <- aggregate(Sigma_bs~Area_F+DepthMean_m, dfSED_int, mean) # median

# Check what percentage of the NASC_UHC does not have a sigma_bs values

#24% NASC missing!
sum(dfAcousticData_Integration_int$NASC_UHC[is.na(dfAcousticData_Integration_int$Sigma_bs)&(dfAcousticData_Integration_int$NASC_UHC>0)])/sum(dfAcousticData_Integration_int$NASC_UHC)
dim(dfAcousticData_Integration_int[is.na(dfAcousticData_Integration_int$Sigma_bs)&(dfAcousticData_Integration_int$NASC_UHC>0),])
#1128/3638

# Iterate over indices of cells that do not have a sigma_bs
indices <- which(is.na(dfAcousticData_Integration_int$Sigma_bs) & dfAcousticData_Integration_int$NASC_UHC > 0)

# Add sigma_s from the same area and 10 m depth layer
for (i in indices) {
  area1 <- dfAcousticData_Integration_int$Area[i]
  depth1 <- dfAcousticData_Integration_int$DepthMean_m[i]
  
  # Find matching Sigma_bs value
  new_sigma_bs <- Sigma_bs_Area_Depth10m$Sigma_bs[Sigma_bs_Area_Depth10m$Area == area1 & Sigma_bs_Area_Depth10m$DepthMean_m == depth1]
  
  # Check if new_sigma_bs has a matching value before assignment
  if (length(new_sigma_bs) > 0) {
    dfAcousticData_Integration_int$Sigma_bs[i] <- new_sigma_bs
  }
}
# Check how many cells are left
# 5 cells left

# Fill in last 5 holes with data from adjacent depth layers, i.e. +- 10 m
indices <- which(is.na(dfAcousticData_Integration_int$Sigma_bs) & dfAcousticData_Integration_int$NASC_UHC > 0)
for (i in indices) {
  area1 <- dfAcousticData_Integration_int$Area[i]
  depth1 <- dfAcousticData_Integration_int$DepthMean_m[i]
  
  # Find matching Sigma_bs value
  new_sigma_bs <- mean(Sigma_bs_Area_Depth10m$Sigma_bs[Sigma_bs_Area_Depth10m$Area == area1 & Sigma_bs_Area_Depth10m$DepthMean_m >= depth1-10 & Sigma_bs_Area_Depth10m$DepthMean_m <= depth1+10])
  
  # Check if new_sigma_bs has a matching value before assignment
  if (length(new_sigma_bs) > 0) {
    dfAcousticData_Integration_int$Sigma_bs[i] <- new_sigma_bs
  }
}
# No cells missing

### Derive fish density

# NASC divided by Sigma_sp - #fish/(nmi^2*10m)
dfAcousticData_Integration_int$Density_N_nmi2x10m <- dfAcousticData_Integration_int$NASC_UHC*1.05704/(dfAcousticData_Integration_int$Sigma_bs*4*pi) #exp(sigma^2*log(10)^2/200) -> sd=1.446563 -> 1.05704
dfAcousticData_Integration_int$Density_N_nmi2x10m[dfAcousticData_Integration_int$NASC_UHC==0] <- 0

# Split cod and redfish following model (1+tanh(-k1*pi*(depth-d0)))/2 - model found in different script
dfAcousticData_Integration_int$Density_N_nmi2x10m_cod <- dfAcousticData_Integration_int$Density_N_nmi2x10m*(1+tanh(-0.00235*pi*(dfAcousticData_Integration_int$DepthMean_m-227)))/2

dfAcousticData_Integration_AreaDensity <- aggregate(Density_N_nmi2x10m_cod ~ Rect_Lat + Rect_Lon, dfAcousticData_Integration_int, sum)
names(dfAcousticData_Integration_AreaDensity)[3] = c("Density_N_nmi2_cod")


##########################################################################################
#####                                 CREATE MAPS                                    ##### 
##########################################################################################


# Define box around Nuuk Fjord
#nuuk_bbox <- st_bbox(c(xmin = min(dfAcousticData_Integration_AreaDensity$Rect_Lon)-0.1, ymin = min(dfAcousticData_Integration_AreaDensity$Rect_Lat)-0.025, xmax = max(dfAcousticData_Integration_AreaDensity$Rect_Lon)+0.1, ymax = max(dfAcousticData_Integration_AreaDensity$Rect_Lat)+0.025),crs = st_crs(4326))
#nuuk_bbox <- as.array(nuuk_bbox)

# Get the countries from the `rnaturalearth` package
#nuuk_countries <- ne_countries(scale = "large", returnclass = "sf") %>% st_make_valid() %>%
#  st_intersection(st_as_sfc(nuuk_bbox))

# Get the landmasses from the `rnaturalearthdata` package
#landmasses <- ne_download(type = "land", category = "physical", returnclass = "sf")

#
setwd("C:/Giacomo Gardella/Comparisons LSSS")
true_rect_area = read.csv("Prop_Area_Rect_Sea.csv")

# Plot
abundance_map_FM <- ggplot() +
  geom_tile(data = true_rect_area, aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_Integration_AreaDensity, aes(x = Rect_Lon, y = Rect_Lat, fill = Density_N_nmi2_cod), color="black") +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 220600),  
                      breaks = seq(0, 200000, by = 50000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", title="Echo-integration - 2024", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=13995992", size = 4, hjust = 0) + 
  theme_classic()
abundance_map_FM
#ggsave(filename="new abundance map FM final.png", plot=abundance_map_FM, bg="white", width=5.3, height=4, dpi=500)


##########################################################################################
####                   ESTIMATE THE TOTAL NUMBER OF COD IN THE FJORD                  ####
##########################################################################################

### Simple mean
#mean_AreaDensity_Integration_fishnmi2 <- mean(dfAcousticData_Integration_AreaDensity$Density_N_nmi2_cod)



### Weighted average - weights are the proportion of the area of each rectangle that is at sea

# Load "Prop_Area_Rect_Sea" file
setwd("C:/Giacomo Gardella/Comparisons LSSS")
true_rect_area = read.csv("Prop_Area_Rect_Sea.csv")

# Merge
dfAcousticData_Integration_AreaDensity = merge(dfAcousticData_Integration_AreaDensity, true_rect_area, by=c("Rect_Lat","Rect_Lon"), all.x=TRUE)

# Multiply fish densities by the weights
dfAcousticData_Integration_AreaDensity$Prop_Area_Sea_Norm = dfAcousticData_Integration_AreaDensity$Prop_Area_Sea/sum(dfAcousticData_Integration_AreaDensity$Prop_Area_Sea, na.rm=TRUE)
dfAcousticData_Integration_AreaDensity$Density_N_nmi2_cod_weighted = dfAcousticData_Integration_AreaDensity$Density_N_nmi2_cod*dfAcousticData_Integration_AreaDensity$Prop_Area_Sea_Norm

# Add up weighted densities to get weighte mean
mean_AreaDensity_Integration_fishnmi2 <- sum(dfAcousticData_Integration_AreaDensity$Density_N_nmi2_cod_weighted, na.rm=TRUE) # fish per nmi^2
 
#Scale up to entire fjord
dfBathymetryData = read.csv("Godthåbsfjord_Ameralik_strata_20200414.csv", sep=";")
numArea_nm2 = (sum(dfBathymetryData$areakm2)*(1000^2))/(1852^2)

Cod_total_Integration <- mean_AreaDensity_Integration_fishnmi2*numArea_nm2
Cod_total_Integration #13995992 if weighted with .csv file

# Proportion of WISC stock (approximation)
WISC_total_Integration <- Cod_total_Integration*0.55
WISC_total_Integration

#SAM only Nuuk
(2949+3329+4449+1516+1894+723+241+64+14)*1000

##########################################################################################
####                        INVESTIGATE VERTICAL PROFILE                              ####
##########################################################################################


# Aggregate by layer
dfAcousticData_integration_layer = aggregate(Density_N_nmi2x10m_cod~DepthMean_m, dfAcousticData_Integration_int, mean)
dfAcousticData_integration_layer$Density_N_per_m3 = dfAcousticData_integration_layer$Density_N_nmi2x10m_cod/((1852^2)*10)


dim(dfAcousticData_integration_layer)
head(dfAcousticData_integration_layer)

# Comparison with counting in other script




