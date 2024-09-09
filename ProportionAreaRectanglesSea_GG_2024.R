##########################################################################################
###########            SCRIPT TO CALCULATE TRUE AREA OF RECTANGLES            ############
##########################################################################################

## Copyright 2024 Technical University of Denmark
## Author: Giacomo Gardella <giaga@aqua.dtu.dk>

# Load relevant packages           
library("ggplot2")
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

##########################################################################################
## In this script the proportion of the area of each rectangle that is at sea is        ##
## estimated by randomly generating 10 000 points within each rectangle and checking    ##
## whether they are at sea or on land by using the sf and naturalearth packages.        ##
##########################################################################################

# Load in data from other scripts that is aggregated to rectangles
head(dfAcousticData_Counting_merged_Area)

# Make new data file for this analysis
true_rect_area <- dfAcousticData_Counting_merged_Area[,c("Rect_Lat","Rect_Lon")]

# Make new rectangle ID
true_rect_area$Rect_ID <- paste(
  "N", sprintf("%05.2f", true_rect_area$Rect_Lat),
  "W", sprintf("%05.2f", abs(true_rect_area$Rect_Lon)), sep="")
head(true_rect_area)

# Define box around Nuuk Fjord
nuuk_bbox <- st_bbox(c(xmin = -52.2, ymin = 63.8, xmax = -49.5, ymax = 65.2),crs = st_crs(4326))
nuuk_bbox <- as.array(nuuk_bbox)

# Get the countries from the `rnaturalearth` package
nuuk_countries <- ne_countries(scale = "large", returnclass = "sf") %>% st_make_valid() %>%
  st_intersection(st_as_sfc(nuuk_bbox))

# Make empty column
true_rect_area$N_in_Sea <- NA
true_rect_area$N_on_Land <- NA
true_rect_area$Prop_Area_Sea <- NA

### Add all other possible rectangles within the nuuk box

# Take limits of nuuk box
lat_range <- seq(63.8, 65.2, by = 0.05)
lon_range <- seq(-52.2, -49.5, by = 0.1)

extended_grid <- expand.grid(Rect_Lat = lat_range, Rect_Lon = lon_range)

# Create Rect_ID for the full grid
extended_grid$Rect_ID <- paste(
  "N", sprintf("%05.2f", extended_grid$Rect_Lat),
  "W", sprintf("%05.2f", abs(extended_grid$Rect_Lon)), sep="")
dim(extended_grid)
head(extended_grid)

# Merge with data
true_rect_area = merge(extended_grid, true_rect_area, by=c("Rect_Lat", "Rect_Lon", "Rect_ID"), all.x=TRUE) 
dim(true_rect_area)


### For loop to get proportion of area at sea

for (i in 1:length(true_rect_area$Rect_ID)){
  
  processed_rect_ID <- true_rect_area$Rect_ID[i]
  
  # Generate random latitudes and longitudes
  random_latitudes <- runif(10000, min = true_rect_area$Rect_Lat[true_rect_area$Rect_ID==processed_rect_ID]-0.025, max = true_rect_area$Rect_Lat[true_rect_area$Rect_ID==processed_rect_ID]+0.025)
  random_longitudes <- runif(10000, min = true_rect_area$Rect_Lon[true_rect_area$Rect_ID==processed_rect_ID]-0.05, max = true_rect_area$Rect_Lon[true_rect_area$Rect_ID==processed_rect_ID]+0.05)
  
  # Combine into a data frame and convert to sf object
  random_points <- data.frame(longitude = random_longitudes, latitude = random_latitudes)
  random_points_sf <- st_as_sf(random_points, coords = c("longitude", "latitude"), crs = st_crs(4326))
  
  
  # Check if the points are inside nuuk_countries
  inside_nuuk_countries <- st_within(random_points_sf, nuuk_countries, sparse = FALSE)
  true_rect_area$N_on_Land[true_rect_area$Rect_ID==processed_rect_ID] <- sum(inside_nuuk_countries)
  
  # Determine points in the sea (not inside nuuk_countries)
  in_sea <- !(inside_nuuk_countries)
  true_rect_area$N_in_Sea[true_rect_area$Rect_ID==processed_rect_ID] = sum(in_sea)
  
  # Percentage of points in sea
  true_rect_area$Prop_Area_Sea[true_rect_area$Rect_ID==processed_rect_ID]<- sum(in_sea)/10000
  
}
head(true_rect_area)

# Standardized proportions
true_rect_area$Prop_Area_Sea_Norm = true_rect_area$Prop_Area_Sea/sum(true_rect_area$Prop_Area_Sea)
head(true_rect_area)


# Remove all rectangles that are at sea
true_rect_area = subset(true_rect_area, Prop_Area_Sea!=0)
dim(true_rect_area)


### Make map
ggplot() +
  geom_tile(data = true_rect_area, aes(x = Rect_Lon, y = Rect_Lat, fill = Prop_Area_Sea)) +
  scale_fill_gradient(low = "lightblue1", high = "blue4") +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", fill="Proportion of
  area in sea")+
  theme_classic()
#ggsave(filename="entire area grid.png", bg="white", width=5.3, height=4, dpi=500)

#write.csv(true_rect_area, "Prop_Area_Rect_Sea.csv")
head(true_rect_area)





