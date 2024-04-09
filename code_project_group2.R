# Import libraries
library(readxl)
library(dplyr)
library(readr)
easypackages::packages(
    "sf", "sp", "spdep", "spatialreg", "GWmodel", "tmap", "mapview", "car", "RColorBrewer",
    "cowplot", "leafsync", "leaflet.extras2", "mapview", "tidyverse"
)

# Import municipality data
data_landuse <- read_excel("CBS_landuse_2017.xlsx")
data_muni <- read_excel("CBS_data_2017.xlsx")

# Transpose
t_data_landuse <- t(data_landuse)

# Get column names
colnames(t_data_landuse) <- data_landuse$Variable

# Remove rows with missing data
clean_data_muni <- na.omit(data_muni)

# Merge datasets to obtain all CBS data per municipality
data_CBS <- merge(clean_data_muni, t_data_landuse)

# Import NO2 data
data_NO2 <- read_sf("data_NO2_2017.gpkg")

# Merge datasets
data <- merge(data_CBS, data_NO2, by.x = "municipality", by.y = "naam")
# Check which municipalities were not included in the merge
not_in_merge <- anti_join(data_CBS, data_NO2, by = c("municipality" = "naam"))

# Change column names
colnames(data)[which(colnames(data) == "mean_NO2me")] <- "mean_NO2"

# Remove open water variable due to NA's
data <- subset(data, select = -c(open_water))

# Change variable type from character to numeric
data$traffic_area <- as.numeric(data$traffic_area)
data$building_area <- as.numeric(data$building_area)
data$semi_building_area <- as.numeric(data$semi_building_area)
data$recreational_area <- as.numeric(data$recreational_area)
data$agricultural_area <- as.numeric(data$agricultural_area)
data$forest_area <- as.numeric(data$forest_area)
data$inland_water <- as.numeric(data$inland_water)

##### Start modelling

# Define equation
equation <- mean_NO2 ~ tot_population + pop_density + house_density + liquid_manure + solid_manure + private_car + road_length + traffic_area + building_area + semi_building_area + recreational_area + agricultural_area + forest_area + inland_water

## Linear regression
linear_model <- lm(equation, data = data)
summary(linear_model)
AIC(linear_model)

## Check for spatial autocorrelation
data_sf <- st_as_sf(data)

# Creating adjacency matrix,
data_nbq <- poly2nb(data_sf, queen = TRUE) # Queen’s Contiguity neighborhood
data_nbq_w <- nb2listw(data_nbq, style = "W", zero.policy = TRUE) # Queen’s neighborhood weights

# Use Monte Carlo method to bootstrap different polygon distribution.
mc_global <- moran.mc(data_sf$mean_NO2, data_nbq_w, 2999, alternative = "greater") # 2999 simulations

# Plot the  Moran’s I
plot(mc_global)
mc_global

## GWR
data_sp <- as_Spatial(data_sf)

# Adaptive kernel
# Find adaptive kernel using gaussian function
adapt_bandw <- bw.gwr(equation,
    approach = "AIC",
    adaptive = TRUE,
    kernel = "gaussian",
    data = data_sp
)

# Fitting the model with adaptive kernel with optimal bandwidth
agwr_model <- gwr.basic(equation,
    adaptive = TRUE,
    kernel = "gaussian",
    bw = adapt_bandw,
    data = data_sp
)
agwr_model

# Fixed kernel
fixed_bandw <- bw.gwr(equation,
    approach = "AIC",
    adaptive = F,
    kernel = "gaussian",
    data = data_sp
)

fgwr_model <- gwr.basic(equation,
    adaptive = FALSE,
    kernel = "gaussian",
    bw = fixed_bandw,
    data = data_sp
)
fgwr_model

## Produce maps of significant variables

# Get the SDF out of model
fgwr_sf <- st_as_sf(fgwr_model$SDF)

# Produce maps for all significant variables
fgwrM1 <- qtm(fgwr_sf, "pop_density", fill.title = "Population density", legend.outside = TRUE)
fgwrM2 <- qtm(fgwr_sf, "solid_manure", fill.title = "Solid manure", legend.outside = TRUE)
fgwrM3 <- qtm(fgwr_sf, "traffic_area", fill.title = "Traffic", legend.outside = TRUE)
fgwrM4 <- qtm(fgwr_sf, "semi_building_area", fill.title = "Semi built-up", legend.outside = TRUE)
fgwrM5 <- qtm(fgwr_sf, "agricultural_area", fill.title = "Agriculture", legend.outside = TRUE)
fgwrM6 <- qtm(fgwr_sf, "forest_area", fill.title = "Forests", legend.outside = TRUE)
fgwrRrsq <- qtm(fgwr_sf, "Local_R2", fill.title = "Local R2", legend.outside = TRUE)

# Plot each map separately for highest quality
tmap_arrange(fgwrRrsq)
tmap_arrange(fgwrM1)
tmap_arrange(fgwrM2)
tmap_arrange(fgwrM3)
tmap_arrange(fgwrM4)
tmap_arrange(fgwrM5)
tmap_arrange(fgwrM6)
