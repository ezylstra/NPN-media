# Request to find places in the US where spring has arrived earlier than ever
# this year
# 27 Mar 2026

library(rnpn)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(ggplot2)
library(terra)
library(tidyterra)

# Set parameters --------------------------------------------------------------#

# Current year
year <- 2026

# First year to download layers for
year1 <- 1981

# Prior years (for data download)
prior_years <- year1:(year - 1)

# Location to save rasters
rast_folder <- paste0(getwd(), "/rasters")

# Helper calls ----------------------------------------------------------------#

# layers <- npn_get_layer_details()
# filter(layers, grepl("leaf_prism", name)) %>% data.frame()

# Want name = "si-x:average_leaf_prism"

# Get raster layers with annual leaf index values -----------------------------#

ind <- "leaf"
ind_files <- paste0(rast_folder, "/", ind, "/", prior_years, "_si-x_", 
                    ind, "_prism.tif")
  
# Download rasters for prior years (only if they don't already exist)
for (i in 1:length(ind_files)) {
  if (!file.exists(ind_files[i]))
    npn_download_geospatial(paste0("si-x:average_", ind, "_prism"), 
                            date = "2025-01-01",
                            date = paste0(prior_years[i], "01-01"),
                            output_path = ind_files[i])
}

# https://geoserver.usanpn.org/geoserver/wcs?service=WCS&version=2.0.1&request=GetCoverage&coverageId=si-x:average_leaf_prism&SUBSET=time("2025-01-01T00:00:00.000Z")&format=geotiff


# Download rasters for current year. For now, used this (with Jeff's help, to fix resolution issues):
# https://geoserver.usanpn.org/geoserver/wcs?service=WCS&version=2.0.1&request=GetCoverage&coverageId=si-x:average_leaf_ncep&SUBSET=time("2026-03-26T00:00:00.000Z")&format=geotiff&content-disposition=attachment&SCALESIZE=http://www.opengis.net/def/axis/OGC/1/i(1405),http://www.opengis.net/def/axis/OGC/1/j(621)&RESOLUTION=http://www.opengis.net/def/axis/OGC/1/i(0.041666666667),http://www.opengis.net/def/axis/OGC/1/j(0.041666666667)
# Saved it as: 2026_si-x_leaf_prism_4km.tif

# Create multilayer raster ----------------------------------------------------#
rast_files <- list.files(paste0(rast_folder, "/leaf"), full.names = TRUE)

for (i in 1:(length(rast_files) - 1)) {
  rast1 <- rast(rast_files[i])
  if (i == 1) {
    allyrs <- rast1 
  } else {
    allyrs <- c(allyrs, rast1)
  }
}
rm(rast1)

current <- rast("rasters/leaf/2026_si-x_leaf_prism_4km.tif")
# Slightly different extent, resolution
current <- resample(current, allyrs)

allyrs26 <- c(allyrs, current)

# Was going to use the which.min function with terra::app, but that doesn't
# properly deal with ties...

# Create new layer with minimum values:
mins <- terra::app(allyrs26, min, na.rm = TRUE)
# Check if the value in 2026 is equal to the minimum
min26 <- (allyrs26[[nlyr(allyrs26)]] == mins) * 1L
# Propagate NAs from 2026
min26 <- terra::mask(min26, allyrs26[[nlyr(allyrs26)]])
plot(min26)
names(min26) <- "earliest2026"

# Check a few locations
# df <- data.frame(
#   lon = c(-92.424018, -96.911713, -100.165552),
#   lat = c(30.135268, 28.249617, 40.695847)
# )
# check <- terra::extract(allyrs26, y = df) %>% t() %>% data.frame
# ggplot(check[2:46,], aes(x = X3)) +
#   geom_histogram(binwidth = 1) +
#   geom_vline(xintercept = check[47, 3], color = "blue")

# Load layers with state/county boundaries ------------------------------------#

state48 <- state.abb[!state.abb %in% c("HI", "AK")]

states <- vect("shapefiles/us_states.shp")
states <- terra::subset(states, states$postal %in% state48)
states <- terra::project(states, crs(min26))

counties <- vect("shapefiles/us_counties.shp")
counties <- terra::subset(counties, counties$STUSPS %in% state48)
counties <- terra::project(counties, crs(min26))

# Summarize -------------------------------------------------------------------#

min26f <- as.factor(min26)

map <- ggplot() +
  geom_spatraster(data = min26f, aes(fill = earliest2026)) +
  scale_fill_manual(values = c("lightblue", "salmon3"), 
                    na.value = "transparent",
                    na.translate = FALSE,
                    name = "2026 Earliest?",
                    labels = c("No", "Yes")) + 
  geom_spatvector(data = counties, color = "gray90", fill = NA) +
  geom_spatvector(data = states, color = "gray30", fill = NA) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("spring-index/early-20260326.png",
       map, 
       width = 6.5,
       height = 5,
       units= "in")
