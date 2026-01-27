# Groundhog day script
# Adapted from https://github.com/alyssarosemartin/spring-media/blob/main/groundhog/groundhog.R#L1
# 27 Jan 2026

# 2026 notes: 
# We're no longer producing AGDD forecasts for the next 7 days because we're 
# just using PRISM data. So the maps will reflect AGDD values as the end of 
# January.
# I added code to fill the state code when missing, which added ~10% more sites.

# General reminder: like in the past, we haven't restricted observations to 
# those with a prior no. Map includes all first yeses, regardless of whether 
# there was a prior no or not.

library(here)
library(rnpn)
library(dplyr)
library(lubridate)
library(ggplot2)
library(terra)
library(tidyterra)
library(ggimage)

# Set parameters --------------------------------------------------------------#

# Year
year <- 2026

# Define breaks that will be used to delineate areas cooler/warmer than 30-year
# normals. Breaks = AGDD anomolies on Groundhog day (deg F).
breaks <- c(-500, -20, 20, 500)

# Colors for map (blue = cooler, offwhite = normal, orange = warmer)
cols <- c("#6fa8d6", "#f7f0da", "#fa824d")
# plot(1:3, rep(1, 3), pch = 19, cex = 5, col = cols)

# Helper calls ----------------------------------------------------------------#

layers <- npn_get_layer_details()
# layers$abstract

phenoclasses <- npn_pheno_classes() %>% data.frame()
# phenoclasses

# Get states layer ------------------------------------------------------------#

# File location
states_shp <- here("shapefiles", "us_states.shp")

# If file exists, load it. Otherwise download it first.
if (!file.exists(states_shp)) {
  library(rnaturalearth)
  states <- rnaturalearth::ne_states(country = "united states of america", 
                                     returnclass = "sv")
  states <- states[, "postal"]
} else {
  states <- terra::vect(states_shp)
}

# Aquire, format phenometric data ---------------------------------------------#

# Download individual phenometrics (first yeses) so far this year: 
# breaking leaf bud / initial growth (class 1) and open flowers (class 7)
df <- npn_download_individual_phenometrics(
  request_source = 'erinz', 
  years = year,
  pheno_class_ids = c(1, 7)
)

# Some observations missing state ID. Will fill in using states layer from the
# tigris package adn then restrict observations to continental US
states48 <- state.abb[! state.abb %in% c("AK", "HI")]
state_fill <- filter(df, is.na(state)) %>%
  select(site_id, longitude, latitude) %>%
  rename(lon = longitude, lat = latitude) %>%
  distinct()
state_fillv <- vect(state_fill, crs = "epsg:4326")
state_new <- terra::extract(states, state_fillv)
state_fill <- cbind(state_fill, state_new = state_new$postal)
df <- df %>%
  left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
  mutate(state = ifelse(!is.na(state), state, state_new)) %>%
  select(-state_new) %>%
  # Restrict observations to the continental US
  filter(state %in% states48)

# Simplify first observation data
df <- df %>%
  select(site_id, latitude, longitude, elevation_in_meters, state, 
         common_name, individual_id, phenophase_description, pheno_class_id,
         first_yes_doy, numdays_since_prior_no, last_yes_doy, 
         numdays_until_next_no) %>%
  rename(lat = latitude,
         lon = longitude, 
         elev = elevation_in_meters, 
         id = individual_id, 
         phenophase = phenophase_description,
         first_yes = first_yes_doy,
         prior_no = numdays_since_prior_no,
         last_yes = last_yes_doy,
         next_no = numdays_until_next_no) %>%
  data.frame()

# Look at species included
count(df, common_name)
# Delete any 'ohi'a lehua records, if present
df <- df %>%
  filter(common_name != "'ohi'a lehua")

# Delete any observations after groundhog day (in case we're running this 
# after the fact)
groundhog_date <- ymd(paste0(year, "-02-02"))
groundhog_doy <- yday(groundhog_date)
df <- df %>%
  filter(first_yes <= groundhog_doy)

# Distribution of first observation dates
ggplot(df) +
  geom_histogram(aes(x = first_yes)) +
  facet_wrap(~pheno_class_id, ncol = 1)

# How many observations and individuals?
df %>%
  group_by(pheno_class_id) %>%
  summarize(n_obs = n(),
            n_indiv = n_distinct(id))

# Keep just one "first yes" for each plant and phenophase class 
df <- df %>%
  arrange(pheno_class_id, common_name, id, first_yes) %>%
  distinct(pheno_class_id, id, .keep_all = TRUE)

# How many observations, and what proportion have prior no?
df %>%
  group_by(pheno_class_id) %>%
  summarize(n_obs = n(),
            with_prior_no = sum(is.na(prior_no))) %>%
  data.frame() %>%
  mutate(prop_prior_no = with_prior_no / n_obs)

leaf <- subset(df, pheno_class_id == 1)
flower<- subset(df, pheno_class_id == 7)

# Aquire AGDD anomalies -------------------------------------------------------#

# Acquire raster forecast anomaly for Feb 2 (32 base) or most recent date when
# AGDD data are available
layers %>%
  filter(name == "gdd:agdd_anomaly") %>%
  select(name, title, abstract, dimension.name)

today <- ymd(Sys.Date())
agdd_day <- min(c(today, groundhog_date))

groundhog <- npn_download_geospatial(coverage_id = "gdd:agdd_anomaly", 
                                     date = agdd_day)  
# If max value in raster is 0, then there's no data for that date. Go back one
# day
max_value <- terra::global(groundhog, fun = max, na.rm = TRUE)
if (max_value == 0) {
  groundhog <- npn_download_geospatial(coverage_id = "gdd:agdd_anomaly", 
                                       date = agdd_day - 1)  
}

# plot(groundhog)
# hist(groundhog)

# Convert map to WGS84, which is datum that NPN data are in
groundhog <- terra::project(groundhog, "epsg:4326")

# Create map ------------------------------------------------------------------#

# Note: flower and leaf icons already downloaded into resources folder
# (leaf.svg, flower.ico)

# Bin data
ghd_class <- terra::classify(x = groundhog, 
                             rcl = breaks, 
                             include.lowest = TRUE)
freq(ghd_class)
plot(ghd_class)

ghd_plot <- ggplot() +
  geom_spatraster(data = ghd_class, maxcell = Inf) +
  scale_fill_manual(values = cols, na.value = "transparent") +
  geom_image(data = leaf, 
             aes(x = lon, y = lat, image = here("icons", "leaf.png"))) +
  geom_image(data = flower, 
             aes(x = lon, y = lat, image = here("icons", "flower.ico"))) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) 
ghd_plot

ghd_filename <- here("groundhog-day/output",
                     paste0("groundhog-day-", year, ".png"))

# Save to file (commented out for now, so we don't accidentally overwrite anything)
# ggsave(filename = ghd_filename,
#        plot = ghd_plot,
#        width = 6,
#        height = 5,
#        units = "in",
#        dpi = 600,
#        bg = "transparent")
