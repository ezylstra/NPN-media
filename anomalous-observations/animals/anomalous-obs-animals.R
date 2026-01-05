# Anomalous observations script - for animals
# Originally developed by A. Rosemartin and T. Crimmins for plants
# Adapted from https://github.com/alyssarosemartin/spring-media/anomalous_obs
# 12 Feb 2025

library(here)
library(rnpn)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(scales) # To make integer axis labels
library(sf)
library(terra)
library(tidyterra)
library(elevatr)
library(leaflet)

# library(rnaturalearth) 
  # Only needed if us_states shapefile not already present in shapefiles folder


# This code identifies anomalous observations of animals in the current year 
# compared to the long term record (2009-XXXX). First observation date in 
# current year is compared to the distribution of dates from prior years. Early 
# (lower than Xth percentile) and outliers (earlier than 1.5x the IQR) are 
# flagged.

# Need individual phenometrics for current year
# Need individual phenometrics for prior years (for current year species)

# Set parameters --------------------------------------------------------------#

# Current year
year <- 2025

# Prior years
prior_years <- 2009:(year - 1)

# Define early observations (below X percentile)
lowerq <- 0.05

# Logical indicating whether to denote Tukey outliers in addition to early 
# observations (those below X percentile) in histograms
outliers <- FALSE

# Set minimum number of years and observations needed in prior years to 
# evaluate whether a current observation is anomalous
min_yrs <- 8
min_obs <- 20

# Set radius (in km) defining region within which to extract prior observations
# for comparison
radius <- 100

# Set elevational buffer (in m) defining band within which to extract prior
# observations for comparison
elev_buffer <- 1000

# Select taxanomic classes (so we get just animals)
species <- npn_species() %>% data.frame()
class_ids <- unique(species$class_id[species$kingdom == "Animalia"])

# Look at phenophases?
phenophases <- npn_phenophases() %>% 
  # Eliminate plant phenophases (classes 1-13)
  filter(pheno_class_id > 13) %>%
  data.frame()
# arrange(phenophases, phenophase_category, phenophase_id)
# Development and Reproduction categories include observations of young/dead animals

# Set mapping parameters ------------------------------------------------------#

# Set date to extract AGDD anomalies (for map). If don't want today, then set
# agdd_date to desired date.
agdd_today <- TRUE
if (!agdd_today) {
  agdd_date <- ymd("2025-02-03")
} else {
  agdd_date <- today()
}

# Define breaks that will be used to delineate areas cooler/warmer than 30-year
# normals. Breaks = AGDD anomolies on agdd_date.
agdd_breaks <- c(-2000, -40, 40, 2000)

# Colors for map (blue = cooler, offwhite = normal, orange = warmer)
cw_cols <- c("#6fa8d6", "#f7f0da", "#fa824d")
# plot(1:3, rep(1, 3), pch = 19, cex = 5, col = cols)

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

# Aquire, format phenometric data in current year -----------------------------#

# current_dl <- npn_download_individual_phenometrics(
#   request_source = 'erinz', 
#   years = year,
#   class_ids = class_ids,
#   additional_fields = c("observed_status_conflict_flag",
#                         "species_functional_type",
#                         "phenophase_category", 
#                         "pheno_class_id")
# )
# 
# # Simplify data (remove unnecessary columns; replace -9999 with NA)
# current <- current_dl %>%
#   select(site_id, latitude, longitude, elevation_in_meters, state, 
#          common_name, species_id, species_functional_type, phenophase_category,
#          phenophase_description, pheno_class_id, first_yes_year, first_yes_doy, 
#          numdays_since_prior_no, observed_status_conflict_flag) %>%
#   rename(lat = latitude,
#          lon = longitude, 
#          elev = elevation_in_meters, 
#          func_type = species_functional_type,
#          phenophase = phenophase_description,
#          year = first_yes_year,
#          first_yes = first_yes_doy,
#          prior_no = numdays_since_prior_no,
#          flag = observed_status_conflict_flag) %>%
#   mutate(across(c(elev, prior_no), 
#                 ~ ifelse(. == -9999, NA, .))) %>%
#   mutate(across(c(state, flag), 
#                 ~ ifelse(. == "-9999", NA, .))) %>%
#   mutate(state = ifelse(state == "", NA, state)) %>%
#   data.frame()
# 
# # Remove observations with conflicts or those without a prior no in last 30 days
# #TODO: should we relax the prior no constraint?
# current <- current %>%
#   filter(is.na(flag)) %>%
#   filter(!is.na(prior_no)) %>%
#   filter(prior_no <= 30) %>%
#   select(-flag)
# 
# # Evaluate phenophases to aggregate and/or eliminate some
# count(current, phenophase_category, phenophase)
#   # category = Activity (can probably aggregate all)
#   # category = Reproduction (can probably aggregate all)
#   # category = Development (delete all dead, but maybe keep larvae/young/etc?)
#   # category = Method (include with activity)
# 
# # For now will ignore any non-dead Development phenophases
# current <- current %>%
#   filter(phenophase_category != "Development") %>%
#   # filter(str_detect(phenophase, "Dead|dead", negate = TRUE)) %>%
#   mutate(phenophase_category = ifelse(phenophase_category == "Method", 
#                                       "Activity",
#                                       phenophase_category))
# 
# # Keep only first observation of species at a site (for a pheno-category) in
# # current year
# current <- current %>%
#   arrange(common_name, site_id, phenophase_category, first_yes) %>%
#   distinct(common_name, site_id, phenophase_category, .keep_all = TRUE)
# 
# # Some observations missing state ID. Will fill in using states layer from the 
# # tigris package
# states48 <- state.abb[! state.abb %in% c("AK", "HI")]
# state_fill <- filter(current, is.na(state)) %>%
#   select(site_id, lon, lat) %>%
#   distinct()
# state_fillv <- vect(state_fill, crs = "epsg:4326")
# state_new <- terra::extract(states, state_fillv)
# state_fill <- cbind(state_fill, state_new = state_new$postal)
# current <- current %>%
#   left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
#   mutate(state = ifelse(!is.na(state), state, state_new)) %>%
#   select(-state_new) %>%
#   # Restrict observations to the continental US
#   filter(state %in% states48) 
# 
# # Some observations missing elevation. Will fill in using the elevatr package
# elev_fill <- filter(current, is.na(elev)) %>%
#   distinct(site_id, lat, lon) %>%
#   st_as_sf(coords = c("lon", "lat"), crs = 4326)
# elev_fill <- get_elev_point(locations = elev_fill, src = "epqs")
# elev_fill <- data.frame(elev_fill) %>%
#   mutate(elev_new = round(elevation)) %>%
#   select(site_id, elev_new)
# current <- current %>%
#   left_join(elev_fill, by = "site_id") %>%
#   mutate(elev = ifelse(!is.na(elev), elev, elev_new)) %>%
#   select(-elev_new)

# Aquire, format phenometric data in previous years ---------------------------#

# spp <- unique(current$species_id)
# pheno_class_ids <- unique(current$pheno_class_id)
# prior_dl <- npn_download_individual_phenometrics(
#   request_source = 'erinz', 
#   years = prior_years,
#   species_ids = spp,
#   pheno_class_ids = pheno_class_ids,
#   additional_fields = c("observed_status_conflict_flag",
#                         "species_functional_type",
#                         "phenophase_category")
# )
# 
# # Simplify data (remove unnecessary columns; replace -9999 with NA)
# prior <- prior_dl %>%
#   select(site_id, latitude, longitude, elevation_in_meters, state, 
#          common_name, species_id, species_functional_type, phenophase_category,
#          phenophase_description, first_yes_year, first_yes_doy, 
#          numdays_since_prior_no, observed_status_conflict_flag) %>%
#   rename(lat = latitude,
#          lon = longitude, 
#          elev = elevation_in_meters, 
#          func_type = species_functional_type,
#          phenophase = phenophase_description,
#          year = first_yes_year,
#          first_yes = first_yes_doy,
#          prior_no = numdays_since_prior_no,
#          flag = observed_status_conflict_flag) %>%
#   mutate(across(c(elev, prior_no), 
#                 ~ ifelse(. == -9999, NA, .))) %>%
#   mutate(across(c(state, flag), 
#                 ~ ifelse(. == "-9999", NA, .))) %>%
#   mutate(state = ifelse(state == "", NA, state)) %>%
#   data.frame()
# 
# # Remove observations with conflicts or those without a prior no in last 30 days
# prior <- prior %>%
#   filter(is.na(flag)) %>%
#   filter(!is.na(prior_no)) %>%
#   filter(prior_no <= 30) %>%
#   select(-flag)
# 
# # Modify phenophase categories
# prior <- prior %>%
#   filter(phenophase_category != "Development") %>%
#   # filter(str_detect(phenophase, "Dead|dead", negate = TRUE)) %>%
#   mutate(phenophase_category = ifelse(phenophase_category == "Method", 
#                                       "Activity",
#                                       phenophase_category))
# 
# # Keep only first observation of species at a site (for a pheno-category) in
# # a year
# prior <- prior %>%
#   arrange(common_name, site_id, phenophase_category, year, first_yes) %>%
#   distinct(common_name, site_id, phenophase_category, year, .keep_all = TRUE)
# 
# # Fill in state or elevation where needed
# state_fill <- filter(prior, is.na(state)) %>%
#   select(site_id, lon, lat) %>%
#   distinct()
# state_fillv <- vect(state_fill, crs = "epsg:4326")
# state_new <- terra::extract(states, state_fillv)
# state_fill <- cbind(state_fill, state_new = state_new$postal)
# prior <- prior %>%
#   left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
#   mutate(state = ifelse(!is.na(state), state, state_new)) %>%
#   select(-state_new) %>%
#   # Restrict observations to the continental US
#   filter(state %in% states48) 
# 
# elev_fill <- filter(prior, is.na(elev)) %>%
#   distinct(site_id, lat, lon) %>%
#   st_as_sf(coords = c("lon", "lat"), crs = 4326)
# elev_fill <- get_elev_point(locations = elev_fill, src = "epqs")
# elev_fill <- data.frame(elev_fill) %>%
#   mutate(elev_new = round(elevation)) %>%
#   select(site_id, elev_new)
# prior <- prior %>%
#   left_join(elev_fill, by = "site_id") %>%
#   mutate(elev = ifelse(!is.na(elev), elev, elev_new)) %>%
#   select(-elev_new)

# write.csv(current, 
#           "anomalous-observations/current_animals_20250212.csv", 
#           row.names = FALSE)
# write.csv(prior, 
#           "anomalous-observations/prior_animals_20250212.csv", 
#           row.names = FALSE)

# Compare current year observations with distribution of first observation dates
# for plants within xx km and xx m elevation ----------------------------------#

current <- read.csv(here("anomalous-observations", "current_animals_20250212.csv"))
prior <- read.csv(here("anomalous-observations", "prior_animals_20250212.csv"))

current$n_obs_r <- NA
current$n_site_r <- NA
current$n_yrs_r <- NA

for (i in 1:nrow(current)) {
  
  current1 <- current[i, ]
  
  # Get locations for all observations of that species, phenophase
  locs <- prior %>%
    filter(common_name == current1$common_name, 
           phenophase_category == current1$phenophase_category) %>%
    select(site_id, lat, lon) %>%
    distinct()
  
  # Calculate distance between 2025 site and all others
  if (nrow(locs) == 0) {next}
  dists <- terra::distance(x = as.matrix(current1[,c("lon", "lat")]),
                           y = as.matrix(locs[,c("lon", "lat")]), 
                           lonlat = TRUE,
                           unit = "km")
  locs_r <- locs[which(dists <= radius), ] %>%
    mutate(withinr = 1)
  
  # Extract prior observations within radius and elevational band
  prior_r <- prior %>%
    left_join(select(locs_r, site_id, withinr), by = "site_id") %>%
    filter(withinr == 1) %>%
    filter(common_name == current1$common_name) %>%
    filter(phenophase_category == current1$phenophase_category) %>%
    filter(elev > (current1$elev - elev_buffer) & elev < (current1$elev + elev_buffer))
  
  # Add sample size info (for region defined by loc, elev) to current dataframe
  current$n_obs_r[i] <- nrow(prior_r)
  current$n_site_r[i] <- n_distinct(prior_r$site_id)
  current$n_yrs_r[i] <- n_distinct(prior_r$year)
  
  if (current$n_obs_r[i] < min_obs | current$n_yrs_r[i] < min_yrs) {next}
  
  # Create dataframe with distributional summaries from prior years for 
  # current-year observations that have sufficient data in prior years 
  quants_temp <- data.frame(sppsite = paste0(current1$common_name, "_", 
                                             current1$site_id), 
                            common_name = current1$common_name,
                            func_type = current1$func_type,
                            site_id = current1$site_id,
                            lon = current1$lon,
                            lat = current1$lat, 
                            state = current1$state,
                            phenophase_category = current1$phenophase_category,
                            first_yes = current1$first_yes,
                            min = min(prior_r$first_yes),
                            qearly = quantile(prior_r$first_yes, lowerq),
                            q0.25 = quantile(prior_r$first_yes, 0.25),
                            q0.75 = quantile(prior_r$first_yes, 0.75),
                            IQR = IQR(prior_r$first_yes)) %>%
    mutate(whisker = q0.25 - 1.5 * IQR,
           outlier_threshold = ifelse(whisker < min, min, whisker),
           early = ifelse(first_yes < qearly, 1, 0),
           outlier = ifelse(first_yes <= outlier_threshold, 1, 0))
  
  if (exists("quants_r")) {
    quants_r <- rbind(quants_r, quants_temp)
  } else {
    quants_r <- quants_temp
  }
  
  if (quants_temp$early == 1 | quants_temp$outlier == 1) {
    # Create prior dataframe for ggplot
    prior_temp <- prior_r %>%
      select(common_name, func_type, phenophase_category, site_id, lat, lon, 
             elev, state, year, first_yes) %>%
      mutate(sppsite = quants_temp$sppsite) %>%
      mutate(panel = paste0(sppsite, " (", quants_temp$state, ")"))
      
    if (exists("prior_plot_r")) {
      prior_plot_r <- rbind(prior_plot_r, prior_temp)
    } else {
      prior_plot_r <- prior_temp
    }
  }
}  

# Extract just early/outlier observations
quants_r_eo <- quants_r %>%
  filter(early == 1 | outlier == 1) %>%
  select(sppsite, common_name, func_type, state, site_id, phenophase_category, 
         first_yes, min, qearly, q0.25, q0.75, early, outlier) %>%
  mutate(panel =  paste0(sppsite, " (", state, ")"),
         eo = ifelse(outlier == 1, "Outlier", "Early"),
         eo = factor(eo, levels = c("Early", "Outlier"))) %>%
  arrange(func_type, sppsite, phenophase_category)
# write.table(quants_r_eo, "clipboard", sep = "\t", row.names = FALSE)

# Loop through phenophase category and taxonomic group and create ggplot objects 
# with histograms. Each panel represents distribution associated with one 
# current-year observation
taxa_groups <- unique(quants_r_eo$func_type)
for (tg in taxa_groups) {
  
  quants_r_tg <- filter(quants_r_eo, func_type == tg)
  prior_plot_r_tg <- filter(prior_plot_r, func_type == tg)
  phenophase_categories <- unique(quants_r_tg$phenophase_category)
  
  for (phc in phenophase_categories) {
    
    quants_r_p <- filter(quants_r_tg,  phenophase_category == phc)
    prior_plot_r_p <- filter(prior_plot_r_tg, phenophase_category == phc)
  
    plot_name <- paste0(tg, " ", tolower(phc), 
                        ", first observation of the year")
    
    plot_temp <- prior_plot_r_p %>%
      ggplot(aes(x = first_yes)) +
      geom_histogram(bins = 60, fill = "steelblue3") +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      scale_y_continuous(breaks = scales::breaks_extended(Q = c(1, 5, 2, 4, 3)))
    
    if (outliers) {
      plot_temp <- plot_temp +
        geom_vline(data = quants_r_p,
                   aes(xintercept = first_yes, linetype = eo),
                   color = "red") +
        scale_linetype_manual(values = c("dashed", "solid")) +
        labs(x = "Day of year", y = "Count", title = plot_name,
             linetype = paste0(year, " observations")) +
        theme_bw() +
        theme(legend.position = "bottom")
    } else {
      plot_temp <- plot_temp +
        geom_vline(data = quants_r_p, aes(xintercept = first_yes), 
                   color = "red") +
        labs(x = "Day of year", y = "Count", title = plot_name) +
        theme_bw()
    }
    
    # Save ggplot objects with name = plot_r_tg_phc (eg, plot_r_Bird_Activity)
    assign(paste0("plot_r_", tg, "_", phc), plot_temp)
  }
}

plot_r_Bird_Activity
  # For 2025, TRSW in CA is the only species that isn't present year-round in 
  # the area it was observed. However, there are multiple observations in other
  # years in the days immediately after the 2025 observation.
plot_r_Insect_Activity
  # All 2025 species have observations year-round in the regions where they were
  # observed.

