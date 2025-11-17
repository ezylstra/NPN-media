# Anomalous observations script - dashboard
# 17 November 2025

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
library(pins)

board <- board_connect(
  auth = "manual",
  server = Sys.getenv("CONNECT_SERVER"),
  key = Sys.getenv("CONNECT_API_KEY")
)

# This code identifies anomalous observations of plants in the current year 
# compared to the long term record (2009-XXXX). First observation date in 
# current year is compared to the distribution of dates from prior years. Early 
# (lower than Xth percentile) and outliers (earlier than 1.5x the IQR) are 
# flagged.

rm(list = ls())

source(here("anomalous-observations", "functions.R"))

# Set parameters --------------------------------------------------------------#

# Current year
year <- year(Sys.Date())

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

# Set the number of days prior to the first yes when a "no" must have been 
# reported 
prior_no_mindays <- 30

# Select phenophase classes we're interested in
pheno_class_ids <- c(1, 3, 6, 7)
  # 1: initial shoot/leaf growth
  # 3: leaves/needles
  # 6: flowers/cones
  # 7: open flowers/cones 
# Note: for this script, not going to differentiate among phenophases in 
# each phenophase class. Mostly phenophases differ among species or functional
# type, but there are also some old phenophases that are no longer used.

# Postal codes for lower 48 states
states48 <- state.abb[! state.abb %in% c("AK", "HI")]

# Set mapping parameters ------------------------------------------------------#

# Set date to extract AGDD anomalies (for map). Latest possible date = 30 April. 
# If don't want today, then set agdd_date to desired date.
agdd_today <- TRUE
if (!agdd_today) {
  agdd_date <- ymd("2025-02-03")
} else {
  agdd_date <- min(today(), ymd(paste0(year, "-04-30")))
}

# Define breaks that will be used to delineate areas cooler/warmer than 30-year
# normals. Breaks = AGDD anomolies on agdd_date.
agdd_breaks <- c(-2000, -40, 40, 2000)

# Colors for map (blue = cooler, offwhite = normal, orange = warmer)
cw_cols <- c("#6fa8d6", "#f7f0da", "#fa824d")
# plot(1:3, rep(1, 3), pch = 19, cex = 5, col = cols)

# Load states layer -----------------------------------------------------------#

states_shp <- here("shapefiles", "us_states.shp")
states <- vect(states_shp)

# Aquire, format phenometric data in current year -----------------------------#

# For current year, it's probably worth downloading all data every time because
# the individual phenometrics are summarized data, not just raw observations.
# Also, the current year doesn't take too long....

current_dl <- npn_download_individual_phenometrics(
  request_source = 'erinz',
  years = year,
  period_start = "01-01",
  period_end = "04-30",
  pheno_class_ids = pheno_class_ids,
  additional_fields = c("observed_status_conflict_flag",
                        "species_functional_type")
)

current <- ind_ph_cleanup(current_dl, prior_no_mindays = prior_no_mindays)
current <- state_fill_in(current)
current <- elev_fill_in(current)
rm(current_dl)

# Load, merge phenometric data in previous years ------------------------------#

# Existing dataset example filename: obs_phc1_prior.csv
# Existing dataset example filename: obs_phc1_2009-2024.csv ##### WILL SWITCH

for (pheno_class in pheno_class_ids) {
  # Name of existing file for that phenophase class
  prior_file <- list.files(here("anomalous-observations"), 
                           pattern = paste0("obs-phc", pheno_class, "-2009"),
                           full.names = TRUE)
  
  # Load file
  prior_data <- read.csv(prior_file)
  
  # Merge 
  if (pheno_class == pheno_class_ids[1]) {
    prior <- prior_data
  } else {
    prior <- rbind(prior, prior_data)
  }
}
rm(prior_data)

# Compare current year observations with distribution of first observation dates
# for plants within xx km and xx m elevation ----------------------------------#

current <- read.csv(here("anomalous-observations", "current_20250204.csv"))
prior <- read.csv(here("anomalous-observations", "prior_20250204.csv"))

current$n_obs_r <- NA
current$n_indiv_r <- NA
current$n_yrs_r <- NA

for (i in 1:nrow(current)) {
  
  current1 <- current[i, ]
  
  # Get locations for all observations of that species, phenophase
  locs <- prior %>%
    filter(common_name == current1$common_name, 
           pheno_class_id == current1$pheno_class_id) %>%
    select(site_id, lat, lon) %>%
    distinct()
  
  # Calculate distance between current-year plant and all others
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
    filter(pheno_class_id == current1$pheno_class_id) %>%
    filter(elev > (current1$elev - elev_buffer) & elev < (current1$elev + elev_buffer))
  
  # Add sample size info (for region defined by loc, elev) to current dataframe
  current$n_obs_r[i] <- nrow(prior_r)
  current$n_indiv_r[i] <- n_distinct(prior_r$id)
  current$n_yrs_r[i] <- n_distinct(prior_r$year)
  
  if (current$n_obs_r[i] < min_obs | current$n_yrs_r[i] < min_yrs) {next}
  
  # Create dataframe with distributional summaries from prior years for 
  # current-year observations that have sufficient data in prior years 
  quants_temp <- data.frame(indiv = paste0(current1$common_name, "_", current1$id), 
                            common_name = current1$common_name,
                            site_id = current1$site_id,
                            lon = current1$lon,
                            lat = current1$lat, 
                            state = current1$state,
                            id = current1$id, 
                            pheno_class_id = current1$pheno_class_id,
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
      select(common_name, func_type, pheno_class_id, site_id, lat, lon, elev, 
             state, year, first_yes) %>%
      mutate(indiv = quants_temp$indiv) %>%
      mutate(panel = paste0(indiv, " (", quants_temp$state, ")"))
      
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
  select(indiv, common_name, state, site_id, pheno_class_id, 
         first_yes, min, qearly, q0.25, q0.75, early, outlier) %>%
  mutate(panel =  paste0(indiv, " (", state, ")"),
         eo = ifelse(outlier == 1, "Outlier", "Early"),
         eo = factor(eo, levels = c("Early", "Outlier")))

# Loop through phenophase classes and create ggplot objects with histograms
# Each panel represents distribution associated with one current-year observation
for (phc in pheno_class_ids) {
  
  quants_r_p <- filter(quants_r_eo,  pheno_class_id == phc)
  prior_plot_r_p <- filter(prior_plot_r, pheno_class_id == phc)

  phc_name <- paste0(phenoclasses$name[phc], ", first observation of the year")
  
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
      labs(x = "Day of year", y = "Count", title = phc_name,
           linetype = paste0(year, " observations")) +
      theme_bw() +
      theme(legend.position = "bottom")
  } else {
    plot_temp <- plot_temp +
      geom_vline(data = quants_r_p, aes(xintercept = first_yes), 
                 color = "red") +
      labs(x = "Day of year", y = "Count", title = phc_name) +
      theme_bw()
  }
  
  # Save ggplot objects with name = plot_r_Phenophase class (eg, plot_r_1)
  assign(paste0("plot_r_", phc), plot_temp)
}

# Print ggplot object for each phenophase class
# for (phc in pheno_class_ids) {
#   print(get(paste0("plot_r_", phc)))
# }

# Create map with early observations

# Acquire raster with forecasted AGDD anomalies (in F, with 32-deg base)
agdd_anom <- npn_download_geospatial(coverage_id = "gdd:agdd_anomaly", 
                                     date = agdd_date)  

# Prepare binned color palette
pal <- colorBin (palette = cw_cols, 
                 domain = agdd_anom, 
                 bins = agdd_breaks,
                 na.color = "transparent")

# Set background color to white (skipping this for now)
  # backg <- htmltools::tags$style(".leaflet-container { background: white; }" )
# To implement this, would need to add the following to leaflet call:
  # htmlwidgets::prependContent(backg)

# Create flower icon
iflower <- makeIcon(iconUrl = here("icons", "flower.ico"),
                    iconWidth = 15, iconHeight = 15)
# Create leaf icon
ileaf <- makeIcon(iconUrl = here("icons", "leaf.png"),
                  iconWidth = 20, iconHeight = 20)

# Create sf polygon feature with 48 states
states48v <- terra::subset(states, !states$postal %in% c("HI", "AK"))
states48v <- st_as_sf(states48v)

# Radius for comparison observations, in meters
radius_m <- radius * 1000

# Extract early observations for each phenophase class in dataset
for (i in pheno_class_ids) {
  assign(paste0("early", i),
         quants_r_eo %>%
           filter(early == 1 & pheno_class_id == i) %>%
           left_join(select(current, site_id, lon, lat) %>% distinct, 
                     by = c("site_id")))
}
# Now we can have dataframes named: early1, early3, early6, early 7

# Short names for each phenophase class
class1 <- "initial vegetative growth"
class3 <- "leaves or needles"
class6 <- "flowers or cones"
class7 <- "open flowers or cones"

# Create dataframe with mapping parameters for each phenophase class
class_df <- data.frame(
  id = pheno_class_ids,
  title = str_to_sentence(unlist(mget(paste0("class", pheno_class_ids)))))
class_df$col <- ifelse(class_df$id %in% c(1, 3), "green", "purple")
class_df$icon <- ifelse(class_df$id %in% c(1, 3), "ileaf", "iflower")

# Map observations that are earlier than 95% of all previous first-of-the-year 
# observations (and show radius used to delineate area for comparison)
map <- leaflet(states48v) %>%
  addRasterImage(agdd_anom, colors = pal, opacity = 0.6) %>%
  addPolygons(data = states48v, color = "gray", weight = 1, fill = FALSE) %>%
  addLegend("bottomright", pal = pal, values = ~agdd_breaks,
            title = "AGDD anomaly", opacity = 1,
            labFormat = labelFormat(prefix = c("Cool (", "Normal (", "Warm ("),
                                    suffix = ")"))
for (i in 1:nrow(class_df)) {
  early_sub <- get(paste0("early", class_df$id[i]))
  map <- map %>%
    addMarkers(data = early_sub, lng = ~lon, lat = ~lat, 
               group = class_df$title[i], icon = get(class_df$icon[i]),
               popup = ~paste0(common_name, "<br>", 
                               "ID: ", str_split_fixed(indiv, "_", 2))) %>%
    addCircles(data = early_sub, lng = ~lon, lat = ~lat, radius = radius_m,
               color = class_df$col[i], weight = 1, fill = FALSE,
               group = class_df$title[i])
}
map <- map %>%
  addLayersControl(overlayGroups = class_df$title[class_df$id %in% pheno_class_ids],
                   options = layersControlOptions(collapsed = FALSE), 
                   position = "bottomleft")
map

# Create table with early and outlier observations for all phenophase classes
current_join <- current %>%
  mutate(indiv = paste0(common_name, "_", id)) %>%
  select(common_name, indiv, id, pheno_class_id, lat, lon, elev, func_type,
         n_obs_r, n_yrs_r)

eo_table <- quants_r_eo %>%
  select(indiv, state, pheno_class_id, first_yes, min, early, outlier) %>%
  left_join(current_join, by = c("indiv", "pheno_class_id")) %>%
  mutate(early = ifelse(early == 1, "Yes", "No"),
         outlier = ifelse(outlier == 1, "Yes", "No"), 
         earliest = ifelse(first_yes <= min, "Yes", "No"),
         firstyes = parse_date_time(x = paste(year, first_yes), 
                                    orders = "yj")) %>%
  mutate(func_type2 = case_when(
    func_type == "Deciduous broadleaf" ~ "DB",
    func_type == "Deciduous conifer" ~ "DC",
    func_type == "Drought deciduous broadleaf" ~ "DDB",
    func_type == "Evergreen broadleaf" ~ "EB",
    func_type == "Evergreen conifer" ~ "EC",
    func_type == "Graminoid" ~ "Gram",
    func_type == "Semi-evergreen broadleaf" ~ "SEB",
    func_type == "Forb" ~ "Forb",
    .default = NA
  )) %>%
  select(pheno_class_id, common_name, func_type, id, state, lat, lon, elev, 
         first_yes, firstyes, early, earliest, outlier, n_yrs_r, n_obs_r) %>%
  arrange(pheno_class_id, func_type, common_name, state, first_yes) %>%
  select(-first_yes)


# Not sure if we'll use stuff below -------------------------------------------#


# Compare current year observations with distribution of first observation dates
# for plants in the same state ------------------------------------------------#

current <- current %>%
  mutate(state_spp = paste0(state, "_", common_name),
         state_spp_ph = paste0(state_spp, "_", pheno_class_id))
prior <- prior %>%
  mutate(state_spp = paste0(state, "_", common_name),
         state_spp_ph = paste0(state_spp, "_", pheno_class_id))

# Identify early or outlier observations based on previous observations in 
# same state
prior_st <- prior %>%
  filter(state_spp_ph %in% unique(current$state_spp_ph))
state_spp_n <- prior_st %>%
  group_by(state_spp, pheno_class_id, state_spp_ph) %>%
  summarize(n_obs = n(),
            n_indiv = n_distinct(id),
            n_years = n_distinct(year),
            .groups = "keep") %>%
  data.frame()

# Only evaluate current-year anomalies for state-species combinations that have
# sufficient observations in prior years
state_spp_n <- state_spp_n %>%
  filter(n_years >= min_yrs & n_obs >= min_obs)

prior_st <- prior_st %>%
  filter(state_spp_ph %in% state_spp_n$state_spp_ph)
quantiles_st <- prior_st %>%
  group_by(state_spp_ph) %>%
  summarize(min = min(first_yes),
            qearly = quantile(first_yes, lowerq),
            q0.25 = quantile(first_yes, 0.25),
            q0.75 = quantile(first_yes, 0.75),
            IQR = IQR(first_yes)) %>%
  data.frame() %>%
  mutate(whisker = q0.25 - 1.5 * IQR,
         outlier_threshold = ifelse(whisker < min, min, whisker))
current_st <- current %>%
  left_join(quantiles_st, by = "state_spp_ph") %>%
  mutate(early = ifelse(first_yes < qearly, 1, 0),
         outlier = ifelse(first_yes < outlier_threshold, 1, 0)) %>%
  # Remove state_spp_ph combinations that don't have enough data in prior years
  filter(!is.na(early))

# Look at what's left
current_st %>%
  select(pheno_class_id, state_spp, func_type, first_yes, min, 
         outlier_threshold, q0.25, q0.75, early, outlier) %>%
  filter(early == 1 | outlier == 1) %>%
  arrange(pheno_class_id, func_type, state_spp, first_yes)

prior %>%
  filter(state_spp == "CA_coyotebrush" & pheno_class_id == 7) %>%
  ggplot(aes(first_yes)) +
  geom_histogram(bins = 60)
# Coyotebrush histogram highlights why boxplots might be deceiving. It's a 
# multimodal distribution, with the highest frequency in fall, but a distinct
# peak at the beginning of the year. 

# Extract just early/outlier observations
current_steo <- current_st %>%
  filter(early == 1 | outlier == 1) %>%
  select(common_name, state, state_spp, pheno_class_id, state_spp_ph, func_type, 
         first_yes, min, qearly, q0.25, q0.75, early, outlier) %>%
  mutate(panel = paste0(common_name, " (", state, ")"),
         eo = ifelse(outlier == 1, "Outlier", "Early"),
         eo = factor(eo, levels = c("Early", "Outlier")))

# Loop through phenophase classes
for (phc in pheno_class_ids) {
  current_steo_p <- filter(current_steo,  pheno_class_id == phc)
  current_eo_stspp <- sort(unique(current_steo_p$state_spp))
  current_eo_spp <- str_split_i(current_eo_stspp, "_", 2)
  current_eo_st <- str_split_i(current_eo_stspp, "_", 1)
  
  # Put prior data into better format for plotting
  for (i in 1:length(current_eo_stspp)) {
    prior_spp <- prior %>%
      filter(common_name == current_eo_spp[i] & pheno_class_id == phc) %>%
      select(common_name, state_spp, state, func_type, year, first_yes) %>%
      mutate(group = "All") %>%
      mutate(panel = paste0(current_eo_spp[i], " (", current_eo_st[i], ")"))
    prior_stspp <- prior %>%
      filter(state_spp == current_eo_stspp[i] & pheno_class_id == phc) %>%
      select(common_name, state_spp, state, func_type, year, first_yes) %>%
      mutate(group = "State") %>%
      mutate(panel = paste0(current_eo_spp[i], " (", current_eo_st[i], ")"))
    prior_plot_new <- rbind(prior_spp, prior_stspp)
    if (i == 1) {
      prior_plot <- prior_plot_new
    } else {
      prior_plot <- rbind(prior_plot, prior_plot_new)
    }
  }
  
  phc_name <- paste0(phenoclasses$name[phc], ", first observation of the year")
  
  plot_temp <- prior_plot %>%
    ggplot(aes(x = first_yes, fill = group)) +
    geom_histogram(bins = 60, position = "identity") +
    scale_fill_manual(values = c("gray", "steelblue3")) +
    facet_wrap(~panel, scales = "free_y", ncol = 2) +
    geom_vline(data = current_steo_p,
               aes(xintercept = first_yes, linetype = eo), 
               color = "red") +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Day of year", y = "Count", title = phc_name,
         fill = "Prior observations", linetype = paste0(year, " observations")) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.spacing.x = unit(1, "cm"))
  assign(paste0("plot_", phc), plot_temp)
}

for (phc in pheno_class_ids) {
  print(get(paste0("plot_", phc)))
}

# Evaluate how much data we have at each site ---------------------------------#
  site_spp_ph <- prior %>%
    group_by(site_id, lat, lon, elev, common_name, pheno_class_id) %>%
    summarize(n_obs = n(),
              n_indiv = n_distinct(id),
              n_years = n_distinct(year),
              .groups = "keep") %>%
    mutate(spp_site_ph = paste0(common_name, "_", site_id, "_", pheno_class_id)) %>%
    data.frame()
  
  # Calculate proportion of site-species combos that have 5/10 years of data for 
  # each phenophase class
  site_spp_ph %>%
    group_by(pheno_class_id) %>%
    summarize(n_sites = n_distinct(site_id),
              # Number of sites with 10+ years of data
              n_sites_10 = n_distinct(site_id[n_years >= 10]),
              # Number of sites with 5+ years and 10+ observations
              n_sites_5_10 = n_distinct(site_id[n_years >= 5 & n_obs >= 10])) %>%
    mutate(prop_sites_10 = n_sites_10 / n_sites,
           prop_sites_5_10 = n_sites_5_10 / n_sites) %>%
    data.frame()
