# Anomalous observations script - dashboard
# 20 November 2025

library(dplyr)
library(lubridate)
library(stringr)
library(data.table)
library(ggplot2)
library(scales) # To make integer axis labels
# library(sf)
library(terra)
# library(tidyterra)
library(leaflet)
library(pins)

# This code identifies anomalous observations of plants in the current year 
# compared to the long term record (2009-XXXX). First observation date in 
# current year is compared to the distribution of dates from prior years. Early 
# (lower than Xth percentile) and outliers (earlier than 1.5x the IQR) are 
# flagged.


# Establish link to board
board <- board_connect(
  auth = "manual",
  server = Sys.getenv("CONNECT_SERVER"),
  key = Sys.getenv("CONNECT_API_KEY")
)

# Load data from pins
# usstates <- pin_read(board = board, "ezylstra/us-states")
current <- pin_read(board = board, "ezylstra/anom-data-current-yr")
prior <- pin_read(board = board, "ezylstra/anom-data-prior-yrs")
dist_matrix <- pin_read(board = board, "ezylstra/anom-data-dist-matrix")

# Set parameters --------------------------------------------------------------#

# Define early observations (below X percentile)
lowerq <- 0.05

# Logical indicating whether to denote Tukey outliers in addition to early 
# observations (those below X percentile) in histograms
# outliers <- FALSE

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
# reported (must be < 30 since we already filtered out observations with 
# longer gaps)
prior_no_max <- 30

# Select phenophase classes we're interested in
pheno_class_ids <- c(1, 3, 6, 7)
  # 1: initial shoot/leaf growth
  # 3: leaves/needles
  # 6: flowers/cones
  # 7: open flowers/cones 

# Set mapping parameters ------------------------------------------------------#

# # Set date to extract AGDD anomalies (for map). Latest possible date = 30 April. 
# # If don't want today, then set agdd_date to desired date.
# agdd_today <- TRUE
# if (!agdd_today) {
#   agdd_date <- ymd("2025-02-03")
# } else {
#   agdd_date <- min(today(), ymd(paste0(year(Sys.Date()), "-04-30")))
# }
# 
# # Define breaks that will be used to delineate areas cooler/warmer than 30-year
# # normals. Breaks = AGDD anomolies on agdd_date.
# agdd_breaks <- c(-2000, -40, 40, 2000)
# 
# # Colors for map (blue = cooler, offwhite = normal, orange = warmer)
# cw_cols <- c("#6fa8d6", "#f7f0da", "#fa824d")
# # plot(1:3, rep(1, 3), pch = 19, cex = 5, col = cols)

# Compare current year observations with distribution of first observation dates
# for plants within xx km and xx m elevation ----------------------------------#

# Convert to data.table for performance
current_dt <- as.data.table(current)
prior_dt <- as.data.table(prior)

# Add unique row IDs to current observations
current_dt[, row_id := .I]

# Extract 
current_sites <- as.numeric(rownames(dist_matrix))
prior_sites <- as.numeric(colnames(dist_matrix))
# Checks:
# all.equal(unique(current_dt$site_id), current_sites)
# all.equal(unique(prior_dt$site_id), prior_sites)

# Create lookup table of all site pairs within selected radius (in km)
within_radius <- which(dist_matrix <= radius, arr.ind = TRUE)
site_pairs <- data.table(
  current_site = current_sites[within_radius[, 1]],
  prior_site = prior_sites[within_radius[, 2]]
)

# Expand current_dt to include all nearby prior sites
current_expanded <- merge(
  current_dt,
  site_pairs,
  by.x = "site_id",
  by.y = "current_site",
  allow.cartesian = TRUE
)
# This will have as many rows for each observation in current as there are 
# sites that are within the given radius (does NOT depend on spp, class, elev)

# Join with prior observations (matching spp, phenophase class, nearby sites)
matches <- merge(
  current_expanded,
  prior_dt,
  by.x = c("common_name", "pheno_class_id", "prior_site"),
  by.y = c("common_name", "pheno_class_id", "site_id"),
  suffixes = c("", "_prior"),
  allow.cartesian = TRUE
)
# For any row that has prior observations for the site within xx km and the
# current species and phenophase class, attach prior observation info (incl 
# elev, year, first_yes). IF there are no prior observations for any site 
# within xx km, then rows for that current observation are removed. 

# Filter by elevation range
matches <- matches[
  elev_prior > (elev - elev_buffer) & elev_prior < (elev + elev_buffer)
]

# Calculate statistics for each current observation
stats <- matches[, .(
  prior_nobs = .N,
  prior_nyrs = uniqueN(year_prior),
  prior_min = min(first_yes_prior),
  prior_qearly = quantile(first_yes_prior, 0.05)
), by = row_id]

# Join statistics back to current_df
setkey(current_dt, row_id)
setkey(stats, row_id)
result <- stats[current_dt, on = "row_id"]
# In data.table syntax: X[Y] means "keep all rows from Y, add matching data from X"
# So: stats[current_dt] = Keep ALL current_dt rows, add stats where available

# Fill in NAs for observations with no matches
result[is.na(prior_nobs), `:=`(
  prior_nobs = 0L,
  prior_nyrs = 0L,
  prior_min = NA_real_,
  prior_qearly = NA_real_
)]

# Remove row_id column
# result[, row_id := NULL]

# Convert back to data frame
current_df <- as.data.frame(result)

# Final calculations
current_df <- current_df %>%
  mutate(suff_obs = ifelse(prior_nobs >= min_obs & prior_nyrs >= min_yrs, 
                           1, 0)) %>%
  mutate(early = ifelse(suff_obs == 1 & first_yes <= prior_qearly, 1, 0),
         earliest = ifelse(suff_obs == 1 & first_yes <= prior_min, 1, 0))

# Checks....
# count(current_df, suff_obs, early, earliest)
# filter(current_df, early == 1 & earliest == 0) %>% head(10)

# Create histogram for a particular observation:
matches %>% 
  filter(row_id == 212) %>%
  ggplot(aes(x = first_yes_prior)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = filter(matches, row_id == 212) %>% 
               pull(first_yes) %>% 
               first(),
             color = "blue", linetype = "dashed") +
  labs(x = "First yes day of year", y = "Count")


# -----------------------------------------------------------------------------#
# Comparison of different approaches (loops w dplyr vs data.table) ------------#
current <- pin_read(board = board, "ezylstra/anom-data-current-yr")
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
           early = ifelse(first_yes <= qearly, 1, 0),
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
         firstyes = parse_date_time(x = paste(2025, first_yes), 
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

# Compare results with what's above:
count(eo_table, early, earliest)
count(current_df, suff_obs, early, earliest)
# Counts of early, earliest are the same
