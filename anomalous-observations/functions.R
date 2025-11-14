# Function to clean up individual phenometrics dataset after downloading
  # Remove unnecessary columns
  # Remove observations of 'ohi'a lehua  
  # Remove observations with conflicts or those without prior no in last X days
  # Keep only first observation of individual for a phenophase in a year
ind_ph_cleanup <- function(df, prior_no_mindays = 30) {
  df <- df %>%
    select(site_id, latitude, longitude, elevation_in_meters, state,
           common_name, species_id, species_functional_type, individual_id,
           phenophase_description, pheno_class_id, first_yes_year, first_yes_doy,
           numdays_since_prior_no, observed_status_conflict_flag) %>%
    rename(lat = latitude,
           lon = longitude,
           elev = elevation_in_meters,
           func_type = species_functional_type,
           id = individual_id,
           phenophase = phenophase_description,
           year = first_yes_year,
           first_yes = first_yes_doy,
           prior_no = numdays_since_prior_no,
           flag = observed_status_conflict_flag) %>%
    filter(common_name != "'ohi'a lehua") %>%
    filter(is.na(flag)) %>%
    filter(!is.na(prior_no)) %>%
    filter(prior_no <= prior_no_mindays) %>%
    select(-flag) %>%
    arrange(id, pheno_class_id, year, first_yes) %>%
    distinct(id, pheno_class_id, year, .keep_all = TRUE) %>%
    data.frame()
  df
}

# Function to fill in state postal code if missing
state_fill_in <- function(df, states_shp = states) {
  state_fill <- filter(df, is.na(state)) %>%
    select(site_id, lon, lat) %>%
    distinct()
  state_fillv <- vect(state_fill, crs = "epsg:4326")
  state_new <- terra::extract(states_shp, state_fillv)
  state_fill <- cbind(state_fill, state_new = state_new$postal)
  df <- df %>%
    left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
    mutate(state = ifelse(!is.na(state), state, state_new)) %>%
    select(-state_new) %>%
    # Restrict observations to the continental US
    filter(state %in% states48)
  df
}

# Function to fill in elevation if missing (using elevatr package)
elev_fill_in <- function(df) {
  elev_fill <- filter(df, is.na(elev)) %>%
    distinct(site_id, lat, lon) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
  elev_fill <- get_elev_point(locations = elev_fill, src = "epqs")
  elev_fill <- data.frame(elev_fill) %>%
    mutate(elev_new = round(elevation)) %>%
    select(site_id, elev_new)
  df <- df %>%
    left_join(elev_fill, by = "site_id") %>%
    mutate(elev = ifelse(!is.na(elev), elev, elev_new)) %>%
    select(-elev_new)
  df
}