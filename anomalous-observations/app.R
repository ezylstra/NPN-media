library(data.table)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(shiny)
library(bslib)
library(leaflet)
library(DT)
library(pins)

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


# ui --------------------------------------------------------------------------#

ui <- page_navbar(
  title = "Anomalous observations",
  sidebar = sidebar(
    width = "20%",
    checkboxGroupInput(inputId = "inPhps",
                       label = "Phenophases",
                       choices = c("Breaking leaf buds" = 1,
                                   "Leaves" = 3,
                                   "Flowers" = 6,
                                   "Open flowers" = 7),
                       selected = c(1, 3, 6, 7)),
    sliderInput(
      inputId = "inRadius",
      label = "Geographic radius (km)",
      min = 20,
      max = 1000,
      step = 20,
      value = 100
    ),
    sliderInput(
      inputId = "inElev",
      label = "Elevational buffer (m)",
      min = 100,
      max = 3000,
      step = 100,
      value = 1000
    ),
    sliderInput(
      inputId = "inYrs",
      label = "Min. prior years",
      min = 5,
      max = 10,
      step = 1,
      value = 8
    ),
    sliderInput(
      inputId = "inObs",
      label = "Min. prior observations",
      min = 10,
      max = 50,
      step = 1,
      value = 20
    ),
    sliderInput(
      inputId = "inEarly",
      label = "Early percentile",
      min = 0.025,
      max = 0.25,
      step = 0.025,
      value = 0.05
    )
  ),
  nav_panel(title = "Data", DTOutput(outputId = "table")),
  nav_panel(title = "Map", leafletOutput(outputId = "map"))
)

server <- function(input, output) {

  # Convert to data.table for performance
  current_dt <- as.data.table(current)
  prior_dt <- as.data.table(prior)
  
  # Add unique row IDs to current observations
  current_dt[, row_id := .I]
  
  # Extract 
  current_sites <- as.numeric(rownames(dist_matrix))
  prior_sites <- as.numeric(colnames(dist_matrix))

  output$table <- renderDT({
    
    qearly <- input$inEarly
    radius <- input$inRadius
    elev_buffer <- input$inRadius
    min_obs <- input$inObs
    min_yrs <- input$inYrs
    # blb <- input$inBLB
    # leaves <- input$inLeaves
    # flowers <- input$inFlowers
    # open <- input$inOpenFlowers
    phps <- input$inPhps
    
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
      prior_qearly = quantile(first_yes_prior, qearly)
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
    
    # Convert back to data frame
    current_df <- as.data.frame(result)
    current_df <- current_df %>%
      filter(pheno_class_id %in% phps)
    # if (blb == FALSE) {
    #   current_df <- filter(current_df, pheno_class_id != 1)
    # }
    # if (leaves == FALSE) {
    #   current_df <- filter(current_df, pheno_class_id != 3)
    # }
    # if (flowers == FALSE) {
    #   current_df <- filter(current_df, pheno_class_id != 6)
    # }
    # if (open == FALSE) {
    #   current_df <- filter(current_df, pheno_class_id != 7)
    # }
    
    # Final calculations
    current_df <- current_df %>%
      mutate(suff_obs = ifelse(prior_nobs >= min_obs & prior_nyrs >= min_yrs, 
                               1, 0)) %>%
      mutate(early = ifelse(suff_obs == 1 & first_yes <= prior_qearly, 1, 0),
             earliest = ifelse(suff_obs == 1 & first_yes <= prior_min, 1, 0))
    
    current_df <- current_df %>%
      separate_wider_delim(phenophase, 
                           delim = " (", 
                           names = c("php_simple", NA),
                           too_few = "align_start",
                           cols_remove = TRUE) %>%
      mutate(first_date = parse_date_time(x = paste(year, first_yes), 
                                          orders = "yj")) %>%
      mutate(first_date = format(first_date, "%Y-%m-%d")) %>%
      mutate(early_cat = ifelse(earliest == 1, "Earliest",
                                ifelse(early == 1, "Early", "Not early")))
    
    # Create table with summary of results (that can be filtered by user)
    outtable <- current_df %>%
      filter(early == 1 | earliest == 1) %>%
      arrange(pheno_class_id, func_type, common_name, state, site_id, id) %>%
      select(php_simple, common_name, func_type, state, site_id, id, 
             first_date, prior_no, early_cat, prior_nyrs, prior_nobs) %>%
      mutate(php_simple = factor(php_simple),
             common_name = factor(common_name),
             func_type = factor(func_type),
             state = factor(state),
             early_cat = factor(early_cat))

    datatable(outtable,
              colnames = c("Phenophase", 
                           "Species",
                           "Functional type",
                           "State",
                           "Site ID",
                           "Plant ID",
                           "First yes",
                           "Days since no",
                           "Early",
                           "Num. comparison yrs",
                           "Num. comparison obs"),
              extensions = 'Buttons',
              options = list(
                server = FALSE,
                scrollX = TRUE,
                scrollY = TRUE,
                dom = 'lfrtipB',
                buttons = c("copy", "csv"),
                pageLength = 10,
                lengthMenu = c(10, 50, 100),
                autoWidth = TRUE,
                # Hiding some columns, but keeping in table for downloads
                # Centering some columns
                columnDefs = list(list(targets = c(5, 8, 10:11), visible = FALSE),
                                  list(targets = 4:11, className = 'dt-center'))),
              filter = "top")
  })
  
  # Create base map
  output$map <- renderLeaflet({
    leaflet() %>%
      fitBounds(lng1 = -125, lat1 = 26, lng2 = -65, lat2 = 47) %>%
      addTiles()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
