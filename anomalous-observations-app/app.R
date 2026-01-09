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
library(sf)
library(glue)
library(ggplot2)
library(shinycssloaders)

# Establish link to board
board <- board_connect(
  auth = "manual",
  server = Sys.getenv("CONNECT_SERVER"),
  key = Sys.getenv("CONNECT_API_KEY")
)

# Load data from pins
# current <- pin_read(board = board, "ezylstra/anom-data-current-yr")
# prior <- pin_read(board = board, "ezylstra/anom-data-prior-yrs")
# dist_matrix <- pin_read(board = board, "ezylstra/anom-data-dist-matrix")

# Load data from pins (2025 dashboard)
current <- pin_read(board = board, "ezylstra/anom-data-current-yr-25")
prior <- pin_read(board = board, "ezylstra/anom-data-prior-yrs-25")
dist_matrix <- pin_read(board = board, "ezylstra/anom-data-dist-matrix-25")

# ui --------------------------------------------------------------------------#

ui <- page_navbar(
  title = "Anomalous observations",
  sidebar = sidebar(
    p(strong(em("Testing with 2025 data"))),
    width = "20%",
    div(style = "font-size:80%",
    navset_tab(
      nav_panel(
        title = "Settings",
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
      nav_panel(
        title = "Info",
        br(),
        p("This is where we'll include a brief description of the app and
          descriptions of the various options in the settings.")
      )
    )
    )
  ),
  fillable = "Map",
  nav_panel(title = "Data", 
            downloadButton("download", "Download filtered data"),
            div(withSpinner(DTOutput(outputId = "table"), type = 6, 
                            color = "#0dc5c1"), style = "font-size:80%")),
  nav_panel(title = "Map", 
            leafletOutput(outputId = "map"),
            em("Locations in map have been adjusted slightly to limit overlap")),
)

server <- function(input, output, session) {

  # Convert to data.table for performance
  current_dt <- as.data.table(current)
  prior_dt <- as.data.table(prior)
  
  # Add unique row IDs to current observations
  current_dt[, row_id := .I]
  
  # Extract 
  current_sites <- as.numeric(rownames(dist_matrix))
  prior_sites <- as.numeric(colnames(dist_matrix))

  # Get input values
  qearly <- reactive(input$inEarly)
  radius <- reactive(input$inRadius)
  elev_buffer <- reactive(input$inElev)
  min_obs <- reactive(input$inObs)
  min_yrs <- reactive(input$inYrs)
  phps <- reactive(input$inPhps)
  
  # Create lookup table of all site pairs within selected radius (in km)
  within_radius <- reactive({which(dist_matrix <= radius(), arr.ind = TRUE)})
  site_pairs <- reactive({
    data.table(
      current_site = current_sites[within_radius()[, 1]],
      prior_site = prior_sites[within_radius()[, 2]]
    )
  })
  
  # Expand current_dt to include all nearby prior sites
  current_expanded <- reactive({
    merge(
      current_dt,
      site_pairs(),
      by.x = "site_id",
      by.y = "current_site",
      allow.cartesian = TRUE
    )
  })
  # This will have as many rows for each observation in current as there are 
  # sites that are within the given radius (does NOT depend on spp, class, elev)
    
  # Join with prior observations (matching spp, phenophase class, nearby sites)
  matches_all <- reactive({
    merge(
      current_expanded(),
      prior_dt,
      by.x = c("common_name", "pheno_class_id", "prior_site"),
      by.y = c("common_name", "pheno_class_id", "site_id"),
      suffixes = c("", "_prior"),
      allow.cartesian = TRUE
    )
  })
  # For any row that has prior observations for the site within xx km and the
  # current species and phenophase class, attach prior observation info (incl 
  # elev, year, first_yes). IF there are no prior observations for any site 
  # within xx km, then rows for that current observation are removed. 
  
  # Filter by elevation range
  matches <- reactive({
    matches_all()[
      elev_prior > (elev - elev_buffer()) & elev_prior < (elev + elev_buffer())
    ]
  })
  
  # Calculate statistics for each current observation
  stats <- reactive({ 
    matches()[, .(
      prior_nobs = .N,
      prior_nyrs = uniqueN(year_prior),
      prior_min = min(first_yes_prior),
      prior_qearly = quantile(first_yes_prior, qearly())
    ), by = row_id]
  })  
    
  # Join statistics back to current_df
  result_temp <- reactive({
    setkey(current_dt, row_id)
    setkey(stats(), row_id)
    stats()[current_dt, on = "row_id"]
  })  
    
  # Fill in NAs for observations with no matches
  result <- reactive({
    result_temp()[is.na(prior_nobs), `:=`(
      prior_nobs = 0L,
      prior_nyrs = 0L,
      prior_min = NA_real_,
      prior_qearly = NA_real_
    )]
  })  
  # In data.table syntax: X[Y] means "keep all rows from Y, add matching data from X"
  # So: stats[current_dt] = Keep ALL current_dt rows, add stats where available
    
  # Get anom obs dataframe
  current_df <- reactive({
    result() %>%
      as.data.frame() %>%
      filter(pheno_class_id %in% phps()) %>%
      mutate(suff_obs = ifelse(prior_nobs >= min_obs() & prior_nyrs >= min_yrs(), 
                               1, 0)) %>%
      mutate(early = ifelse(suff_obs == 1 & first_yes <= prior_qearly, 1, 0),
             earliest = ifelse(suff_obs == 1 & first_yes <= prior_min, 1, 0)) %>%
      separate_wider_delim(phenophase, 
                           delim = " (", 
                           names = c("php_simple", NA),
                           too_few = "align_start",
                           cols_remove = TRUE) %>%
      mutate(first_date = parse_date_time(x = paste(year, first_yes), 
                                          orders = "yj")) %>%
      mutate(first_date = format(first_date, "%Y-%m-%d")) %>%
      mutate(early_cat = ifelse(earliest == 1, "Earliest",
                                ifelse(early == 1, "Early", "Not early"))) %>%
      filter(early == 1 | earliest == 1) %>%
      arrange(pheno_class_id, func_type, common_name, state, site_id, id) %>%
      select(php_simple, common_name, func_type, state, 
             first_date, early_cat, prior_nyrs, prior_nobs, prior_no, 
             id, site_id, lat, lon, pheno_class_id) %>%
      mutate(php_simple = factor(php_simple),
             common_name = factor(common_name),
             func_type = factor(func_type),
             state = factor(state),
             id = factor(id),
             first_date = as.Date(first_date),
             early_cat = factor(early_cat),
             early_percentile = qearly(),
             elevation_buffer_m = elev_buffer(),
             geographic_radius_km = radius(),
             id_class = paste0(id, "_", pheno_class_id))
  })

  output$table <- renderDT({
    newtable <- current_df() %>%
      mutate(
        actionable = glue('<button id="custom_btn" onclick="Shiny.onInputChange(\'button_id\', \'{id_class}\')">Histogram</button>')
      ) %>%
      relocate(actionable)
    
    datatable(newtable,
              escape = FALSE,
              selection = "single",
              colnames = c("Histogram",
                           "Phenophase", 
                           "Species",
                           "Functional type",
                           "State",
                           "First yes",
                           "Early/Earliest",
                           "Comparison yrs",
                           "Comparison obs",
                           "Prior no",
                           "Plant ID",
                           "Site ID",
                           "Latitude",
                           "Longitude",
                           "Phenophase class ID",
                           "Early percentile",
                           "Elevation buffer (m)",
                           "Geographic radius (km)",
                           "Id_PhenoClass"),
              extensions = 'Buttons',
              options = list(
                server = FALSE,
                scrollX = TRUE,
                scrollY = TRUE,
                dom = 'lrtip',
                pageLength = 10,
                lengthMenu = c(10, 50, 100),
                autoWidth = TRUE,
                # Hiding some columns, but keeping in table for downloads
                # Centering some columns
                columnDefs = list(list(targets = 15:19, visible = FALSE),
                                  list(targets = 4:11, className = 'dt-center'))),
              filter = "top") %>%
      formatStyle("php_simple", "white-space"="nowrap") %>%
      formatStyle("common_name", "white-space"="nowrap") %>%
      formatStyle("func_type", "white-space"="nowrap") %>%
      formatStyle("first_date", "white-space"="nowrap")
  })

  observeEvent(input$button_id, {
    
    select_id <- str_split_1(input$button_id, pattern = "_")[1]
    select_php <- str_split_1(input$button_id, pattern = "_")[2]
    
    # Clear the old histogram output
    output$histogram <- renderPlot({
      NULL
    })
    
    # Show a modal with spinner immediately
    showModal(modalDialog(
      title = "Loading histogram...",
      withSpinner(plotOutput("histogram"), type = 6, color = "#0dc5c1"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
    
    hist_data <- matches() %>%
      filter(id == select_id & pheno_class_id == select_php)
    
    output$histogram <- renderPlot({
      ggplot(hist_data, aes(x = first_yes_prior)) +
        geom_histogram(binwidth = 1, color = "gray",
                       aes(fill = "Comparison observations (prior years)")) +
        geom_vline(aes(xintercept = hist_data$first_yes[1],
                       color = "Focal plant (current year)"), 
                   linetype = "dashed") +
        scale_fill_manual("", values = "gray30") +
        scale_color_manual("", values = "dodgerblue4") +
        labs(x = "First yes day of year", y = "No. comparison observations",
             title = paste0(hist_data$common_name, ", ", hist_data$phenophase,
                            "\nPlant ID: ", hist_data$id, "; Site ID: ", 
                            hist_data$site_id)) +
        theme(plot.title = element_text(size = 15),
              axis.text = element_text(size = 15),
              axis.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              legend.position = "bottom")
    })
    
    # Update modal title after plot is created
    showModal(modalDialog(
      title = NULL,
      plotOutput("histogram"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  filtered_frame <- reactive({
    frame <- req(current_df())
    indexes <- req(input$table_rows_all)
    frame[indexes,]
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste("anom_obs_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(filtered_frame(), file)
    }
  )
  
  # Jittering locations (slightly) in map data to see things better
  map_frame <- reactive({
    frame <- filtered_frame() %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
      st_jitter(factor = 0.0001) 
    frame %>%
      cbind(st_coordinates(frame)) %>%
      rename(lon = X, lat = Y) %>%
      st_drop_geometry()
  })
  
  # Create base map
  output$map <- renderLeaflet({
    leaflet(data = map_frame()) %>%
      fitBounds(lng1 = -125, lat1 = 26, lng2 = -65, lat2 = 47) %>%
      addTiles(options = tileOptions(opacity = 0.7)) %>%
      addCircleMarkers(
        lng = ~lon, 
        lat = ~lat,
        data = filter(map_frame(), pheno_class_id == 1),
        group = "Breaking leaf buds",
        radius = 5,
        fillColor = "green", 
        fillOpacity = 0.7,
        stroke = FALSE,
        popup = ~paste0(common_name, "<br>",
                        "ID: ", id, "<br>",
                        "First yes: ", first_date, "<br>",
                        "No. comparison years: ", prior_nyrs, "<br>",
                        "No. comparison obs: ", prior_nobs)
        ) %>%
      addCircleMarkers(
        lng = ~lon, 
        lat = ~lat,
        data = filter(map_frame(), pheno_class_id == 3),
        group = "Leaves",
        radius = 5,
        fillColor = "blue", 
        fillOpacity = 0.7,
        stroke = FALSE,
        popup = ~paste0(common_name, "<br>",
                        "ID: ", id, "<br>",
                        "First yes: ", first_date, "<br>",
                        "No. comparison years: ", prior_nyrs, "<br>",
                        "No. comparison obs: ", prior_nobs)
      ) %>%
      addCircleMarkers(
        lng = ~lon, 
        lat = ~lat,
        data = filter(map_frame(), pheno_class_id == 6),
        group = "Flowers",
        radius = 5,
        fillColor = "orange", 
        fillOpacity = 0.7,
        stroke = FALSE,
        popup = ~paste0(common_name, "<br>",
                        "ID: ", id, "<br>",
                        "First yes: ", first_date, "<br>",
                        "No. comparison years: ", prior_nyrs, "<br>",
                        "No. comparison obs: ", prior_nobs)
      ) %>%
      addCircleMarkers(
        lng = ~lon, 
        lat = ~lat,
        data = filter(map_frame(), pheno_class_id == 7),
        group = "Open flowers",
        radius = 5,
        fillColor = "red", 
        fillOpacity = 0.7,
        stroke = FALSE,
        popup = ~paste0(common_name, "<br>",
                        "ID: ", id, "<br>",
                        "First yes: ", first_date, "<br>",
                        "No. comparison years: ", prior_nyrs, "<br>",
                        "No. comparison obs: ", prior_nobs)
      ) %>%
      addLayersControl(overlayGroups = c("Breaking leaf buds",
                                         "Leaves",
                                         "Flowers",
                                         "Open flowers"),
                       options = layersControlOptions(collapse = FALSE)) %>%
      addLegend(position = "bottomright",
                colors = c("green", "blue", "orange", "red"),
                labels = c("Breaking leaf buds", "Leaves",
                           "Flowers", "Open flowers"),
                opacity = 1)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
