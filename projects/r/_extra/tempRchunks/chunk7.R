load(".r_session_state.RData")
for(i in 1:length(packages_required)) {
  suppressPackageStartupMessages(library(packages_required[i], character.only = T))
}

# Chunk 25
tmap_mode("plot")
map2.1 <- tm_basemap() +
  tm_shape(bayarea_zip_vars1) +    # shapefile for map 2
  tm_fill(col = "estimate_MedianHhI", palette = "YlGn", n = 5, style = "fisher", title = "Median Household Income (in $ USD)") + # estimate_MedianHhI
  tm_compass(position = c(0.85, 0.85)) +  # Adding Compass
  tm_scalebar(position = c(0.625, 0.02)) + # Adding Scalebar
  tm_borders(lwd = 0.3)+      # Borders for zcta's polygons
  tm_layout(
    legend.title.size = 0.715,
    legend.text.size = 0.6, 
    inner.margins = c(0.22, 0.13, 0.02, 0.08), 
    legend.position = c(0.05, 0.02), 
    legend.width = 0.5, 
    bg.color = "#eaf5fa",
    main.title = "Median Household Income per ZIP Code",  # Add the map title within tm_layout
    main.title.size = 0.7,
    main.title.position = "center"
  )


map2.2 <- tm_basemap() +
  tm_shape(bayarea_zip_vars) +
  tm_fill(col = "estimate_InternetUse", palette = "YlGn", n = 5, style = "fisher", title = "Broadband Internet Connection (People)") +
  tm_compass(position = c(0.85, 0.85)) +
  tm_scalebar(position = c(0.625, 0.02)) +
  tm_borders(lwd = 0.3)+
  tm_layout(
    legend.title.size = 0.715,
    legend.text.size = 0.6, 
    inner.margins = c(0.185, 0.13, 0.02, 0.08), 
    legend.position = c(0.05, 0.02), 
    legend.width = 0.5, 
    bg.color = "#eaf5fa",
    main.title = "No. of People aged (18-64) using Broadband Internet per ZIP Code",  # Add the map title within tm_layout
    main.title.size = 0.7,
    main.title.position = "center"
  )

# Chunk 26
print(tmap_arrange(map2.1, map2.2))

# Chunk 27
bayarea_zip_vars_map3 <- bayarea_zip_vars %>%
  select(estimate_InternetUse, everything())  # Changing the order of attributes for interactive map

# Chunk 28
listings_bayarea_sf_transform_copy <- listings_bayarea_sf_transform
listings_bayarea_sf_transform_copy$nl_price <- log(listings_bayarea_sf_transform_copy$price) # calculating natural logarithm of price ( log to base e)
print("CRS check for map 3:")
print(st_crs(bayarea_zip_vars_map3) == st_crs(listings_bayarea_sf_transform_copy)) # checking CRS of both shapefiles

# Chunk 29
listings_bayarea_sf_transform_copy <- listings_bayarea_sf_transform_copy %>% 
  select(nl_price, everything()) %>% mutate(nl_price = round(nl_price, 1))

# Chunk 30
tmap_mode("view") # interactive mode activating

map3 <- tm_shape(bayarea_zip_vars_map3) +
  tm_fill(col = "estimate_InternetUse", alpha = 1, palette = "YlGn", n = 5, style = "fisher",title = "Broadband Internet Connection(People)") +
  tm_borders() +
  tm_text("ZIP", size = 0.3)+
  tm_shape(listings_bayarea_sf_transform_copy) +
  tm_dots(col = "nl_price", palette = "Reds", n = 3, style = "fisher", title = "Natural Logarithm of Airbnb's Price") +
  tm_view(set.view = c(center_lon, center_lat, zoom_level))# Set legend position for interactive view

print(map3)

save.image(file=".r_session_state.RData")
