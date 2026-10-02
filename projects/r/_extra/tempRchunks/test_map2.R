load(".r_session_state.RData")
for(i in 1:length(packages_required)) {
  suppressPackageStartupMessages(library(packages_required[i], character.only = T))
}

tmap_mode("plot")
map2.1 <- tm_basemap() +
  tm_shape(bayarea_zip_vars1) +    # shapefile for map 2
  tm_fill(col = "estimate_MedianHhI", palette = "YlGn", n = 5, style = "fisher", title = "Median Household Income\n(in $ USD)") + # estimate_MedianHhI
  tm_compass(position = c(0.85, 0.85)) +  # Adding Compass
  tm_scalebar(position = c("right", "bottom")) + # Adding Scalebar
  tm_borders(lwd = 0.3)+      # Borders for zcta's polygons
  tm_layout(
    legend.title.size = 0.715,
    legend.text.size = 0.6, 
    inner.margins = c(0.22, 0.13, 0.02, 0.08), 
    legend.position = c("left", "bottom"), 
    bg.color = "#eaf5fa",
    main.title = "Median Household Income per ZIP Code",  # Add the map title within tm_layout
    main.title.size = 0.7,
    main.title.position = "center"
  )

map2.2 <- tm_basemap() +
  tm_shape(bayarea_zip_vars) +
  tm_fill(col = "estimate_InternetUse", palette = "YlGn", n = 5, style = "fisher", title = "Broadband Internet Connection\n(People)") +
  tm_compass(position = c(0.85, 0.85)) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_borders(lwd = 0.3)+
  tm_layout(
    legend.title.size = 0.715,
    legend.text.size = 0.6, 
    inner.margins = c(0.185, 0.13, 0.02, 0.08), 
    legend.position = c("left", "bottom"), 
    bg.color = "#eaf5fa",
    main.title = "No. of People aged (18-64) using Broadband Internet per ZIP Code",  # Add the map title within tm_layout
    main.title.size = 0.7,
    main.title.position = "center"
  )

tmap_save(tmap_arrange(map2.1, map2.2), "tempRchunks/map2_combined.png", width=12, height=6, dpi=300)
