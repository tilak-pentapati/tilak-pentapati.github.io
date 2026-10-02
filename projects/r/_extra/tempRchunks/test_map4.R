load(".r_session_state.RData")
for(i in 1:length(packages_required)) {
  suppressPackageStartupMessages(library(packages_required[i], character.only = T))
}

tmap_mode("plot")

map_cluster1 <- tm_basemap() +
  tm_shape(overlay_listings_zip) +
  tm_fill(col = "k6cls", palette = viridis(6), style = "cat", title = "Type of\nCluster") +
  tm_compass(position = c(0.85, 0.88)) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_borders(col = "white", lwd = 0.2)+
  tm_layout(
    legend.title.size = 0.7,
    legend.text.size = 0.6, 
    inner.margins = c(0.2, 0.13, 0.02, 0.08), 
    legend.position = c("left", "bottom"), 
    bg.color = "#eaf5fa",
    main.title = "Type of Clusters in Bay Area",  # Add the map title within tm_layout
    main.title.size = 0.7,
    main.title.position = "center"
  )

map_veh <- tm_basemap() +
  tm_shape(overlay_listings_zip) +
  tm_fill(col = "No_of_Light_Duty_Vehicles", palette = "YlGn", n = 5, style = "fisher", title = "Vehicles\nper ZIP") +
  tm_compass(position = c(0.85, 0.88)) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_borders(col = "black", lwd = 0.2)+
  tm_layout(
    legend.title.size = 0.7,
    legend.text.size = 0.6, 
    inner.margins = c(0.185, 0.13, 0.02, 0.08), 
    legend.position = c("left", "bottom"), 
    bg.color = "#eaf5fa",
    main.title = "No. of Vehicles(Light-Duty) in Bay Area",  # Add the map title within tm_layout
    main.title.size = 0.7,
    main.title.position = "center"
  )


map_electric <- tm_basemap() +
  tm_shape(overlay_listings_zip) +
  tm_fill(col = "Combined_Avg_Comm_Rate", palette = "YlGn", n = 5, style = "fisher", title = "Commerical\nElectrical Rate") +
  tm_compass(position = c(0.85, 0.88)) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_borders(col = "black", lwd = 0.2)+
  tm_layout(
    legend.title.size = 0.7,
    legend.text.size = 0.6, 
    inner.margins = c(0.185, 0.13, 0.02, 0.08), 
    legend.position = c("left", "bottom"), 
    bg.color = "#eaf5fa",
    main.title = "Commercial Electrical Rate in Bay Area",  # Add the map title within tm_layout
    main.title.size = 0.7,
    main.title.position = "center"
  )


map4_combined <- tmap_arrange(map_cluster1, map_veh, map_electric, nrow = 1, ncol = 3)
tmap_save(map4_combined, "tempRchunks/map4_combined.png", width = 18, height = 6, units = "in", dpi = 300)
