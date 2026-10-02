load(".r_session_state.RData")
for(i in 1:length(packages_required)) {
  suppressPackageStartupMessages(library(packages_required[i], character.only = T))
}

tmap_mode("plot")

lisa_map_cluster1 <- borders1 +
  hh_map1 + ll_map1 + hl_map1  +lh_map1  + ns_map1 +   # combining all quadrant maps
  tm_compass(position = c(0.85, 0.85)) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_add_legend(type = "fill", col = c("royalblue2", "red2", "darkgreen", "gold", "lightgrey"), 
                labels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant"), title = "LISA cluster\n(Median Household Income)") +
  tm_layout(
    legend.title.size = 0.715,
    legend.text.size = 0.6, 
    inner.margins = c(0.22, 0.13, 0.02, 0.08), 
    legend.position = c("left", "bottom"), 
    bg.color = "#eaf5fa",
    main.title = "LISA Clusters for Median Household Income Variable", 
    main.title.size = 0.7,
    main.title.position = "center"
  )


lisa_map_cluster2 <- borders2 +
  hh_map2 + ll_map2 + hl_map2   + ns_map2 +  #  Combining All quadrant maps
  tm_compass(position = c(0.85, 0.85)) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_add_legend(type = "fill", col = c("royalblue2", "red2", "darkgreen", "lightgrey"), 
                labels = c("High-High", "Low-Low", "High-Low", "Not significant"), title = "LISA cluster\n(Broadband Internet Use)") +
  tm_layout(
    legend.title.size = 0.715,
    legend.text.size = 0.6, 
    inner.margins = c(0.185, 0.13, 0.02, 0.08), 
    legend.position = c("left", "bottom"), 
    bg.color = "#eaf5fa",
    main.title = "LISA Clusters for Broadband Internet Use Variable",  
    main.title.size = 0.7,
    main.title.position = "center"
  )

map3_combined <- tmap_arrange(lisa_map_cluster1, lisa_map_cluster2)

tmap_save(map3_combined, "tempRchunks/map3_combined.png", width = 12, height = 6, units = "in", dpi = 300)
