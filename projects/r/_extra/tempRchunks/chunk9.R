load(".r_session_state.RData")
for(i in 1:length(packages_required)) {
  suppressPackageStartupMessages(library(packages_required[i], character.only = T))
}

# Chunk 46
overlay_listings <- st_join(bayarea_zip_vars, listings_bayarea_sf_transform,  join = st_intersects) # spatial join to create Overlay Shapefiles

overlay_listings_zip <- overlay_listings  %>% 
  group_by(ZIP) %>% 
  summarise(
    count_airbnb =    sum(!is.na(price)),  # count non-NA price rows for count
    mean_price =      mean(price, na.rm = TRUE), # Mean price 
    mean_review  =    mean(number_of_reviews, na.rm = TRUE), # No. of Reviews
    mean_min_nights = mean(minimum_nights, na.rm = TRUE),   # No. of Minimum Nights stay
    Median_HhI =      mean(estimate_MedianHhI, na.rm = TRUE), # Median Household income
    Internet_Use =    mean(estimate_InternetUse, na.rm = TRUE) # Internet Usage
  )%>%
  filter(count_airbnb != 0) 

overlay_listings_zip <- overlay_listings_zip %>% left_join(vc_califorina, by = "ZIP") # Join overlay listings with Light-Duty Vehicles count
overlay_listings_zip <- overlay_listings_zip %>% left_join(electric_u, by = "ZIP") #Joining Overlay Listings with Commerical Electricity Rate.

overlay_listings_zip <- na.omit(overlay_listings_zip) # removing NA values

# Chunk 47
airbnb_ra <- c('count_airbnb', 'mean_price', 'mean_review', 'mean_min_nights', 'Median_HhI', 'Internet_Use', 'No_of_Light_Duty_Vehicles', 'Combined_Avg_Comm_Rate') # geodemographic arguments - (8 - tuple) , multivariate classification

# Chunk 48
overlay_listings_zip_NoGeo <- st_drop_geometry(overlay_listings_zip[, airbnb_ra,  drop=FALSE]) # we are dropping geometry col in shp

overlay_listings_zip_NoGeo <- overlay_listings_zip_NoGeo[complete.cases(overlay_listings_zip_NoGeo), ] # removing NA Values

# Chunk 49
set.seed(12345)  # setting seed

k6cls <- kmeans(overlay_listings_zip_NoGeo, centers=6, iter.max = 1000) # KNN clusteing - 6 clusters (Each Observation will be 8-tuple )
overlay_listings_zip$k6cls <- k6cls$cluster   # Type of cluster assigned bacnk to overlay listings

# Chunk 50
print("Centers for K-Means Clustering:")
print(k6cls$centers) # centres of 6 clusters, six 8-tuple vectors defines centre of each cluster

# Chunk 51
tmap_mode("plot")

map_cluster1 <- tm_basemap() +
  tm_shape(overlay_listings_zip) +
  tm_fill(col = "k6cls", palette = viridis(256), style = "cont", title = " Type of Cluster") +
  tm_compass(position = c(0.85, 0.88)) +
  tm_scalebar(position = c(0.58, 0.02)) +
  tm_borders(col = "white", lwd = 0.2)+
  tm_layout(
    legend.title.size = 0.7,
    legend.text.size = 0.6, 
    inner.margins = c(0.2, 0.13, 0.02, 0.08), 
    legend.position = c(0.05, 0.02), 
    legend.width = 0.5, 
    bg.color = "#eaf5fa",
    main.title = "Type of Clusters in Bay Area",  # Add the map title within tm_layout
    main.title.size = 0.7,
    main.title.position = "center"
  )

map_veh <- tm_basemap() +
  tm_shape(overlay_listings_zip) +
  tm_fill(col = "No_of_Light_Duty_Vehicles", palette = "YlGn", n = 5, style = "fisher", title = "Vehicles per ZIP") +
  tm_compass(position = c(0.85, 0.88)) +
  tm_scalebar(position = c(0.58, 0.02)) +
  tm_borders(col = "black", lwd = 0.2)+
  tm_layout(
    legend.title.size = 0.7,
    legend.text.size = 0.6, 
    inner.margins = c(0.185, 0.13, 0.02, 0.08), 
    legend.position = c(0.05, 0.02), 
    legend.width = 0.5, 
    bg.color = "#eaf5fa",
    main.title = "No. of Vehicles(Light-Duty) in Bay Area",  # Add the map title within tm_layout
    main.title.size = 0.7,
    main.title.position = "center"
  )


map_electric <- tm_basemap() +
  tm_shape(overlay_listings_zip) +
  tm_fill(col = "Combined_Avg_Comm_Rate", palette = "YlGn", n = 5, style = "fisher", title = "Commerical Electrical Rate") +
  tm_compass(position = c(0.85, 0.88)) +
  tm_scalebar(position = c(0.58, 0.02)) +
  tm_borders(col = "black", lwd = 0.2)+
  tm_layout(
    legend.title.size = 0.7,
    legend.text.size = 0.6, 
    inner.margins = c(0.185, 0.13, 0.02, 0.08), 
    legend.position = c(0.05, 0.02), 
    legend.width = 0.5, 
    bg.color = "#eaf5fa",
    main.title = "Commercial Electrical Rate in Bay Area",  # Add the map title within tm_layout
    main.title.size = 0.7,
    main.title.position = "center"
  )

# Chunk 52
print(tmap_arrange(map_cluster1, map_veh, map_electric, nrow =1, ncol= 3))

save.image(file=".r_session_state.RData")
