load(".r_session_state.RData")
for(i in 1:length(packages_required)) {
  suppressPackageStartupMessages(library(packages_required[i], character.only = T))
}

# Chunk 13
listings_bayarea_sf <- listings_bayarea %>%
  st_as_sf(coords = c("longitude", "latitude")) %>% # creating point shapefile from coordinates
  st_set_crs(4326)                                  # setting CRS to geodetic system of earth (4326, WGS84) as they are long, lat coords

listings_bayarea_sf_transform <- st_transform(listings_bayarea_sf, crs = 2227)   # Now, Changing CRS to mentioned Project CRS
overlay_listings <- st_join(bayarea_zipcodes, listings_bayarea_sf_transform)     # Performing Spatial Join to create overlay Shapefile of Airbnb Listings and Bay area ZIP shapefile

print("CRS Match Check:")
print(st_crs(overlay_listings) == st_crs(listings_bayarea_sf_transform))                # Verifying CRS of New shapefile

# Chunk 14
listings_bayarea_zip <- overlay_listings  %>%  # creating new dataframe with no.of Airbnb's per ZIP and Mean Price per ZIP
  group_by(ZIP) %>%                     # group by ZIP
  summarise(
    count_airbnb = sum(!is.na(price)),      # counting non-NA price rows for listings
    mean_price = mean(price, na.rm = TRUE)  # calculate mean, ignoring NAs
  )%>%
  filter(count_airbnb != 0) 

# Chunk 15
listings_bayarea_zip1 <- listings_bayarea_zip %>%
  select(count_airbnb, everything())  # Changing the order of attributes for interactive map

listings_bayarea_zip2 <- listings_bayarea_zip %>%
  select(mean_price, everything())%>% # Changing the order of attributes for interactive map
  mutate(mean_price = round(mean_price, 1))  # rounding mean price to one decimal place.

# Chunk 16
cat("Summary of Airbnb's Count per Zipcode Variable(count_airbnb): ")
cat("\n")
print(summary(listings_bayarea_zip1$count_airbnb))

cat("\nSummary of Mean Price per Zipcode Variable(mean_price): ")
cat("\n")
print(summary(listings_bayarea_zip1$mean_price))
cat("\n")
# Calculate standard deviation
std_dev_count_airbnb <- sd(listings_bayarea_zip1$count_airbnb, na.rm = TRUE)
std_dev_mean_price<- sd(listings_bayarea_zip1$mean_price, na.rm = TRUE)
# Display the result
cat("Standard Deviation of Variable(count_airbnb):", std_dev_count_airbnb, "\n")
cat("Standard Deviation of Variable(mean_price):  ", std_dev_mean_price, "\n")

# Chunk 17
center_lat <- 37.703253  # Replace with the latitude of your area's center
center_lon <- -122.310913  # Replace with the longitude of your area's center
zoom_level <- 10 

save.image(file=".r_session_state.RData")
