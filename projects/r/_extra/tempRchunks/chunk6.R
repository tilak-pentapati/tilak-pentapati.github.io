load(".r_session_state.RData")
for(i in 1:length(packages_required)) {
  suppressPackageStartupMessages(library(packages_required[i], character.only = T))
}

# Chunk 18
tmap_mode("view")   # Interactive mode
map1.11 <- tm_shape(listings_bayarea_zip1) +  # listings_bayarea_zip1 is just re-arrangement of listings_bayareaa_zip with first column as count_airbnb
  tm_fill( "count_airbnb",style="fixed", title = "No. of Airbnb's per ZIP", alpha = 1, breaks=c( 1, 72, 172, 359, 500, Inf), palette = "-viridis")+ 
  tm_borders()+
  tm_text("ZIP", size = 0.3)+
  tm_view(set.view = c(center_lon, center_lat, zoom_level))

map1.12 <- tm_shape(listings_bayarea_zip2) +  # listings_bayarea_zip2 is just re-arrangement of listings_bayareaa_zip with first column as Mean Price
  tm_fill( "mean_price",style="fixed", title = "Mean Price of Airbnb per ZIP", alpha = 1, breaks=c( 74, 155, 287, 501,1000, 1800, Inf), palette = "-viridis")+
  tm_borders()+
  tm_text("ZIP", size = 0.3)+
  tm_view(set.view = c(center_lon, center_lat, zoom_level))

# Chunk 19
print(map1.11)

# Chunk 20
print(map1.12)

# Chunk 21
bay_area_zipcodes <- bayarea_zipcodes$ZIP   # Creating new var with Bayarea ZIP codes

bayarea_zip_vars_data <- us_zip %>%        # Filtering ACS Data with Two variables, by Bay Area ZIP codes
  filter(GEOID %in% bay_area_zipcodes)%>%  # Using Bay Area ZIP codes variable, created above
  rename(ZIP = GEOID)                     # renaming GEOID with ZIP

bayarea_zip_vars <- bayarea_zipcodes %>% left_join(bayarea_zip_vars_data, by = "ZIP") %>% filter(!(is.na(estimate_InternetUse))) # Joining two vars from ACS with bayarea_zipcodes shapefile which is downloaded from berkeley library and also removing NA values of Internet Use variable.

# Chunk 22
print("Head of bayarea_zip_vars:")
print(head(bayarea_zip_vars))

# Chunk 23
bayarea_zip_vars1 <- bayarea_zip_vars %>%
  filter(!(is.na(estimate_MedianHhI)))  #  Removing rows of bayarea_zip_vars where 'estimate_MedianHhI' is NA

# Chunk 24
cat("Summary of Median Household Income per Zipcode (estimate_MedianHhI): ")
cat("\n")
print(summary(bayarea_zip_vars1$estimate_MedianHhI))   # Different Standard statistics for estimate_MedianHhI
cat("\n")
cat("\nSummary of No.of People aged 18-64 who has Computer and Internet Subscription ")
cat("per Zipcode (estimate_InternetUse): ")
cat("\n")
print(summary(bayarea_zip_vars$estimate_InternetUse))   # Different Standard statistics for estimate_MedianHhI
cat("\n")
# Calculate standard deviation
std_dev_MedianHhI<- sd(bayarea_zip_vars1$estimate_MedianHhI, na.rm = TRUE)
std_dev_iu<- sd(bayarea_zip_vars$estimate_InternetUse, na.rm = TRUE)
# Display the result
cat("Standard Deviation of Variable(estimate_MedianHhI):    ", std_dev_MedianHhI, "\n")
cat("\n")
cat("Standard Deviation of Variable(estimate_InternetUse):  ", std_dev_iu, "\n")


save.image(file=".r_session_state.RData")
