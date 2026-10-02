# Chunk 1
packages_required <- c("tidycensus", "tigris","tidyverse", "data.table","sf", "tmap", "tmaptools","readr", "geojsonsf", "osmdata","devtools","logger" ,"RColorBrewer","classInt","R.utils","dplyr", "ggplot2", "viridis","raster","terra","exactextractr", "tidyterra", "tibble", "rosm", "spdep", "cluster", "GGally", "tidyr", "gridExtra","patchwork") 
# Include any additional packages needed, in above variable .
 
new_req_packages <- packages_required[!(packages_required %in% installed.packages()[,"Package"])] # checking which packages are already installed and filtering out installed packages, new variable is assigned for names of uninstalled packages.
if(length(new_req_packages)) install.packages(new_req_packages) # packages(not installed yet) are being installed. 

for(i in 1:length(packages_required)) {
  library(packages_required[i], character.only = T)        # loading all packages named in packages_required
}
# Below message is differently coded output for checking installed packages, not the output of above line

# Chunk 2
# Checking if all required packages are installed and loaded
if(all(packages_required %in% installed.packages()[,"Package"])) {
   message_1 <- "All required packages for Assignment-2 are installed and loaded properly."
} else {
   message_1 <- "Some packages could not be installed or loaded."
}

# Chunk 3
print(message_1)

# Chunk 4
J_API_KEY <- "77a4dfc5c8bcf14004b7c739aa036700d2fd4b2d"   # this is api key, specific to each individual

# Chunk 5
census_api_key(J_API_KEY , install = TRUE, overwrite = TRUE) # Please replace 'J_API_KEY' with API Key obtained from census bureau
# Below message is differently coded output for checking system environment variable(CENSUS_API_KEY), not the output of above line

# Chunk 6
options(tigris_use_cache = TRUE) # setting this option, helps in a way that download shapefiles will be in our Local R environment.

# Checking if the API key is successfully set up
if (nzchar(Sys.getenv("CENSUS_API_KEY"))) {
  message_text <-  "API key to get Census data is successfully set up."
} else {
  message_text <- "Census API key setup failed."
}

# Chunk 7
print(message_text)

# Chunk 8
median_hh_income <- get_acs(  # get_acs is a function, we are using to get ACS Suvrey data from US Census Bureau
  geography = "zcta",         # we are specifying geography unit for variable we getting, like, ZCTA, Tract, County, State etc.
  variables = "B19013_001E",  # we are geeting MEDIAN HOUSEHOLD INCOME, "B19013_001E" is variable name in ACS
  year = 2020,      # we are getting census data collected from 2016-2020
  geometry = FALSE  # we use geometry of bay area zipcode shapefile from berkeley library, not ACS Data
) %>%
  rename(estimate_MedianHhI = estimate) %>%   # we are renaming retrieved variable(estimate to appropriate name ( estimate_Median_HhI)
  select(GEOID, estimate_MedianHhI)           # we are selecting only GEOID(ZCTA's CODE, also ZIP CODE), estimate Value (ignoring MOE)

# Retrieving Internet Use Variable
internet_use <- get_acs(      # get_acs is a function, we are using to get ACS Suvrey data from US Census Bureau 
  geography = "zcta",         # we are specifying geography unit for variable we getting, like, ZCTA, Tract, County, State etc.
  variables = "B28005_011E",  # we are geeting No.of People aged 18-64, who has computer and Broadband Internet Subscription, "B28005_011E" is variable name in ACS
  year = 2020,                # we are getting census data collected from 2016-2020
  geometry = FALSE            #  we use geometry of bay area zipcode shapefile from berkeley library, not ACS Data
) %>%
  rename(estimate_InternetUse = estimate) %>%  # we are renaming retrieved variable(estimate to appropriate name ( estimate_Internet_Use)
  select(GEOID, estimate_InternetUse)          # we are selecting only GEOID(ZCTA's CODE, also ZIP CODE), estimate Value (ignoring MOE)

us_zip <- left_join(median_hh_income, internet_use, by = "GEOID") # simply join two datasets by GEOID to get Dataset in req format  

# Chunk 9
head(us_zip)  #showing first six rows of Dataset

# Chunk 10
# Loading airbnb Listings downloaded from AIRBNB website for areas(san francisco, oakland and san mateo)
listings_oakland <- read.csv("/Users/jay/Desktop/websiteV02/data/data_python/airbnb_sf/tables/Oakland/listings.csv") 


listings_sanfrancisco <- read.csv("/Users/jay/Desktop/websiteV02/data/data_python/airbnb_sf/tables/SanFrancisco/listings.csv")
listings_sanmateo <- read.csv("/Users/jay/Desktop/websiteV02/data/data_python/airbnb_sf/tables/SanMateo/listings.csv")

# Row Binding all three datasets to produce specified bay area's Airbnb Listings
listings_bayarea <- rbind(listings_sanfrancisco, listings_sanmateo, listings_oakland)

# Loading Bay area ZCTA's or ZIPCODE Areas for geometry which are downloaded from Berkeley Library
bayarea_zipcodes <- read_sf("/Users/jay/Desktop/websiteV02/data/data_python/airbnb_sf/bayarea_zipcodes/bayarea_zipcodes.shp")

# Chunk 11
# We are using Commerical Electricity Rate Dataset for Both investor-owned and non-investor owned utilities per ZIP Code
iou_data <- read.csv("/Users/jay/Desktop/websiteV02/data/data_python/airbnb_sf/tables/iou_zipcodes_2020.csv")
non_iou_data <- read.csv("/Users/jay/Desktop/websiteV02/data/data_python/airbnb_sf/tables/non_iou_zipcodes_2020.csv")

# We are using No.of Vehicles per ZIP Code dataset
vehicle_count_2020 <- read.csv("/Users/jay/Desktop/websiteV02/data/data_python/airbnb_sf/tables/vehicle-count-2020.csv")

# Chunk 12
non_iou_data_c <- non_iou_data %>%  # Data cleaning process for additional datasets, 
  select(zip, comm_rate) %>%       # we select only zip, commerical rate from non-investor based utilities
  group_by(zip) %>%                # we group them by ZIP
  summarise(
    Avg_Non_IOU_Comm_Rate = mean(comm_rate, na.rm = TRUE)  # Mean of different Non-Investor Based utilities per ZIP
  )
  

iou_data_c <- iou_data %>%
  select(zip, comm_rate) %>%      # we select only zip, commerical rate from Investor based utilities
  group_by(zip) %>%               # we group them by ZIP
  summarise(
    Avg_IOU_Comm_Rate = mean(comm_rate, na.rm = TRUE)  # Mean of different Investor Based utilities per ZIP
  )

electric_u <- non_iou_data_c %>%
  full_join(iou_data_c, by = "zip")%>%   # Join both Investor-owned and Non-Investor utility datasets by zipcode
  mutate(
    Avg_IOU_Comm_Rate = replace_na(Avg_IOU_Comm_Rate, 0),    # replacing any NA values with 0
    Avg_Non_IOU_Comm_Rate = replace_na(Avg_Non_IOU_Comm_Rate, 0 ),  # replacing any NA values with 0
    Combined_Avg_Comm_Rate = (Avg_Non_IOU_Comm_Rate + Avg_IOU_Comm_Rate ) / 2  # using average as statistic for combined Commerical electricity Rate
  )%>%
  select(zip, Combined_Avg_Comm_Rate) %>% rename(ZIP = zip) %>% mutate(ZIP = as.character(ZIP)) # selecting only ZIP and New Statistic, ignoring rest

vc_califorina <- vehicle_count_2020 %>%
  select(-Date, -Model.Year, -Fuel, -Make) %>%     # we are only selecting Zip.code, Vehicles and Duty cols, and ignoring rest
  filter(Duty == "Light" & Zip.Code != "OOS") %>%  # we are filter only Light-Duty vehicles and removing Out of Service(OOS) locations from data
  group_by(Zip.Code) %>%                           # group by ZIP Code
  summarise(No_of_Light_Duty_Vehicles = sum(Vehicles)) %>% rename( ZIP = Zip.Code)  # summing up Light-Duty Vehicles per ZIP and then rename it

# Chunk 13
listings_bayarea_sf <- listings_bayarea %>%
  st_as_sf(coords = c("longitude", "latitude")) %>% # creating point shapefile from coordinates
  st_set_crs(4326)                                  # setting CRS to geodetic system of earth (4326, WGS84) as they are long, lat coords

listings_bayarea_sf_transform <- st_transform(listings_bayarea_sf, crs = 2227)   # Now, Changing CRS to mentioned Project CRS
overlay_listings <- st_join(bayarea_zipcodes, listings_bayarea_sf_transform)     # Performing Spatial Join to create overlay Shapefile of Airbnb Listings and Bay area ZIP shapefile

st_crs(overlay_listings) == st_crs(listings_bayarea_sf_transform)                # Verifying CRS of New shapefile

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
summary(listings_bayarea_zip1$count_airbnb)

cat("\nSummary of Mean Price per Zipcode Variable(mean_price): ")
cat("\n")
summary(listings_bayarea_zip1$mean_price)
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
map1.11

# Chunk 20
map1.12

# Chunk 21
bay_area_zipcodes <- bayarea_zipcodes$ZIP   # Creating new var with Bayarea ZIP codes

bayarea_zip_vars_data <- us_zip %>%        # Filtering ACS Data with Two variables, by Bay Area ZIP codes
  filter(GEOID %in% bay_area_zipcodes)%>%  # Using Bay Area ZIP codes variable, created above
  rename(ZIP = GEOID)                     # renaming GEOID with ZIP

bayarea_zip_vars <- bayarea_zipcodes %>% left_join(bayarea_zip_vars_data, by = "ZIP") %>% filter(!(is.na(estimate_InternetUse))) # Joining two vars from ACS with bayarea_zipcodes shapefile which is downloaded from berkeley library and also removing NA values of Internet Use variable.

# Chunk 22
head(bayarea_zip_vars)

# Chunk 23
bayarea_zip_vars1 <- bayarea_zip_vars %>%
  filter(!(is.na(estimate_MedianHhI)))  #  Removing rows of bayarea_zip_vars where 'estimate_MedianHhI' is NA

# Chunk 24
cat("Summary of Median Household Income per Zipcode (estimate_MedianHhI): ")
cat("\n")
summary(bayarea_zip_vars1$estimate_MedianHhI)   # Different Standard statistics for estimate_MedianHhI
cat("\n")
cat("\nSummary of No.of People aged 18-64 who has Computer and Internet Subscription ")
cat("per Zipcode (estimate_InternetUse): ")
cat("\n")
summary(bayarea_zip_vars$estimate_InternetUse)   # Different Standard statistics for estimate_MedianHhI
cat("\n")
# Calculate standard deviation
std_dev_MedianHhI<- sd(bayarea_zip_vars1$estimate_MedianHhI, na.rm = TRUE)
std_dev_iu<- sd(bayarea_zip_vars$estimate_InternetUse, na.rm = TRUE)
# Display the result
cat("Standard Deviation of Variable(estimate_MedianHhI):    ", std_dev_MedianHhI, "\n")
cat("\n")
cat("Standard Deviation of Variable(estimate_InternetUse):  ", std_dev_iu, "\n")

# Chunk 25
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

# Chunk 26
tmap_arrange(map2.1, map2.2)

# Chunk 27
bayarea_zip_vars_map3 <- bayarea_zip_vars %>%
  select(estimate_InternetUse, everything())  # Changing the order of attributes for interactive map

# Chunk 28
listings_bayarea_sf_transform_copy <- listings_bayarea_sf_transform
listings_bayarea_sf_transform_copy$nl_price <- log(listings_bayarea_sf_transform_copy$price) # calculating natural logarithm of price ( log to base e)
st_crs(bayarea_zip_vars_map3) == st_crs(listings_bayarea_sf_transform_copy) # checking CRS of both shapefiles

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

map3

# Chunk 31
nb_q1 <- poly2nb(bayarea_zip_vars1, queen = TRUE) # shapefile used, have no NA values of Median Household Income variable

w_queen1 <- nb2listw(nb_q1, style = "B", zero.policy=TRUE) # queen-contiguity based Weights for list of neighbors for all observations 
isolates1 <- which(w_queen1$neighbours == "0")
bayarea_zip_vars1_NoZero <- bayarea_zip_vars1[-c(isolates1),] # shapefile where there are no observation with zero neighbours based on queen-contiguity

# Chunk 32
nb_q2 <- poly2nb(bayarea_zip_vars, queen = TRUE) # shapefile used, have no NA values of Broadband Internet Use

w_queen2 <- nb2listw(nb_q2, style = "B", zero.policy=TRUE)
isolates2 <- which(w_queen2$neighbours == "0")

bayarea_zip_vars_NoZero <- bayarea_zip_vars[-c(isolates2),] # shapefile where there are no observation with zero neighbours based on queen-contiguity

# Chunk 33
nb_q1 <- poly2nb(bayarea_zip_vars1_NoZero, queen = TRUE) # Constructing neighbours list from filtering observations which have no zero neighbours and no NA values of Median Household Income
w_queen_std1 <- nb2listw(nb_q1, style = "W") # creating spatial weights matrix using queen contiguity and row-standardardised weights

nb_q2 <- poly2nb(bayarea_zip_vars_NoZero, queen = TRUE) # Constructing neighbours list from filtering out observations which have zero neighbours and no NA values of Broadband Internet Use
w_queen_std2 <- nb2listw(nb_q2, style = "W") # creating spatial weights matrix using queen contiguity and row-standardardised weights

# Chunk 34
bayarea_zip_vars1_NoZero$sl_MedianHhI <- lag.listw(w_queen_std1, bayarea_zip_vars1_NoZero$estimate_MedianHhI) # calculating spatial lag for Medain Household Income Variable in Bay Area 

bayarea_zip_vars_NoZero$sl_InternetUse <-  lag.listw(w_queen_std2, bayarea_zip_vars_NoZero$estimate_InternetUse) # calculating spatial lag for Internet Use Variable in Bay Area 

# Chunk 35
bayarea_zip_vars1_NoZero$MedianHhI_std <- (bayarea_zip_vars1_NoZero$estimate_MedianHhI - mean(bayarea_zip_vars1_NoZero$estimate_MedianHhI))/sd(bayarea_zip_vars1_NoZero$estimate_MedianHhI)   
# Spatial Lag for Median Household Income using standardized weights 
bayarea_zip_vars1_NoZero$sl_Median_HhI_std <- lag.listw(w_queen_std1, bayarea_zip_vars1_NoZero$MedianHhI_std)

bayarea_zip_vars_NoZero$InternetUse_std <- (bayarea_zip_vars_NoZero$estimate_InternetUse - mean(bayarea_zip_vars_NoZero$estimate_InternetUse))/sd(bayarea_zip_vars_NoZero$estimate_InternetUse)
# Spatial Lag for Internet Use variable using standardized weights 
bayarea_zip_vars_NoZero$sl_InternetUse_std <- lag.listw(w_queen_std2, bayarea_zip_vars_NoZero$InternetUse_std)

# Chunk 36
moran.mc(bayarea_zip_vars1_NoZero$estimate_MedianHhI, w_queen_std1, nsim=1000, alternative="greater") # Moran's I statistic for Median Household Income Variable in Bay Area

# Chunk 37
moran.mc(bayarea_zip_vars_NoZero$estimate_InternetUse, w_queen_std2, nsim=1000, alternative="greater") # Moran's I statistic for Interent Use Variable in Bay Area

# Chunk 38
lisa_perm1 <- localmoran_perm(bayarea_zip_vars1_NoZero$estimate_MedianHhI, w_queen_std1, nsim=1000, alternative="two.sided") # Lisa Clustering for Median Household Income in Bay Area

lisa_perm2 <- localmoran_perm(bayarea_zip_vars_NoZero$estimate_InternetUse, w_queen_std2, nsim=1000, alternative="two.sided") # Lisa Clustering for Internet Use Variable in Bay Area

# Chunk 39
quadrants1 <- hotspot(lisa_perm1, Prname="Pr(z != E(Ii)) Sim", cutoff=0.2) # Creating quadrants from significance values of local Moran's I statistic calculated from Median Household Income

quadrants2 <- hotspot(lisa_perm2, Prname="Pr(z != E(Ii)) Sim", cutoff=0.2) # Creating quadrants from significance values of local Moran's I statistic calculated from Internet Use Variable

# Chunk 40
cat("Levels of Quadrant 1 :\n")
cat(levels(quadrants1), sep = ", ")  # Significant Levels of quadrant 1
cat("\n\nLevels of Quadrant 2 :\n")
cat(levels(quadrants2), sep = ", ") # Significant Levels of Quadrant 2

# Chunk 41
bayarea_zip_vars1_NoZero$quadrant1 <- as.character(quadrants1)  %>% replace_na("Not significant") # Adding created quadrants to shapefile(no NA values for Median Household Income) where no observation with zero neighbours based on queen-contiguity


bayarea_zip_vars_NoZero$quadrant2 <- as.character(quadrants2)  %>% replace_na("Not significant") # Adding created quadrants to shapefile(no NA values for Internet Use Variable) where no observation with zero neighbours based on queen-contiguity

# Chunk 42
borders1 <- tm_shape(bayarea_zip_vars1_NoZero) + 
  tm_fill() +
  tm_borders(col = "black", lwd = 0.2)

borders2 <- tm_shape(bayarea_zip_vars_NoZero) + 
  tm_fill() +
  tm_borders(col = "black", lwd = 0.2)

# Chunk 43
hh1 <- bayarea_zip_vars1_NoZero %>% dplyr::filter(quadrant1 == "High-High") # filtering only 'High-High' quadrant to map
hh_map1 <- tm_shape(hh1) +  
  tm_fill(col = "royalblue2", alpha=0.8)+
  tm_borders(col = "black", lwd = 0.3)

ll1 <- bayarea_zip_vars1_NoZero %>% dplyr::filter(quadrant1 == "Low-Low") # filtering only 'High-High' quadrant to map
ll_map1 <- tm_shape(ll1) +  
  tm_fill(col = "red2", alpha=0.8)+
  tm_borders(col = "black", lwd = 0.3)

lh1 <- bayarea_zip_vars1_NoZero %>% dplyr::filter(quadrant1 == "Low-High") # filtering only 'Low-High' quadrant to map
lh_map1 <- tm_shape(lh1) +   
  tm_fill(col = "gold", alpha=0.8)+
  tm_borders(col = "black", lwd = 0.3)

hl1 <- bayarea_zip_vars1_NoZero %>% dplyr::filter(quadrant1 == "High-Low") # filtering only 'High-Low' quadrant to map
hl_map1 <- tm_shape(hl1) +  
  tm_fill(col = "darkgreen", alpha=0.8)+
  tm_borders(col = "black", lwd = 0.3)

ns1 <- bayarea_zip_vars1_NoZero %>% dplyr::filter(quadrant1 == "Not significant") # filtering only 'Not significant' quadrant to map
ns_map1 <- tm_shape(ns1) +  
  tm_fill(col = "lightgrey", alpha=0.8)+
  tm_borders(col = "black", lwd = 0.3)

hh2 <- bayarea_zip_vars_NoZero %>% dplyr::filter(quadrant2 == "High-High")   # filtering only 'High-High' quadrant to map
hh_map2 <- tm_shape(hh2) +  
  tm_fill(col = "royalblue2", alpha=0.8)+
  tm_borders(col = "black", lwd = 0.3)

ll2 <- bayarea_zip_vars_NoZero %>% dplyr::filter(quadrant2 == "Low-Low")  # filtering only 'Low-Low' quadrant to map
ll_map2 <- tm_shape(ll2) +  
  tm_fill(col = "red2", alpha=0.8)+
  tm_borders(col = "black", lwd = 0.3)

hl2 <- bayarea_zip_vars_NoZero %>% dplyr::filter(quadrant2 == "High-Low") # filtering only 'High-Low' quadrant to map
hl_map2 <- tm_shape(hl2) +  
  tm_fill(col = "darkgreen", alpha=0.8)+
  tm_borders(col = "black", lwd = 0.3)

ns2 <- bayarea_zip_vars_NoZero %>% dplyr::filter(quadrant2 == "Not significant") # filtering only 'Not significant quadrant to map
ns_map2 <- tm_shape(ns2) +  
  tm_fill(col = "lightgrey", alpha=0.8)+
  tm_borders(col = "black", lwd = 0.3)

# Chunk 44
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

# Chunk 45
tmap_arrange(lisa_map_cluster1, lisa_map_cluster2)

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
k6cls$centers # centres of 6 clusters, six 8-tuple vectors defines centre of each cluster

# Chunk 51
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

# Chunk 52
tmap_arrange(map_cluster1, map_veh, map_electric, nrow =1, ncol= 3)
