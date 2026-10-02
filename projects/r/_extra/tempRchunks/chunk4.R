load(".r_session_state.RData")
for(i in 1:length(packages_required)) {
  suppressPackageStartupMessages(library(packages_required[i], character.only = T))
}

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

print("electric_u summary:")
print(summary(electric_u))
print("vc_califorina summary:")
print(summary(vc_califorina))

save.image(file=".r_session_state.RData")
