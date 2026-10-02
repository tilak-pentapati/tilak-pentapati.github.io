load(".r_session_state.RData")
for(i in 1:length(packages_required)) {
  suppressPackageStartupMessages(library(packages_required[i], character.only = T))
}
readRenviron("~/.Renviron") # load the API key into this session

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
print("Head of us_zip:")
print(head(us_zip))

save.image(file=".r_session_state.RData")
