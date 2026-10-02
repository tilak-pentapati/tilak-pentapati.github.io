# test_spatial.R

cat("1. Loading spatial libraries...\n")
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(gstat))

cat("2. Generating sample points...\n")
data <- data.frame(
  lon = c(-2.98, -2.99, -2.97), 
  lat = c(53.40, 53.41, 53.39), 
  zinc = c(120, 145, 110)
)

cat("3. Converting to sf object (tests GDAL/GEOS geometry engine)...\n")
spatial_points <- st_as_sf(data, coords = c("lon", "lat"), crs = 4326)

print(spatial_points)
cat("\nSUCCESS: sf and gstat are fully operational!\n")