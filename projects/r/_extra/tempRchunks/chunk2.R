load(".r_session_state.RData")
for(i in 1:length(packages_required)) {
  suppressPackageStartupMessages(library(packages_required[i], character.only = T))
}

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
census_api_key(J_API_KEY , install = TRUE, overwrite = TRUE)

# Chunk 6
options(tigris_use_cache = TRUE)

if (nzchar(Sys.getenv("CENSUS_API_KEY"))) {
  message_text <-  "API key to get Census data is successfully set up."
} else {
  message_text <- "Census API key setup failed."
}

# Chunk 7
print(message_text)

save.image(file=".r_session_state.RData")
