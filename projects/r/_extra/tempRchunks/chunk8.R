load(".r_session_state.RData")
for(i in 1:length(packages_required)) {
  suppressPackageStartupMessages(library(packages_required[i], character.only = T))
}

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
print(moran.mc(bayarea_zip_vars1_NoZero$estimate_MedianHhI, w_queen_std1, nsim=100, alternative="greater")) # Moran's I statistic for Median Household Income Variable in Bay Area (nsim=100 for speed)

# Chunk 37
print(moran.mc(bayarea_zip_vars_NoZero$estimate_InternetUse, w_queen_std2, nsim=100, alternative="greater")) # Moran's I statistic for Interent Use Variable in Bay Area (nsim=100 for speed)

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
cat("\n")

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
  tm_scalebar(position = c(0.625, 0.02)) +
  tm_add_legend(type = "fill", col = c("royalblue2", "red2", "darkgreen", "gold", "lightgrey"), 
                labels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant"), title = "LISA cluster(Median Household Income)") +
  tm_layout(
    legend.title.size = 0.715,
    legend.text.size = 0.6, 
    inner.margins = c(0.22, 0.13, 0.02, 0.08), 
    legend.position = c(0.05, 0.02), 
    legend.width = 0.5, 
    bg.color = "#eaf5fa",
    main.title = "LISA Clusters for Median Household Income Variable", 
    main.title.size = 0.7,
    main.title.position = "center"
  )


lisa_map_cluster2 <- borders2 +
  hh_map2 + ll_map2 + hl_map2   + ns_map2 +  #  Combining All quadrant maps
  tm_compass(position = c(0.85, 0.85)) +
  tm_scalebar(position = c(0.625, 0.02)) +
  tm_add_legend(type = "fill", col = c("royalblue2", "red2", "darkgreen", "lightgrey"), 
                labels = c("High-High", "Low-Low", "High-Low", "Not significant"), title = "LISA cluster(Broadband Internet Use)") +
  tm_layout(
    legend.title.size = 0.715,
    legend.text.size = 0.6, 
    inner.margins = c(0.185, 0.13, 0.02, 0.08), 
    legend.position = c(0.05, 0.02), 
    legend.width = 0.5, 
    bg.color = "#eaf5fa",
    main.title = "LISA Clusters for Broadband Internet Use Variable",  
    main.title.size = 0.7,
    main.title.position = "center"
  )

# Chunk 45
print(tmap_arrange(lisa_map_cluster1, lisa_map_cluster2))

save.image(file=".r_session_state.RData")
