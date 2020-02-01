library(tidyverse)

setwd("~/R/PROJECTS_MY/CommunityAssemblages/Results/nSites0.2_maxCombo3_At0.6_Bt0.25")

sc <- readRDS("Indicators4Clusters.RDS")
indCombs <- readRDS("IndicatorCombos.RDS")

indCombs <- as.vector(c("LT", "AM", "AC+LT", "AC+NT"))

# Coverage is equal to the cumulative sum of sensitivity (Bt)

# Pseudocode:
# 1. Build df of Bt for each cluster
# 2. Build ?vector? of indicator type: IndType | Count
# 3. Final df: Cluster | Indicator | Indicator Type | Sensitivity (B)
# 3. Build histogram of Indicator Type + Sensitivity by Cluster

# Select one element with a list of lists (nested list element selection)
sens <- lapply(sc, '[',7)
# indNmes <- sapply(indCombs, "[",1:4)

# Build empty dataframe
i = 4
indCov <- tibble("Cluster" = as.integer(),
             "indSp" = as.character(),
             "typeCnt" = as.character(),
             "B" = as.numeric())
new.row <- indCov

# Loop through each cluster and populate
for (i in 4:length(sc)){
  # Create a list of vectors for each new column for this cluster
  col_vectors <- list(
    indSp = indCombs,
    B = (100*sc[[i]]$A$stat) )
  # Bind the list of vectors into new rows in a dataframe
  new.row <- bind_rows(col_vectors)
  # Add Cluster number
  new.row$Cluster <- i
  # Calculate Indicator type
  new.row$typeCnt <- as.character(1 + (str_count(new.row$indSp, pattern='\\+')) )
  new.row$typeCnt <- paste0(new.row$typeCnt, " Species" )
  indCov <- bind_rows(indCov, new.row)
}
indCov
#ind.Type <- data.frame("typeCnt"=1:3,"Groupings"=c("Singletons", ""))

# Nope
#-----
# Group by Cluster
# Order by B
# Calculate cummulative value of B for each cluster and typeCnt
# Subtract typeCnt 1 - 2
# Subtract typeCnt 

# Try
#----
# Calculate coverage from 1 species
# Calculate coverage from 2 species
# I give up!

indCov$cum <- 

# Build histogram
#----------------
# Change y axis to read "Coverage %"
# Limit y axis to 100
# type in increasing order (1,2,3)
# Save as png


#indCov <- as.data.frame(new.row)

ggplot(indCov, aes(fill=typeCnt, y=B, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")
