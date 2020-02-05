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
    A = (100*sc[[i]]$A$stat) )
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
# Calculate coverage from 1 species only
# Calculate coverage from 2 species only

# What elements do we need for coverage()?
# Can they be subset using [,]?

# New strategy
#-------------
# Determine Ind.Type for each Indicator
# Build "melt" table:
# Site | Indicator | Ind.Type
# Determine which sites can be identified with more than one indicator type
# Reassign those sites to smallest Ind.Type
# Calculate frequency/coverage for each cluster and Ind.Type
# Plot!

# Build histogram
#----------------
# Change y axis to read "Coverage %"
# Limit y axis to 100
# type in increasing order (1,2,3)
# Save as png


#indCov <- as.data.frame(new.row)

ggplot(indCov, aes(fill=typeCnt, y=B, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")


# Change in coverage and number of species per indicator
#-------------------------------------------------------
# Select one element with a list of lists (nested list element selection)
sens <- lapply(sc, '[',7)
# indNmes <- sapply(indCombs, "[",1:4)

# Build empty dataframe
i = 4
#indCombs <- as.vector(c("LT", "AM", "AC+LT", "AC+NT"))

# indCov <- tibble("Cluster" = as.integer(),
#                  "indSp" = as.character(),
#                  "typeCnt" = as.numeric(),
#                  "B" = as.numeric())
# new.row <- indCov

new.row <- tibble("Cluster" = as.integer(),
                  "indSp" = as.character(),
                  "typeCnt" = as.numeric(),
                  "B" = as.numeric())
indCov <- vector("list",4)


# Loop through each cluster and populate
for (i in 1:length(sc)){
  # Create a list of vectors for each new column for this cluster
  col_vectors <- list(
    indSp = row.names(indCombs[[i]]),
    B = (100*sc[[i]]$B$stat) )
  # Bind the list of vectors into new rows in a dataframe
  new.row <- bind_rows(col_vectors)
  # Add Cluster number
  new.row$Cluster <- i
  # Calculate Indicator type
  new.row$typeCnt <- (1 + (str_count(new.row$indSp, pattern='\\+')) )
  #new.row$typeCnt <- paste0(new.row$typeCnt, " Species" )
  #indCov <- bind_rows(indCov, new.row)
  indCov[[i]] <- new.row
}
indCov




covAll <- vector("list",4)
#coverageSD <- coverage(sppComb, indvalspcomb, At=0.6, Bt=Bt, alpha = 0.05)
# 1         2         3         4 
# 0.7374582 0.4233171 0.9925373 0.6380449 

for (i in 1:length(sc)){
  # Build logic vectors
  df <- indCov[[i]]
  df <- df %>% 
    mutate(ind.1 = ifelse(df$typeCnt==1, TRUE, FALSE)) %>%
    mutate(ind.2 = ifelse(df$typeCnt==2, TRUE, FALSE)) %>%
    mutate(ind.3 = ifelse(df$typeCnt==3, TRUE, FALSE))
  covAll[[i]]$combo1 <- coverage(sc[[i]], indvalspcomb, selection=as.vector(df$ind.1), At=0.6, Bt=Bt, alpha = 0.05, type="lowerCI")
  covAll[[i]]$combo2 <- coverage(sc[[i]], indvalspcomb, selection=as.vector(df$ind.2), At=0.6, Bt=Bt, alpha = 0.05, type="lowerCI")
  covAll[[i]]$combo3 <- coverage(sc[[i]], indvalspcomb, selection=as.vector(df$ind.3), At=0.6, Bt=Bt, alpha = 0.05, type="lowerCI")
  covAll[[i]]$Sum <- covAll[[i]]$combo3 + (covAll[[i]]$combo3 - covAll[[i]]$combo2) + (covAll[[i]]$combo2 - covAll[[i]]$combo3)
}
# Coverage - proportion of sites where one or another indicator is found
coverageSC <- coverage(sppComb,indvalspcomb)
coverageSD <- coverage(sppComb, indvalspcomb, At=0.6, Bt=Bt, alpha = 0.05)




plotcoverage(sc[[i]], main=label)
plotcoverage(sc[[i]], max.order=1, add=TRUE, lty=2, col="red", alpha = 0.05)
legend(x=0.1, y=20, legend=c("Species combinations","Species singletons"),
       lty=c(1,2), col=c("black","red"), bty="n")


plot(sc[[i]], type="A", main=label)
plot(sc[[i]], type="B")


coverageSC <- coverage(sppComb,indvalspcomb)
coverageSD <- coverage(sppComb, indvalspcomb, At=0.6, Bt=Bt, alpha = 0.05)
covIndicators <- coverage(sc[[2]], indvalspcomb, selection=sc[[2]]$group.vec) #, At=0.6, Bt=Bt, alpha = 0.05)


plotcoverage(indCombs)
plotcoverage(indCombs, max.order=2, add=TRUE, lty=2)


# Subset sppComb by species indicator vector
for (i in 1:length(indCov)){
  sel <- indCov[[i]]$indSp
  newX <- sppComb[,sel]
  coverageX[[i]] <- coverage(sppComb, indvalspcomb, selection=sel) #, At=0.6, Bt=Bt, alpha = 0.05)
}




