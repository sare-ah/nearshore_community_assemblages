#######################################################################################################################
# Cluster analysis of Benthic Habitat Mapping Dive Survey Sites
#
#######################################################################################################################

# Start fresh
rm(list=ls())

# Load packages
#--------------

library(simba)
library(dendroextras)
library(colorspace)
library(scales)
library(mapview)
library(rgdal)
library(sp)
library(sf)
library(dendextend)
library(rstudioapi)
library(tidyverse)

# Get path for this script
#-------------------------
# Set working directory to one above script directory
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
getwd()

# Source functions
#-----------------
source('CommunityAssemblages_functions.R')

# Set up new directory for all results
#-------------------------------------
date <- format(Sys.Date(), "%b_%d")
region <- "All"
outdir <- paste0(date,"_",region)
dsn.dir <- "SHP"

setwd( "../Results") 
if (dir.exists(outdir)){
  setwd(outdir)
} else{
  dir.create(path = outdir)
  setwd(outdir)
  dir.create(dsn.dir)
}
getwd()

# Inputs
#-------
spThreshold <- 0.03 # Proportion of sites that a species must be found in
siteThreshold <- 3 # Minimum number of species / site
distance <- "simpson"
clusterMethod <- "ward.D2" 

# Read in species by regions
#---------------------------
sppByRegion <- readRDS("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/RDS/sppByRegion.rds")
sppDF <- readRDS("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/RDS/sppByRegion_Dataframe.rds")
allCl <- readRDS("C:/Users/daviessa/Documents/R/PROJECTS_MY/CommunityAssemblages/RDS/clusterAllSites.RDS")

# Read in spatialized points
#---------------------------
dsn <- "C:/Users/daviessa/Documents/R/PROJECTS_OTHERS/SpatializeDiveTransects/ByDepthCat"
sp <- readOGR(dsn=dsn, layer="All_SpatPts_DepthCat", stringsAsFactors = F)

# Create spatialized points link
#-------------------------------
sp.pts <- cbind(sp@data,sp@coords)

# Cluster Analysis
#-----------------

# Add all sites
sppDF$Region <- NULL
sppByRegion$ALL <- sppDF

# --- Remove rare species and barren sites from each region --- #
# Create empty lists needed for outputs --- maybe we only need the listOfLists???
listOfLists <- vector("list", 5)
prep4Cluster <- vector("list", 5)
rareSpp <- vector("list", 5)
barrenSites <- vector("list", 5)

# Calculate summaries
for (i in 1:length(sppByRegion)){
  # Get region name
  names(listOfLists)[[i]] <- names(sppByRegion)[i]
  spp <- as.data.frame(sppByRegion[[i]])
  # Build list of all speciesXsites, rare species, barren sites
  listOfLists[[i]]$sppSiteList <- rare( speciesObs=spp, minSp=spThreshold, minSite=siteThreshold )
  # Build new list of species P/A for cluster analysis
  names(prep4Cluster)[[i]] <- names(sppByRegion)[i]
  prep4Cluster[[i]] <- listOfLists[[i]]$sppSiteList[1]
  # Build list of rare species 
  names(rareSpp)[[i]] <- names(sppByRegion)[i]
  rareSpp[[i]] <- listOfLists[[i]]$sppSiteList[2]
  # Build list of barren sites 
  names(barrenSites)[[i]] <- names(sppByRegion)[i]
  barrenSites[[i]] <- listOfLists[[i]]$sppSiteList[3]
}
saveRDS(listOfLists, "ListOfLists.RDS")


# --- Examine rare species and plot barren sites from each region --- #

# Rare species
# Flatten list into a dataframe to assess differences between regions --- warnings bc df's contain factors, not characters
suppressWarnings(rare <- map_dfr(rareSpp,~map_dfr(.x,identity,.id="species"),identity,.id="region"))
head(rare)

rare.wide <- tidyr::pivot_wider(rare, names_from = "region",values_from = "count") 
path <- paste0(getwd(),"/rareSpecies.csv")
write_csv(rare.wide, path=path) 

# Plot barren sites
list <- unlist(barrenSites, recursive = FALSE)
barren <- do.call("rbind", list)
head(barren)
colnames(barren) <- c("TransDepth","spCnt")
barren$TransDepth <- as.character(barren$TransDepth)
barren <- left_join( barren, sp.pts, by="TransDepth" )
barren <- drop_na(barren)
# Make sf and plot
barren_sf <- st_as_sf(barren, coords = c("coords.x1","coords.x2"), crs = 3005) # BC Albers 
# As "Spatial" - if you need to save as a shp
barren_sp <- as(barren_sf, "Spatial")
mapview(barren_sf, zcol = "spCnt", legend = TRUE) 
sp.filename <- "barren_sites"
dsn.brrn <- paste0(getwd(), "/SHP")
writeOGR(barren_sp, dsn=dsn.brrn, sp.filename, driver="ESRI Shapefile",overwrite_layer = T)

rm(barrenSites,rareSpp,listOfLists,rare.wide,dsn.brrn,barren,barren_sp)

# --- Run cluster analysis and build dendrogram --- #
# Create empty list to store results
myCluster <- vector("list", 5)

# Run cluster analysis
for (i in 1:length(prep4Cluster)){
  # Get region name
  names(myCluster)[[i]] <- names(prep4Cluster)[i]
  print(names(myCluster)[[i]])
  spp <- as.data.frame(prep4Cluster[[i]])
  # Calculate B-sim (simpson) distance on the site x species matrix
  myCluster[[i]]$dist <- sim( spp,  method=distance )
  # Create dendrogram, using average grouping
  myCluster[[i]]$benthtree <- hclust( myCluster[[i]]$dist, method=clusterMethod ) 
  #plot( myCluster[[i]]$benthtree, hang=-1 )
  # Calculate cophenetic correlation value
  # how well does the dendrogram correlate with the original data?
  myCluster[[i]]$cophenetic <- cor( myCluster[[i]]$dist, cophenetic(myCluster[[i]]$benthtree) )
  print(myCluster[[i]]$cophenetic)
}

#rm(prep4Cluster)

# Plot and save dendrograms as png and as separate objects to play with cutoffs
for (i in 1:length(myCluster)){
  print(names(myCluster)[[i]])
  region <- names(myCluster)[[i]]
  file_name = paste("dendrogram_", region, ".png", sep="")
  title_name = paste(region,"Cluster Dendrogram", sep=" ")
  png(file_name)
  plot(myCluster[[i]]$benthtree, hang=-1, main=region)
  dev.off()
}

# # # --- Manually choose heights to cut trees --- #
# hg <- myCluster$HG$benthtree # seth = 4
# ncc <- myCluster$NCC$benthtree # seth = 5
# qcs <- myCluster$QCS$benthtree # seth = 2
# sog <- myCluster$SoG$benthtree # seth = 2
# 
# benthtree <- sog
# plot(benthtree, hang=-1)
# 
# # Choose cutoff, play with h until the visual clusters are kept together
# # Set h (height) for the dendrogram and cut the tree
# #seth <- readline(prompt = "Set height for the dendrogram and cut the tree: ")
# seth <- 2
# rect.hclust(benthtree, h=seth, border="red") # Cutoff based on visual inspection of the tree

# Vector to entry heights from plots
hts <- c(4,5,2,2,8) # corresponds to HG,NCC,QCS,SoG

# Add height to large list
for (i in 1:length(myCluster)){
  myCluster[[i]]$height <- hts[i] 
}

# --- Slice trees --- #

# Final h choice to slice tree 
for (i in 1:length(myCluster)){
  # Get region name
  print(names(myCluster)[[i]])
  myCluster[[i]]$sliceTree <- dendroextras::slice(x=myCluster[[i]]$benthtree, h=myCluster[[i]]$height)
}

# Build a new list of sliced trees
cl.list <- vector("list",5)
for (i in 1:length(myCluster)){
  names(cl.list)[i] <- names(myCluster)[i]
  cl.list[[i]] <- myCluster[[i]]$sliceTree
}

# --- Determine number of clusters to capture 90% of samples --- #
# Get table of cluster memberships - number of sites in each cluster
# lapply(list, function)

# Convert list element to a dataframe
make_df <- function(x){
  as.data.frame(table(x))
}
colorcount <- lapply(cl.list, make_df )
c.names <- c("cl","Freq")
colorcount <- lapply(colorcount, setNames, c.names)
colorcount

# Order clusters by frequency order 
order.cl <- function(x){
  order <- data.frame( order=seq( 1:nrow(x) ),x[order(-x$Freq),] ) 
  order$cumsum <- cumsum( order$Freq )
  order$cumpercent <- round( order$cumsum/max(order$cumsum), 2 )
  order$percent <- round( order$Freq/sum(order$Freq),2 )
  plot( order$order, order$Freq, ylab="n Samples" )
  plot( order$order, order$cumpercent )
  abline( h=0.9, col="red" )
  return(order)
}
par(mfrow = c(2,2))
cluster.frq <- lapply(colorcount, order.cl)
cluster.frq
par(mfrow = c(1,1))

# Number of clusters that capture 90% of samples
nTopClusters <- list(5,4,5,4,4)
buildseq <- function(x){
  seq(1,(x),by=1)
}
myTopClusters <- lapply( nTopClusters, buildseq )
myTopClusters

# Assign colours to clusters
for (i in 1:length(cluster.frq)){
  clusters <- cluster.frq[[i]]$cl[cluster.frq[[i]]$order %in% myTopClusters[[i]] ]
  colorscheme <- myColors(max(myTopClusters[[i]]) )
  df <- colorcount[[i]]
  df <- dplyr::arrange(df, desc(Freq)) 
  df <- dplyr::mutate(df, assigned=ifelse(cl%in%clusters,colorscheme,"grey"))
  df <- dplyr::arrange(df, cl)
  myCluster[[i]]$colours <- df
  df <- NULL
}

#myList <- myCluster

# --- Create color-coded dendrogram....slow! --- #
#colourTree <- vector("list",4)
par(mfrow = c(2,2))
# Create coloured trees
for (i in 1:length(myCluster)){
  main <- names(myCluster)[[i]]
  myCluster[[i]]$tree <- colour_branches( myCluster[[i]]$benthtree, h=myCluster[[i]]$height,col=as.character(myCluster[[i]]$colours$assigned) )
  plot(myCluster[[i]]$tree, main=main)
  myCluster[[i]]$legend <- dplyr::arrange(myCluster[[i]]$colours,  desc(myCluster[[i]]$colours$Freq))
}

# Save coloured trees
for (i in 1:length(myCluster)){
  region <- names(myCluster[i])
  print(region)
  ht <- myCluster[[i]]$height
  file <- paste0("Cluster",region,"_h",ht,".pdf")
  pdf(file)  
  plot(myCluster[[i]]$tree, xlab=NA, ylab="Simpson Distance, Ward.D2 Clustering", 
     main=paste(region," Community assemblages"))
  legend("topright", legend = myCluster[[i]]$legend$cl,fill = as.character(myCluster[[i]]$legend$assigned), bty = "n", ncol=7, horiz = FALSE)
  dev.off()
}

# Save cluster objects
saveRDS(myCluster, "myCluster.RDS")

# --- Assign a cluster group to each site --- # 

# Join cl with Prep4Cluster list 
speciesFullCl <- vector("list", 5)

for (i in 1:length(prep4Cluster)){
  # Organize cluster data into a df
  clTbl <- data.frame("TransDepth"=names(cl.list[[i]]),"cl"=(cl.list[[i]]) )
  clTbl$TransDepth <- as.character(clTbl$TransDepth)
  # Organize siteXspp data into a df and create new variable for join
  spp <- as.data.frame(prep4Cluster[[i]])
  spp$TransDepth <- as.character(row.names(spp))
  # Join clusters to siteXspp data
  spp <- full_join(spp, clTbl, by="TransDepth")
  # Save in new list
  names(speciesFullCl)[[i]] <- names(prep4Cluster)[[i]]
  speciesFullCl[[i]] <- spp
}

# Save work
saveRDS( speciesFullCl, "speciesFullCl.RDS" )

# To Do: join with spat
# layer <- paste0("Cluster_treeHt",seth)
# matFullCl <- makeShp( data=matFullCl, dsn=dsn, layer=layer )
# mapview( matFullCl, zcol = "cl" )


