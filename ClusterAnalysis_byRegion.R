#######################################################################################################################
# Cluster analysis of Benthic Habitat Mapping Dive Survey Sites
#
# To do: Change myCluster to $Element$Region
#        Change loops to lapply statements 
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
library(dendextend)
library(rstudioapi)
library(tidyverse)

# Inputs
#-------

date <- format(Sys.Date(), "%b_%d")
region <- "All"
outdir <- paste0(date,"_",region)
dsn <- "SHP"

spThreshold <- 0.03 # Proportion of sites that a species must be found in
siteThreshold <- 3 # Minimum number of species / site
distance <- "simpson"
clusterMethod <- "ward.D2" 
# seth <- 4 # number for cutoff height
# nTopClusters <- 5 #integer for number of clusters


# Site by Species csv file
# myFile <- choose.files(default = "T:/Benthic_Habitat_Mapping/Data",
#                        caption = "Select site X species, with enviromental variables") 
# myFile <- "T:/Benthic_Habitat_Mapping/Data/Species by Site Matrices/qcsbyDepthCat_AllSpp.csv"

# Read in species by regions
#---------------------------
sppByRegion <- readRDS("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/RDS/sppByRegion.rds")
sppDF <- readRDS("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/RDS/sppByRegion_Dataframe.rds")

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
setwd( "../Results") 
if (dir.exists(outdir)){
  setwd(outdir)
} else{
  dir.create(path = outdir)
  setwd(outdir)
  dir.create(dsn)
}
getwd()


# Import data
# -----------
# # Format is 1 row per sample (depth bin), with columns for species and environmental variables 
# matFull <- read.csv(myFile, header=T, sep=",", stringsAsFactors = F)
# head(matFull,3)

# Look-up table to match species codes to latin names
#luTbl <- read.csv( "T:/Benthic_Habitat_Mapping/Data/Look-upTbls/SpeciesLookUpTbl.csv", header=T, sep=",", stringsAsFactors=F )

# Cluster Analysis
#-----------------

# --- Remove rare species and barren sites from each region --- #
# Create empty lists needed for outputs --- maybe we only need the listOfLists???
listOfLists <- vector("list", 4)
prep4Cluster <- vector("list", 4)
rareSpp <- vector("list", 4)
barrenSites <- vector("list", 4)

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

### TO DO ###
# Outputs: Table of rare species by region (spp | NCC | HG | QCS | SoG)
# Outputs: Figure of barren sites
# Save as ? RDS, shp, csv

# # Keep all species
# forCl <- as.data.frame( sppSiteList[1] ) 
# 
# # Remove rare species
# remSp <- as.data.frame( sppSiteList[2] ) 
# as.character( remSp$species )
# write_csv( remSp, paste0(region,"_SpeciesRemoved.csv" ))
# 
# # Remove barren sites
# remSites <- as.data.frame( sppSiteList[3] ) 
# colnames(remSites) <- c( "TransDepth","SpCnt" )
# remSites$TransDepth <- as.character(remSites$TransDepth)
# remSites <- left_join( remSites, matFull, by="TransDepth" )
# remSites[1,1:6] # confirm that position columns are named x,y
# layer <- paste0(region,"_SitesRemoved" )
# remSites <- makeShp( data=remSites, dsn=dsn, layer=layer )
# mapview(remSites, zcol="SpCnt")
# 
# saveRDS( sppSiteList, "sppSiteList.rds" ) 

rm(barrenSites,rareSpp,sppByRegion,sppDF,listOfLists)

# --- Run cluster analysis and build dendrogram --- #
# Create empty list to store results
myCluster <- vector("list", 4)

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

rm(prep4Cluster)

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
# hts <- c(4,5,2,2)

# Add height to large list
myCluster$HG$height <- 4 # Height on y-axis
myCluster$NCC$height <- 5
myCluster$QCS$height <- 2
myCluster$SoG$height <- 2

# Final h choice to slice tree 
for (i in 1:length(myCluster)){
  # Get region name
  print(names(myCluster)[[i]])
  myCluster[[i]]$sliceTree <- dendroextras::slice(x=myCluster[[i]]$benthtree, h=myCluster[[i]]$height)
}

# Build list of sliced trees
cl.list=list(list(c(myCluster$HG$sliceTree)),list(c(myCluster$NCC$sliceTree)),list(c(myCluster$QCS$sliceTree)),list(c(myCluster$SoG$sliceTree)))
names(cl.list) <- names(myCluster)

# --- Determine number of clusters to capture 90% of samples --- #
# Get table of cluster memberships - number of sites in each cluster
# lapply(list, function)

# Convert list element to a dataframe
make_df <- function(x){
  as.data.frame(table(x))
}
colorcount <- lapply(cl.list, make_df )
#names(colorcount) <- names(myCluster)
c.names <- c("cl","Freq")
colorcount <- lapply(colorcount, setNames, c.names)
colorcount

# Order clusters by frequency order 
order.cl <- function(x){
  order <- data.frame( order=seq( 1:nrow(x) ),x[order(-x$Freq),] ) 
  order$cumsum <- cumsum( order$Freq )
  order$cumpercent <- round( order$cumsum/max(order$cumsum), 2 )
  order$percent <- round( order$Freq/sum(order$Freq),2 )
  #plot( order$order, order$Freq, ylab="n Samples" )
  # plot( order$order, order$cumpercent )
  # abline( h=0.9, col="red" ) 
  return(order)
}
par(mfrow = c(2, 2))
cluster.frq <- lapply(colorcount, order.cl)
cluster.frq

# Number of clusters that capture 90% of samples
#nTopClusters <- c(5,4,5,4)

# Could not get %in% to work for elements within a list, so here goes some ugly code!
hg <- cluster.frq$HG
ncc <- cluster.frq$NCC
qcs <- cluster.frq$QCS
sog <- cluster.frq$SoG

#HG
hg.cluster <- hg$cl[hg$order %in% c( 1:5 )] 
colorscheme <- myColors(5)
colorHG <- colorcount$HG
for (i in c(1:5)){
  colorHG$assigned[colorHG$cl==hg.cluster[i]] <- colorscheme[i]
}
table(colorHG$assigned)
colorHG$assigned[is.na(colorHG$assigned)] <- "grey"
colorHG
#NCC
ncc.cluster <- ncc$cl[ncc$order %in% c( 1:4 )] 
colorscheme <- myColors(4)
colorNCC <- colorcount$NCC
for (i in c(1:4)){
  colorNCC$assigned[colorNCC$cl==ncc.cluster[i]] <- colorscheme[i]
}
table(colorNCC$assigned)
colorNCC$assigned[is.na(colorNCC$assigned)] <- "grey"
colorNCC
#QCS
qcs.cluster <- qcs$cl[qcs$order %in% c( 1:5 )] 
colorscheme <- myColors(5)
colorQCS <- colorcount$QCS
for (i in c(1:5)){
  colorQCS$assigned[colorQCS$cl==qcs.cluster[i]] <- colorscheme[i]
}
table(colorQCS$assigned)
colorQCS$assigned[is.na(colorQCS$assigned)] <- "grey"
colorQCS
#SoG
sog.cluster <- sog$cl[sog$order %in% c( 1:4 )] 
colorscheme <- myColors(4)
colorSoG <- colorcount$SoG
for (i in c(1:4)){
  colorSoG$assigned[colorSoG$cl==sog.cluster[i]] <- colorscheme[i]
}
table(colorSoG$assigned)
colorSoG$assigned[is.na(colorSoG$assigned)] <- "grey"
colorSoG

# Add colour assignments back to myCluster list
# There must be a better way than this
myCluster$HG$colours <- colorHG
myCluster$NCC$colours <- colorNCC
myCluster$QCS$colours <- colorQCS
myCluster$SoG$colours <- colorSoG
rm(colorHG,colorNCC,colorQCS,colorSoG,sog,qcs,ncc,hg
   )

# --- Create color-coded dendrogram....slow... --- #
colourTree <- vector("list",4)

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


