#######################################################################################################################
# Cluster analysis of Benthic Habitat Mapping Dive Survey Sites
# 
# Author:     Katie Gale
#             Katie.Gale@dfo-mpo.gc.ca
#             250-363-6411
# Date:       Sept 26, 2018
#
# Edited by:  Sarah Davies
#             Sarah.Davies@dfo-mpo.gc.ca
#             250-756-7124
# Date:       December 6th, 2019
######################################################################################################################

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
region <- "QCS"
outdir <- paste0(date,"_",region)
dsn <- "SHP"
spThreshold <- 0.03 # Proportion of sites that a species must be found in
siteThreshold <- 3 # Minimum number of species / site
distance <- "simpson"
clusterMethod <- "ward.D2" 
seth <- 2.5 # number for cutoff height
nTopClusters <- 5 #integer for number of clusters

# Site by Species csv file
# myFile <- choose.files(default = "T:/Benthic_Habitat_Mapping/Data",
#                        caption = "Select site X species, with enviromental variables") 
myFile <- "T:/Benthic_Habitat_Mapping/Data/Species by Site Matrices/qcsbyDepthCat_AllSpp.csv"


# Get path for this script
#-------------------------
# Set working directory to one above script directory
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# Source functions
#-----------------
source('CommunityAssemblages_functions.R')

# Set up new directory for all results
#-------------------------------------
setwd( "../Results") 
if (file.exists(outdir)){
  setwd(outdir)
} else{
  dir.create(path = outdir)
  setwd(outdir)
  dir.create(dsn)
}
getwd()


# Import data
# -----------
# Format is 1 row per sample (depth bin), with columns for species and environmental variables 
matFull <- read.csv(myFile, header=T, sep=",", stringsAsFactors = F)
head(matFull,3)


# Organize data into species & env
#---------------------------------
# Checks for complete recordset
matFull <- unique(matFull)
matFull <- dplyr::filter(matFull, x!=0)

# Set rownames
rownames(matFull) <- matFull$TransDepth

# Divide columns into species and environmental variabls
species <- names(matFull)[c(17:185)] # specify species columns, if different
species
envVar <- names(matFull)[c(3:16)] # specify enviro columns, if different
envVar
idvars <- (matFull[,1]) # specify unique identifier for each row, if different

# Create matrix just of species columns
matSp <- matFull[,species]
dim(matSp) 

# Build list of all speciesXsites, rare species, barren sites
sppSiteList <- rare( speciesObs=matSp, minSp=spThreshold, minSite=siteThreshold )

# Keep all species
forCl <- as.data.frame( sppSiteList[1] ) 

# Remove rare species
remSp <- as.data.frame( sppSiteList[2] ) 
as.character( remSp$species )
write_csv( remSp, paste0(region,"_SpeciesRemoved.csv" ))

# Remove barren sites
remSites <- as.data.frame( sppSiteList[3] ) 
colnames(remSites) <- c( "TransDepth","SpCnt" )
remSites$TransDepth <- as.character(remSites$TransDepth)
remSites <- left_join( remSites, matFull, by="TransDepth" )
remSites[1,1:6] # confirm that position columns are named x,y
layer <- paste0(region,"_SitesRemoved" )
remSites <- makeShp( data=remSites, dsn=dsn, layer=layer )
mapview(remSites, zcol="SpCnt")

saveRDS( sppSiteList, "sppSiteList.rds" ) 
rm( remSp, remSites, sppSiteList ) 

# List species that will be used in the analysis
spInCl <- names(forCl)
cat("Species included in analysis:",spInCl)

# # Rename "coords_x1","coords_x2" to x, y
# names(matFull)[names(matFull)=="coords.x1"] <- "x"
# names(matFull)[names(matFull)=="coords.x2"] <- "y"


# Cluster analyis
#----------------
# Calculate B-sim (simpson) distance on the site x species matrix
dist <- sim( forCl,  method=distance )

# Create dendrogram, using average grouping
benthtree <- hclust( dist, method=clusterMethod )

# Calculate cophenetic correlation value
# how well does the dendrogram correlate with the original data?
cor( dist, cophenetic(benthtree) )

# Plot clusters
plot( benthtree, hang=-1 )

# Choose cutoff, play with h until the visual clusters are kept together
# Set h (height) for the dendrogram and cut the tree
#seth <- readline(prompt = "Set height for the dendrogram and cut the tree: ")

rect.hclust(benthtree, h=seth, border="red") # Cutoff based on visual inspection of the tree

# Final h choice
cl <- dendroextras::slice(benthtree, h=seth) # Cut tree


# Determine number of clusters to capture 90% of samples  
# ------------------------------------------------------
# Get table of cluster memberships - number of sites in each cluster
colorcount <- as.data.frame(table(cl)) 

# Order clusters by frequency & select the main clusters to be carried foward thru the analysis
# Clusters with very few samples are not very meaningful and are dropped
order <- data.frame( order=seq( 1:nrow(colorcount) ),colorcount[order(-colorcount$Freq),] ) 
order$cumsum <- cumsum( order$Freq )
order$cumpercent <- round( order$cumsum/max(order$cumsum), 2 )
order$percent <- round( order$Freq/sum(order$Freq),2 )
order

# The first few clusters have most of the samples, and then it drops off
plot( order$order, order$Freq, ylab="n Samples" )
# How many clusters captures 90% of the sample points?
plot( order$order, order$cumpercent )
abline( h=0.9, col="red" ) 

# Set cluster choice number
#nTopClusters <- as.integer( readline( prompt = "Set number of cluster that capture 90% of the sample points: " ) )

# Plot coloured dendrogram 
# ------------------------
# Select n clusters to be color-coded
clToInclude <- order$cl[order$order %in% c( 1:nTopClusters )] 
 
# Set colours for each cluster
colorscheme <- myColors(nTopClusters)
colorcount$assigned <- NULL
for (i in c(1:nTopClusters)){
  colorcount$assigned[colorcount$cl==clToInclude[i]] <- colorscheme[i]
}

# Each selected cluster is assigned a color, NA's assigned to grey
table(colorcount$assigned)
colorcount$assigned[is.na(colorcount$assigned)] <- "grey"

# Create color-coded dendrogram....slow...
colortree <- colour_branches( benthtree, h=seth,col=as.character(colorcount$assigned) )
plot(colortree) 

# Create legend for output
colorCntLegend <- colorcount[colorcount$cl %in% clToInclude,] # Select top n clusters
colorCntLegend$cl <- as.character(colorCntLegend$cl)
colorCntLegend <- colorCntLegend[order(-colorCntLegend$Freq),] # Order by frequency
legendcluster <- rbind( colorCntLegend, c("others", NA, "grey") )
legendcluster
saveRDS( legendcluster, "legendcluster.rds" ) 

# Print dendrogram
## small - for publications
#tiff(file=paste0("./",outdir,"/cluster_",nrow(forCl),"x",ncol(forCl),setrem, "spp_Bsim_h",seth,"x.tiff"), res=350, height=8, width=14, units="cm", pointsize=7)
# large - for CSAS
file <- paste0("Cluster",nrow(forCl),"x",ncol(forCl),distance,"_h",seth,".pdf")
pdf(file)  
plot(colortree, xlab=NA, ylab=paste(distance,"Distance,",clusterMethod, "clustering"), 
    main=paste(region, " Community assemblages (",distance," distance) 
               ( h =",seth,") (",nrow(forCl),"sites x ", ncol(forCl),"spp )" ))
  legend("topleft", legend = legendcluster$cl,fill = as.character(legendcluster$assigned), bty = "n", ncol=7)
dev.off()

 
# Assign a cluster group to each site
#------------------------------------
clTable <- data.frame( ID=names(cl), cl=cl )
colnames(clTable) <- c( "TransDepth","cl" )
matFullCl <- merge( matFull, clTable,idvars="TransDepth" )
head(matFullCl, 3)
length( unique(matFullCl$TransDepth) )

# Save work
write_csv( matFullCl, "matFullCl.csv" )
layer <- paste0("Cluster_treeHt",seth)
matFullCl <- makeShp( data=matFullCl, dsn=dsn, layer=layer )
mapview( matFullCl, zcol = "cl" )


 
# library(NCmisc)
# list.functions.in.file("C:/Users/daviessa/Documents/R/PROJECTS_MY/CommunityAssemblages/Scripts/ClusterAnalysis.R")
