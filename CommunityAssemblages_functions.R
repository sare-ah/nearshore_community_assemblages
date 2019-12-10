#######################################################################################################################
# Adapted from Cluster analysis of Benthic Habitat Mapping Dive Survey Sites
# TO DO: EDIT THIS
# Author:     Katie Gale
#             Katie.Gale@dfo-mpo.gc.ca
#             250-363-6411
#
# Date:       Sept 26, 2018
######################################################################################################################

# Apply threshold for "rare" species & sites that have only a few observations, 
# these can skew the cluster analysis
rare <- function(speciesObs, minSp, minSite ){
  # Calculate number of grids each species occurs in
  specsums <- data.frame(species=names(colSums(speciesObs)), count=colSums(speciesObs))
  # Remove species recorded less than X times
  rare <- subset(specsums, count<=(minSp*nrow(speciesObs)))
  remRare <- speciesObs[,!(colnames(speciesObs) %in% rare$species)]
  cat(c("Threshold for the number of observations needed to be included is", minSp*nrow(speciesObs),"\n"))
  cat(c( nrow(rare), " rare species were removed!","\n"))
  # Remove sites that only have one or two species *throws off cluster analysis)
  sitesums <- data.frame(site=rownames(remRare), count=rowSums(remRare))
  # Remove sites with only a couple of species
  barren <- subset(sitesums, count<=minSite)
  cat(c( nrow(barren), " barren sites were removed!","\n"))
  df <- remRare[!(rownames(remRare) %in% barren$site),]
  print(dim(df))
  rare.List <- list(df, rare, barren)
  return(rare.List)
}


# Create a shapefile
makeShp <- function(data, dsn, layer ){
  coordinates(data) <- ~ x + y
  geoCRS <-  CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
  proj4string(data) <- geoCRS
  writeOGR(data, dsn = dsn, layer=layer, driver="ESRI Shapefile", overwrite=T)
  return(data)
}

# Build my own colour scheme: h=hue, l=lightness, s=saturation
myColors <- function(n, h = c(120,400), l = c(.40,.70), s = c(.8,1), alpha = 1){
  return (alpha(hex(HLS(seq(h[1],h[2],length.out = n), seq(l[1],l[2],length.out = n), seq(s[1],s[2],length.out=n))), alpha))
}

# Check which packages are being used
#library(NCmisc)
#list.functions.in.file("C:/Users/daviessa/Documents/R/PROJECTS_MY/CommunityAssemblages/Scripts/CommunityAssemblages_functions.R")
