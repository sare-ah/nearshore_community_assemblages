#--------------------------------------------
# Indicator Species Analysis
#
# To Do: Add description
# To Define: difference btw IndVal & IndVal.g
#--------------------------------------------

# Approach (from Caceres et. al. 2012)
#---------
# 1. Cluster data --- look at fuzzy cluster
# 2. Calculate number of sites within each cluster
# 3. Create new community data matrix with species that occur in > 40% of target group (or a lower number?)
# 4. Indicator Analysis with species combinations
# 5. Estimate 95% bootstrap CI
# 6. Determine relationship between min positive predictor value for valid indicators & resulting coverage of the target site
# 7. Filter indicators by A=0.6 and B=0.25

# start fresh
rm(list=ls())

# Load packages
#--------------

library(tidyverse)
#library(labdsv) # indval() indicator species analysis
library(indicspecies) # multipatt() multi-level pattern analysis
library(rstudioapi)


# Get path for this script
#-------------------------
# Set working directory to one above script directory
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
getwd()

# Source functions
#-----------------
#source('CommunityAssemblages_functions.R')

# Set up new directory for all results
#-------------------------------------
date <- format(Sys.Date(), "%b_%d")
region <- "All"
outdir <- paste0(date,"_",region)
#dsn.dir <- "SHP"

setwd( "../Results") 
if (dir.exists(outdir)){
  setwd(outdir)
} else{
  dir.create(path = outdir)
  setwd(outdir)
  #dir.create(dsn.dir)
}
getwd()

# Inputs
#-------

# pValCutoff <- 0.05
# indvalCutoff <- 0.15


# Read in RDS files --- #1 & #2
#------------------

# Species by region
speciesFullCl <- readRDS("../../RDS/speciesFullCl.RDS")
# Determine the number of elements in list to use when building output lists
n <- length(speciesFullCl)

# Species list
species <- readRDS("../../RDS/species.RDS")

# Cluster frequency
cluster.freq <- readRDS("cluster.freq.RDS")

# Determine which species occur in atleast 40% of target clusters --- #3
#----------------------------------------------------------------
# # Grab just one list element
sppAll <- speciesFullCl$ALL
cl.freq <- cluster.freq$ALL
# 
# # Create a long pivot table of all species present in each trans
sppAll_long <- pivot_longer(sppAll, cols = 1:109, names_to = "Spp", values_to = "Presence")
sppAll_long <- dplyr::filter( sppAll_long, Presence==1)

# Calculate number of sites for each cluster/species combination - reduce size of species combo matrix
df <- as.data.frame( sppAll_long %>%
  group_by(cl, Spp) %>%
  summarise(nSites = length(TransDepth)) )

# Threshold of sites for each cluster
thrshld <- data.frame( cl=as.integer(cl.freq$cl), thrshld = round(0.20 * cl.freq$Freq) )

# Determine number of species that meet frequency threshold
newSpp <- left_join(df, thrshld, by="cl")
newSpp$IN <- (newSpp$nSites - newSpp$thrshld)
newSpp <- filter(newSpp, IN>=0)
targetSpp <- unique(newSpp$Spp)
targetSpp

# Subset original site X species matrix
# Complicated...to do...
# Other ideas, build a list with df for each cluster?
# Remove species that is not found in at least one cluster
# Create new species list
# Subset original matrix by new list

row.names(sppAll) <- sppAll$TransDepth
#sppObs <- sppObs[,1:109]
clusters <- sppAll$cl

targetSpp <- species[species%in% names(sppAll)]
sppObs <- sppAll[,targetSpp]
sppObs[is.na(sppObs)] <- 0
head(sppObs, 3)

# Multipatt() multi-level pattern analysis --- #4
#-----------------------------------------

indval = multipatt(x = sppObs, cluster = clusters, max.order = 4, control = how(nperm=999))
summary(indval)

# Test the association btw species & each group of sites, regardless of whether the association value 
# was the highest or not. For example, test whether the frequency of the species in each site group 
# is higher or lower than random
prefsign = signassoc(sppObs, cluster=clusters, alternative = "two.sided", control = how(nperm=199))
head(prefsign)


# Quantity coverage of the site group
# The proportion of sites of a given site group where one or another indicator is found
indvalori = multipatt(sppObs, clusters, duleg = TRUE, control = how(nperm=999))
summary(indvalori)
# Input community data, object of class multipatt
coverage(sppObs,indvalori)
coverage(sppObs, indvalori, At = 0.8, alpha = 0.05)

par(mfrow = c(1,1))
# Plot how coverage changes with 'A' threshold
plotcoverage(x=sppObs, y=indvalori, group="1", lty=1)
plotcoverage(x=sppObs, y=indvalori, group="2", lty=2, col="blue", add=TRUE)
plotcoverage(x=sppObs, y=indvalori, group="3", lty=3, col="red", add=TRUE)
plotcoverage(x=sppObs, y=indvalori, group="4", lty=3, col="green", add=TRUE)
#plotcoverage(x=sppObs, y=indvalori, group="5", lty=3, col="purple", add=TRUE)
legend(x = 0.01, y=25,legend=c("group 1","group 2","group 3","group 4"),
       lty=c(1,2,3), col=c("black","blue","red","green"), bty="n")

# Species combinations as indicators of site groups
#--------------------------------------------------
# Build matrix with all possible species combinations
# My computer can't handle max.order=4
sppComb <- combinespecies(sppObs, max.order = 3)$XC
dim(sppComb)
# Re-run mulitpatt with species combinations
indvalspcomb = multipatt(sppComb, clusters, duleg = TRUE, control = how(nperm=999))
summary(indvalspcomb, indvalcomp = TRUE)

# Determine indicators for group 2 
## Determine sensitivity of individual species
B=strassoc(sppObs, cluster=clusters,func="B")
## Select species with more than 20% of sensitivity for the first group
sel=which(B[,1]>0.2)
sc= indicators(X=sppObs[,sel], cluster=clusters, group=1, max.order = 2, 
               verbose=TRUE, At=0.7, Bt=0.4)
print(sc) 
## Plots positive predictive power and sensitivity against the order of combinations
plot(sc, type="A")
plot(sc, type="B")
## Run indicator analysis with species combinations for the first group,
## but forcing 'Orysp' to be in all combinations
sc2= indicators(X=sppObs[,sel], cluster=clusters, group=1, verbose=TRUE, At=0.5, Bt=0.2, enableFixed=TRUE)
print(sc2) 
plot(sc2, type="A")
plot(sc2, type="B")


sc= indicators(X=sppComb, cluster=clusters, group=2, max.order = 2, 
               verbose=TRUE, At=0.5, Bt=0.3)
print(sc, sqrtIVt = 0.02) # throws row.names error
# Do combinations improve coverage
plotcoverage(sc)
plotcoverage(sc, max.order=1, add=TRUE, lty=2, col="red")
legend(x=0.1, y=20, legend=c("Species combinations","Species singletons"),
       lty=c(1,2), col=c("black","red"), bty="n")

# prune indicators
# output does not match example in tutorial
sc2=pruneindicators(sc, At=0.5, Bt=0.2, verbose=TRUE)
print(sc2)

# predict indicators
pcv <- predict(sc2, sppObs, cv=TRUE)
pcv1 <- predict(sc, cv=TRUE)

# Compared predicted probabilities for each site
# ???
data.frame(Group1 = as.numeric(speciesFullCl$ALL$cl==1), Prob = pcv, Prob_CV = pcv)


# *** TO Do *** 
# Look at RDS by region file, it appears to contain only one region 4 times!
# # # # 5. Join with species look-up table
# # # spLookup<-read.csv("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/LookupTbls/SpeciesLookUpTbl.csv")
# # # #names(spLookup)[1]<-"Species_Code"
# # # topSpdf<-merge(topSpdf, spLookup, by.x="species", by.y="Sp_cde")
# # # # topSpdf<-topSpdf[,c(ncol(topSpdf), 2:(ncol(topSpdf)-1))]
# # # # topSpdf<-topSpdf[order(topSpdf$maxcl, -topSpdf$indvalInMaxcl),]
# # # head(topSpdf,3)
# # # topSpdf <- topSpdf[c(1,11,12,2:10)]
# # 
# # # Save output
# # saveRDS(indSpp, "IndicatorSpecies.RDS")
# # 
