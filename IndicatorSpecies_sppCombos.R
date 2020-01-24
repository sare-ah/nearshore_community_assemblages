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
thrshld <- data.frame( cl=as.integer(cl.freq$cl), thrshld = round(0.40 * cl.freq$Freq) )

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

sppObs <- sppAll[,targetSpp]
sppObs[is.na(sppObs)] <- 0
head(sppObs, 3)

# Multipatt() multi-level pattern analysis 
#-----------------------------------------
# 
# indval = multipatt(x = sppObs, cluster = clusters, max.order = 4, control = how(nperm=999))
# summary(indval)
# 
# Test the association btw species & each group of sites, regardless of whether the association value
# was the highest or not. For example, test whether the frequency of the species in each site group
# is higher or lower than random
# ? what is psidak?
# What am I comparing this output too?
prefsign <- signassoc(sppObs, cluster=clusters, alternative = "two.sided", control = how(nperm=199))
head(prefsign)

 
# Quantity coverage of the site group
# The proportion of sites of a given site group where one or another indicator is found
indvalori <- multipatt(sppObs, clusters, duleg = TRUE, control = how(nperm=999))
summary(indvalori)
# Input community data, object of class multipatt
coverage(sppObs,indvalori)
coverage(sppObs, indvalori, At = 0.6, alpha = 0.05)
 
par(mfrow = c(1,1))
# Plot how coverage changes with 'A' threshold
plotcoverage(x=sppObs, y=indvalori, group="1", lty=1)
plotcoverage(x=sppObs, y=indvalori, group="2", lty=2, col="blue", add=TRUE)
plotcoverage(x=sppObs, y=indvalori, group="3", lty=3, col="red", add=TRUE)
plotcoverage(x=sppObs, y=indvalori, group="4", lty=3, col="green", add=TRUE)
#plotcoverage(x=sppObs, y=indvalori, group="5", lty=3, col="purple", add=TRUE)
legend(x = 0.01, y=30,legend=c("group 1","group 2","group 3","group 4"),
       lty=c(1,2,3), col=c("black","blue","red","green"), bty="n")

# Species combinations as indicators of site groups --- #4
#--------------------------------------------------
# Build matrix with all possible species combinations
# My computer can't handle max.order=4
sppComb <- combinespecies(sppObs, max.order = 3)$XC
dim(sppComb)
saveRDS(sppComb, "sppComb.RDS")

# Re-run mulitpatt with species combinations
indvalspcomb = multipatt(sppComb, clusters, duleg = TRUE, control = how(nperm=999))
# List species with a significant association to one combination, including indval components
summary(indvalspcomb, indvalcomp = TRUE)
saveRDS(indvalspcomb, "SpeciesComboMultipatt.RDS")

# Create output with A, B, or indval stat > some threshold to reduce output
#...

# Determine indicators for each group (edit group=_) 
#---------------------------------------------------
# Threshold for positive predictive value
At <- 0.6
# Threshold for sensitivity
Bt <- 0.25

# Determine sensitivity of individual species
# Strength of species site-group associations
# Square root of IndVal index from labdsv pkg
B <- strassoc( sppObs, cluster=clusters, func="B" )

# Loop through clusters and select species with more than 20% of sensitivity for the first group
sc <- vector("list",length(unique(clusters)))

for (i in 1:length(unique(clusters))){
  label <- i
  sel <- which(B[,i]>0.2)
  names(sc)[[i]] <- i
  sc[[i]] <- indicators(X=sppObs[,sel], cluster=clusters, group=i, max.order = 3, 
               verbose=TRUE, At=At, Bt=Bt)
  print(sc[[i]]) 
  plotcoverage(sc[[i]])
  plotcoverage(sc[[i]], max.order=1, add=TRUE, lty=2, col="red")
  legend(x=0.1, y=20, title=label, legend=c("Species combinations","Species singletons"),
         lty=c(1,2), col=c("black","red"), bty="n")
}


## Plots positive predictive power and sensitivity against the order of combinations
plot(sc, type="A")
plot(sc, type="B")

## Run indicator analysis with species combinations for the first group,
## but forcing 'Orysp' to be in all combinations
sc2= indicators(X=sppObs[,sel], cluster=clusters, group=1, verbose=TRUE, At=At, Bt=Bt, enableFixed=TRUE)

sc= indicators(X=sppComb[,sel], cluster=clusters, group=2, max.order = 2, verbose=TRUE, At=0.5, Bt=0.3)
print(sc, sqrtIVt = 0.02) # throws row.names error
summary(sc)

# Determine if combinations improve coverage
plotcoverage(sc)
plotcoverage(sc, max.order=1, add=TRUE, lty=2, col="red")
legend(x=0.1, y=20, legend=c("Species combinations","Species singletons"),
       lty=c(1,2), col=c("black","red"), bty="n")

signassoc(X=sppComb[,sel], cluster=clusters, mode=1, control = how(nperm=999))
## Look for species whose abundance is significantly higher in sites belonging
## to one group as opposed to sites not belonging to it.
signassoc(X=sppComb[,sel], cluster=clusters, mode=0, control = how(nperm=999))


# Prune indicators to determine if coverage is changed with a smaller set of indicators
# output does not match example in tutorial
sc2 <- pruneindicators(sc, At=0.5, Bt=0.2, verbose=TRUE)
print(sc2)

# predict indicators
pcv <- predict(sc2, sppObs, cv=TRUE)
pcv1 <- predict(sc, cv=TRUE)

# Compared predicted probabilities for each site --- ???
data.frame(Group1 = as.numeric(speciesFullCl$ALL$cl==1), Prob = pcv, Prob_CV = pcv)


# 7.

# Create summary output
# Table 1. Cl | Name | Sites | Spp | Ind | Valid | Final | Cover
# Cl = cluster number
# Name = descriptive name 
# Sites = number of sites (cluster.freq.RDS)
# Ind = number of candidate species --- what is candidate? from IndVal()?
# Valid = number of valid indicators
# Final = smallest set of valid indicators with the same coverage as the complete set --- how to calculate??
# Cover = percentage coverage of the final set of valid indicators (output from pruneindicators())

# Table 2. Cl | Name | Indicators | A (95CI) | B (95CI) | sprt(IV)(95CI)

# Cluster | Descript name | No.Sites | No. of Candidate Species | Ind