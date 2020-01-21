# Indicator Species Analysis

# start fresh
rm(list=ls())

# Load packages
#--------------

library(tidyverse)
library(labdsv) # indval() indicator species analysis
library(indicspecies) # multipatt() multi-level pattern analysis
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

# Source functions
#-----------------
source('CommunityAssemblages_functions.R')

# Inputs
#-------

pValCutoff <- 0.05
indvalCutoff <- 0.15


# Read in RDS files
#------------------

# All Site X Species 
sppDF <- readRDS("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/RDS/sppByRegion_Dataframe.rds")
summary(sppDF)
sppDF$Region <- NULL
# Species list
species <- readRDS("../../RDS/species.RDS")

# Prep for Indicator species analysis - Cluster analysis
#-------------------------------------------------------
spThreshold <- 0.03 # Proportion of sites that a species must be found in
siteThreshold <- 3 # Minimum number of species / site
distance <- "simpson"
clusterMethod <- "ward.D2" 

# Build list of all speciesXsites, rare species, barren sites
sppSiteList <- rare( speciesObs=sppDF, minSp=spThreshold, minSite=siteThreshold )
# Build new list of species P/A for cluster analysis
prep4Cluster <- as.data.frame(sppSiteList[1])
# Build list of rare species 
rareSpp <- sppSiteList[2]
# Build list of barren sites 
barrenSites <- sppSiteList[3]

dist <- sim( prep4Cluster,  method=distance )
# Create dendrogram, using average grouping
benthtree <- hclust( dist, method=clusterMethod ) 
plot( benthtree, hang=-1 )
# Calculate cophenetic correlation value
# how well does the dendrogram correlate with the original data?
cophenetic <- cor( dist, cophenetic(benthtree) )
# Cut tree
seth <- 8
rect.hclust(benthtree, h=seth, border="red") # Cutoff based on visual inspection of the tree
# Make final cut of tree
sliceTree <- dendroextras::slice(x=benthtree, h=seth)
# Make df for cluster assignments
clusters <- data.frame("TransDepth"=names(sliceTree),"cl"=(sliceTree) )
clusters$TransDepth <- as.character(clusters$TransDepth)
# Join with prep4Cluster df
prep4Cluster$TransDepth <- as.character(row.names(prep4Cluster))
spp <- full_join(prep4Cluster, clusters, by="TransDepth")
# Save output as RDS
saveRDS(spp, "C:/Users/daviessa/Documents/R/PROJECTS_MY/CommunityAssemblages/RDS/clusterAllSites.RDS")

# Indval() analysis
#------------------

# Match species in species list as speciesNew (observations change in different regions)
observed <- species[species%in% names(spp)]
sppObs <- spp[,observed]
# Ensure there are no NA's - will create error in indval()
sppObs[is.na(sppObs)] <- 0
# Build cluster vector
clusters <- spp$cl
# Run indval analysis 
indvalList <- indval(sppObs, clusters)

# Function to organize indval() output into a summary table by species 
indicatorValues <- function(ind){
  # extract the indicator value for each species
  indv <- data.frame(species=row.names(ind$indval), round(ind$indval,3))
  names(indv) <- c("species",paste0("indval_",gsub("X","",names(indv[,2:ncol(indv)]))))
  # extract the relative abundance of species in classes
  abu <- data.frame(species=row.names(ind$relabu), round((ind$relabu*100),1))
  names(abu) <-c ("species",paste0("freq_",gsub("X","",names(abu[,2:ncol(abu)]))))
  # name the clusters
  cls <- gsub("indval_","",names(indv)[2:length(names(indv))])
  # extract the class each species has maximum indicator value for
  maxcls <- data.frame(species=names(ind$maxcls), maxcl=as.numeric(as.character(cls[ind$maxcls]))) 
  # extract the probability of obtaining as high an indicator values as observed over the specified iterations
  pval <- data.frame(species=names(ind$pval), pval=ind$pval) 
  # build summary table
  indtab <- cbind(maxcls, pval=pval[,-1], indv[,-1], abu[,-1])
  return(indtab)
}
# Run function
indval.summary <- indicatorValues(indvalList)
print(head(indval.summary,3))


# # Filter by pval cutoff
# #----------------------
# for (i in unique(intabSub$maxcl)){
# 
# # Filter results by p value
# for (i in unique(indval.summary$maxcl)){
#   if(!is.na(pValCutoff)){
#     indval.pvals[[i]] <- indval.summary[indval.summary[[i]]$pval <= pValCutoff,] 
#   # } else { indval.pvals[[i]] <- indval.summary[[i]]}
#   # names(indval.pvals)[[i]] <- names(indval.summary)[i]
# }
# summary(indval.pvals)
# 
# 
# # Filter by cluster and indval cutoff
# #------------------------------------
# 
# final <- vector("list", 4)
# 
# for (i in 1:length(indval.pvals)){
#   df <- indval.pvals[[i]]
#   nmes <- paste0("indval_",seq(1:max(df$maxcl)))
#   df$indvalInMaxcl <- apply(df[,nmes], 1, max)
#   df <- df %>% 
#     group_by(maxcl) %>% 
#     dplyr::filter_at(vars(starts_with("indval_")), any_vars( . >= indvalCutoff ) ) %>% 
#     dplyr::select(-pval) %>%
#     arrange(maxcl) %>% 
#     dplyr::select(species, maxcl, indvalInMaxcl, everything())
#   final[[i]] <- as.data.frame(df)
#   names(final)[[i]] <- names(indval.pvals)[i]
# }
# final
# 
# # Multipatt() multi-level pattern analysis
# #-----------------------------------------
# 
# # Empty lists for output
# mlpList <- vector("list", 4)
# 
# # Set-up and run indval() analysis
# for (i in 1:length(speciesFullCl)){
#   # Match species in species list as speciesNew (observations change in different regions)
#   observed <- species[species%in% names(speciesFullCl[[i]])]
#   sppObs <- speciesFullCl[[i]][,observed]
#   # Ensure there are no NA's - will create error in indval()
#   sppObs[is.na(sppObs)] <- 0
#   # Build cluster vector
#   clusters <- speciesFullCl[[i]]$cl
#   # Run multipatt analysis
#   mlpList[[i]] <- multipatt(sppObs, clusters, control = how(nperm=999))
#   names(mlpList)[[i]] <- names(speciesFullCl)[i]
# }
# 
# # Summary of output
# #-----------------
# summary(mlpList$HG)
# # To determine which species had high IndVal in all groups &
# # therefore cannot statistically test their association
# #*** To Do: make df of sign with p.value == NA 
# mlpList$HG$sign
# 
# # Correlation Indices
# #--------------------
# 
# # Correlation indices: Used to determine the ecological preferences of species among a set of alternative sites groups
# # correlation-type indices ('r', 'r.g', 'cos' and 'cos.g')
# 
# # Indicator values: Used to assess the predictive values of species as indicators of the conditions prevailing in site groups
# # IndVal-type indices ('IndVal' and 'IndVal.g')
# 
# # Pearson's phi coefficient of association, fidelity
# # func = "r.g" corrects for some groups having more sites than others
# phi = multipatt(sppObs, clusters, func = "r.g", control = how(nperm=999))
# summary(phi)
# 
# # Negative associations or species that tend to 'avoid' particular environmental conditions
# round(head(phi$str), 3)
# round(head(mlpList$SoG$str),3)
# 
# indvalrest = multipatt(sppObs, clusters, max.order = 2, control = how(nperm=999))
# summary(indvalrest)
# 
# # Test the association btw species & each group of sites, regardless of whether the association value 
# # was the highest or not. For example, test whether the frequency of the species in each site group 
# # is higher or lower than random
# prefsign = signassoc(sppObs, cluster=clusters, alternative = "two.sided", control = how(nperm=199))
# head(prefsign)
# 
# 
# # Quantity coverage of the site group
# # The proportion of sites of a given site group where one or another indicator is found
# indvalori = multipatt(sppObs, clusters, duleg = TRUE, control = how(nperm=999))
# summary(indvalori)
# # Input community data, object of class multipatt
# coverage(sppObs,indvalori)
# coverage(sppObs, indvalori, At = 0.8, alpha = 0.05)
# 
# # Plot how coverage changes with 'A' threshold
# plotcoverage(x=sppObs, y=indvalori, group="1", lty=1)
# plotcoverage(x=sppObs, y=indvalori, group="2", lty=2, col="blue", add=TRUE)
# plotcoverage(x=sppObs, y=indvalori, group="3", lty=3, col="red", add=TRUE)
# plotcoverage(x=sppObs, y=indvalori, group="4", lty=3, col="green", add=TRUE)
# plotcoverage(x=sppObs, y=indvalori, group="5", lty=3, col="purple", add=TRUE)
# legend(x = 0.01, y=25,legend=c("group 1","group 2","group 3","group 4","group 5"),
#        lty=c(1,2,3), col=c("black","blue","red","green","purple"), bty="n")
# 
# # Species combinations as indicators of site groups
# #--------------------------------------------------
# # Build matrix with all possible species combinations
# sppComb <- combinespecies(sppObs, max.order = 2)$XC
# dim(sppComb)
# # Re-run mulitpatt with species combinations
# indvalspcomb = multipatt(sppComb, clusters, duleg = TRUE, control = how(nperm=999))
# summary(indvalspcomb, indvalcomp = TRUE)
# 
# # Determine indicators for group 2 
# ## Determine sensitivity of individual species
# B=strassoc(sppObs, cluster=clusters,func="B")
# ## Select species with more than 20% of sensitivity for the first group
# sel=which(B[,1]>0.2)
# sc= indicators(X=sppObs[,sel], cluster=clusters, group=1, max.order = 2, 
#                verbose=TRUE, At=0.7, Bt=0.4)
# print(sc) 
# ## Plots positive predictive power and sensitivity against the order of combinations
# plot(sc, type="A")
# plot(sc, type="B")
# ## Run indicator analysis with species combinations for the first group,
# ## but forcing 'Orysp' to be in all combinations
# sc2= indicators(X=sppObs[,sel], cluster=clusters, group=1, verbose=TRUE, At=0.5, Bt=0.2, enableFixed=TRUE)
# print(sc2) 
# plot(sc2, type="A")
# plot(sc2, type="B")
# 
# 
# sc= indicators(X=sppComb, cluster=clusters, group=2, max.order = 2, 
#                verbose=TRUE, At=0.7, Bt=0.4)
# print(sc, sqrtIVt = 0.4) # throws row.names error
# # Do combinations improve coverage
# plotcoverage(sc)
# plotcoverage(sc, max.order=1, add=TRUE, lty=2, col="red")
# legend(x=0.1, y=20, legend=c("Species combinations","Species singletons"),
#        lty=c(1,2), col=c("black","red"), bty="n")
# 
# # prune indicators
# # output does not match example in tutorial
# sc2=pruneindicators(sc, At=0.8, Bt=0.2, verbose=TRUE)
# print(sc2)
# 
# # predict indicators
# pcv <- predict(sc2, sppObs, cv=TRUE)
# pcv1 <- predict(sc, cv=TRUE)
# 
# # Compared predicted probabilities for each site
# # ???
# data.frame(Group2 = as.numeric(speciesFullCl$SoG$cl==2), Prob = pcv, Prob_CV = pcv)
# 
# 
# # *** TO Do *** 
# # Look at RDS by region file, it appears to contain only one region 4 times!
# # # # # 5. Join with species look-up table
# # # # spLookup<-read.csv("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/LookupTbls/SpeciesLookUpTbl.csv")
# # # # #names(spLookup)[1]<-"Species_Code"
# # # # topSpdf<-merge(topSpdf, spLookup, by.x="species", by.y="Sp_cde")
# # # # # topSpdf<-topSpdf[,c(ncol(topSpdf), 2:(ncol(topSpdf)-1))]
# # # # # topSpdf<-topSpdf[order(topSpdf$maxcl, -topSpdf$indvalInMaxcl),]
# # # # head(topSpdf,3)
# # # # topSpdf <- topSpdf[c(1,11,12,2:10)]
# # # 
# # # # Save output
# # # saveRDS(indSpp, "IndicatorSpecies.RDS")
# # # 
