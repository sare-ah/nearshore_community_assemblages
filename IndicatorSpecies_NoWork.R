# Indicator Species Analysis

# start fresh
rm(list=ls())

# Load packages
#--------------

library(tidyverse)
library(labdsv) # indval() indicator species analysis
library(indicspecies) # multipatt() multi-level pattern analysis

# Inputs
#-------

pValCutoff <- 0.05
indvalCutoff <- 0.15


# Read in RDS files
#------------------

# Species by region
speciesFullCl <- readRDS("speciesFullCl.RDS")
# Species list
species <- readRDS("../../RDS/species.RDS")


# Indval() analysis
#------------------

# Empty lists for output
indvalList <- vector("list", 4)

# Set-up and run indval() analysis
for (i in 1:length(speciesFullCl)){
  # Match species in species list as speciesNew (observations change in different regions)
  observed <- species[species%in% names(speciesFullCl[[i]])]
  sppObs <- speciesFullCl[[i]][,observed]
  # Ensure there are no NA's - will create error in indval()
  sppObs[is.na(sppObs)] <- 0
  # Build cluster vector
  clusters <- speciesFullCl[[i]]$cl
  # Run indval analysis 
  indvalList[[i]] <- indval(sppObs, clusters)
  names(indvalList)[[i]] <- names(speciesFullCl)[i]
}

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
indval.summary <- vector("list", 4)
for (i in 1:length(indvalList)){
  indval.summary[[i]] <- indicatorValues(indvalList[[i]])
  names(indval.summary)[[i]] <- names(indvalList)[i]
  print(head(indval.summary[[i]],3))
}


# Filter by pval cutoff
#----------------------

# List of species with significant p values
indval.pvals <- vector("list", 4)

# Filter results by p value
for (i in 1:length(indval.summary)){
  if(!is.na(pValCutoff)){
    indval.pvals[[i]] <- indval.summary[[i]][indval.summary[[i]]$pval <= pValCutoff,] 
  } else { indval.pvals[[i]] <- indval.summary[[i]]}
  names(indval.pvals)[[i]] <- names(indval.summary)[i]
}
summary(indval.pvals)


# Filter by cluster and indval cutoff
#------------------------------------

final <- vector("list", 4)

for (i in 1:length(indval.pvals)){
  df <- indval.pvals[[i]]
  nmes <- paste0("indval_",seq(1:max(df$maxcl)))
  df$indvalInMaxcl <- apply(df[,nmes], 1, max)
  df <- df %>% 
    group_by(maxcl) %>% 
    dplyr::filter_at(vars(starts_with("indval_")), any_vars( . >= indvalCutoff ) ) %>% 
    dplyr::select(-pval) %>%
    arrange(maxcl) %>% 
    dplyr::select(species, maxcl, indvalInMaxcl, everything())
  final[[i]] <- as.data.frame(df)
  names(final)[[i]] <- names(indval.pvals)[i]
}
final

# Multipatt() multi-level pattern analysis
#-----------------------------------------

# Empty lists for output
mlpList <- vector("list", 4)

# Set-up and run indval() analysis
for (i in 1:length(speciesFullCl)){
  # Match species in species list as speciesNew (observations change in different regions)
  observed <- species[species%in% names(speciesFullCl[[i]])]
  sppObs <- speciesFullCl[[i]][,observed]
  # Ensure there are no NA's - will create error in indval()
  sppObs[is.na(sppObs)] <- 0
  # Build cluster vector
  clusters <- speciesFullCl[[i]]$cl
  # Run multipatt analysis
  mlpList[[i]] <- multipatt(sppObs, clusters, control = how(nperm=999))
  names(mlpList)[[i]] <- names(speciesFullCl)[i]
}

# Summary of output
#-----------------
summary(mlpList$HG)
# To determine which species had high IndVal in all groups &
# therefore cannot statistically test their association
#*** To Do: make df of sign with p.value == NA 
mlpList$HG$sign

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
