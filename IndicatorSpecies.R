# Indicator Species Analysis

# start fresh
rm(list=ls())

library(tidyverse)
#library(foreign)
library(vegan)
#library(rgdal)
library(labdsv) # for indicator species analysis

setwd("~/R/PROJECTS_MY/CommunityAssemblages/Results/Cluster4")

pValCutoff <- 0.05
indvalCutoff <- 0.15


# Read in RDS files
speciesFullCl <- readRDS("speciesFullCl.RDS")
#dat <- speciesFullCl$HG

# Read in species list
species <- readRDS("../../RDS/species.RDS")

# Determine the number of elements within the site x species list
nElements <- length(unique(speciesFullCl$ALL$cl))

# Match species in species list as speciesNew (species list changes from region to region)
#speciesNew <- species[species%in% names(dat)]
speciesTrim <- vector("list", nElements)
indvalList <- vector("list", nElements)
for (i in 1:length(speciesFullCl)){
  print(dim(speciesFullCl[[i]]))
  # Select species of interest in each region
  speciesNew <- species[species%in% names(speciesFullCl[[i]])]
  speciesTrim[[i]] <- speciesFullCl[[i]][,speciesNew]
  speciesTrim[[i]][is.na(speciesTrim[[i]])] <- 0
  names(speciesTrim)[[i]] <- names(speciesFullCl)[i]
  print(dim(speciesTrim[[i]]))
  # Run indval analysis
  indvalList[[i]] <- indval(speciesTrim[[i]], speciesFullCl[[i]]$cl)
  names(indvalList)[[i]] <- names(speciesFullCl)[i]
}


# Organize indval output into a dataframe

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
indval.summary <- vector("list", nElements)
for (i in 1:length(indvalList)){
  indval.summary[[i]] <- indicatorValues(indvalList[[i]])
  names(indval.summary)[[i]] <- names(indvalList)[i]
  print(head(indval.summary[[i]],3))
}


# Select only species that have a pval less than the predetermined cutoff?
# step1 <- vector("list", nElements)
# step2 <- vector("list", nElements)
signifSpp <- vector("list", nElements)

pValCutOffs <- function(x){
  if(!is.na(pValCutoff)){
    signifSpp[[i]] <- indval.summary[[i]][indval.summary[[i]]$pval <= pValCutoff,] # change name of new df
  } else { signifSpp[[i]] <- indval.summary[[i]]}
}

# Apply function
indval.pvals <- lapply(indval.summary, pValCutOffs)


#listElement <- vector("list", nElements)

# Function to select indicator species
indicatorSpecies <- function(df, indvalCutoff){
  # New empty list
  topSp <- list()
  #topSpdf <- list()
  # # df to hold list element i
  #df <- listElement[[i]]
  df$species <- as.character( df$species )
  clusters <- unique(df$maxcl)
  for (i in 1:length(clusters)){
    print(i)
    # Loop thru each cluster...
    df.cl.X <- df[df$maxcl==(i), ]
    #df.cl.X <- dplyr::filter(df, df$maxcl==(i))
    # Select species
    df.cl.X$species <- as.character( df.cl.X$species )
    # Find index for cluster indval
    indCol <- grep(paste0("indval_",i,"$"), names(df.cl.X))
    # Find index for cluster frequency in dataframe (which column) 
    freqCol <- grep(paste0("freq_",i,"$"), names(df.cl.X)) 
    # Order by indval
    df.cl.X <- df.cl.X[order(-df.cl.X[,indCol]), ]
    # Remove NA's, if any 
    df.cl.X <- df.cl.X[!is.na(df.cl.X$species), ]
    # Select species where indval value is greater than cutoff
    df.cl.XLim <- df.cl.X[df.cl.X[,indCol] >= indvalCutoff, ]
    # If no species meet the cutoff, then ...
    if(nrow(df.cl.XLim)==0){
      df.cl.XLim[1,c(1:2)] <- c( "no ind species",i )
      df.cl.XLim$indvalInMaxcl <- NA
      df.cl.XLim$freqInMaxcl <- NA
    # Else, select appropriate values 
    } else {
      df.cl.XLim$indvalInMaxcl <- unlist( c(df.cl.XLim[,indCol]) )
      df.cl.XLim$freqInMaxcl <- unlist( c(df.cl.XLim[,freqCol]) )
    }
    # Subset dataframe
    #df.cl.XLim <- df.cl.XLim[,c("species","maxcl","indvalInMaxcl","freqInMaxcl",names(abu[,2:ncol(abu)]))]
    df.cl.XLim <- df.cl.XLim[,c("species","maxcl","indvalInMaxcl","freqInMaxcl")]# Need to build a new sequence here!
    df.cl.XLim
    # Add to list
    topSp[[i]] <- df.cl.XLim
  }
  # Build dataframe from list
  topSpdf <- do.call("rbind",topSp)
  return(topSpdf)
}

# Pull out ALL dataframe from list and run function
indval.pvalMin <- indval.pvals$ALL
indSpp <- indicatorSpecies(indval.pvalMin, indvalCutoff)
indSpp

saveRDS(indSpp, "IndicatorSpecies.labdsv.RDS")
write_csv(indSpp, "IndicatorSpecies.labdsv.csv")

# Save the entire workspace
#--------------------------
save.image(file="my_work_space_labdsv.RData")
load("my_work_space_labdsv.RData")

# # Can I run it on each element of the list?
# for (i in 1:length(indval.pvals)){
#   df.pvalMin <- indval.pvals[[i]]
#   cl <- unique(indval.pvals[[i]]$maxcl)
#   indSpp <- lapply(df.pvalMin, indicatorSpecies)
#   newList[[i]] <- indSpp
# }


# # TESTING ## TESTING ## TESTING ## TESTING ## TESTING ## TESTING ## TESTING #
# #888888888888888888888888888888888888888888888888888888888888888888888888888#
# # TESTING ## TESTING ## TESTING ## TESTING ## TESTING ## TESTING ## TESTING #
# 
# str(indval.pvalMin) # starting point
# 
# listElement <- vector("list", 4)
# 
# # Function to select indicator species
# indicatorSpecies <- function(df, clusters = cl, indvalCutoff){
#   # New empty list
#   topSp <- list()
#   # df to hold list element i
#   # df <- listElement[[i]]
#   # cl.cnt <- length(unique(df$maxcl))
#   for (i in 1:length(clusters)){
#     print(i)
#     # Loop thru each cluster...
#     df.cl.X <- df[df$maxcl==(i), ]
#     #df.cl.X <- dplyr::filter(df, df$maxcl==(i))
#     # Select species
#     df.cl.X$species <- as.character( df.cl.X$species )
#     # Find index for cluster indval
#     indCol <- grep(paste0("indval_",i,"$"), names(df.cl.X))
#     # Find index for cluster frequency  
#     freqCol <- grep(paste0("freq_",i,"$"), names(df.cl.X)) 
#     # Order by indval
#     df.cl.X <- df.cl.X[order(-df.cl.X[,indCol]), ]
#     # Remove NA's, if any 
#     df.cl.X <- df.cl.X[!is.na(df.cl.X$species), ]
#     # Select species where indval value is less than cutoff
#     df.cl.XLim <- df.cl.X[df.cl.X[,indCol] >= indvalCutoff, ]
#     # If no species are less than the cutoff, then ...
#     if(nrow(df.cl.XLim)==0){
#       df.cl.XLim[1,c(1:2)] <- c( "no ind species",i )
#       df.cl.XLim$indvalInMaxcl <- NA
#       df.cl.XLim$freqInMaxcl <- NA
#       # Else, select appropriate values 
#     } else {
#       df.cl.XLim$indvalInMaxcl <- unlist( c(df.cl.XLim[,indCol]) )
#       df.cl.XLim$freqInMaxcl <- unlist( c(df.cl.XLim[,freqCol]) )
#       df.cl.XLim[,freqCol] <- NA
#     }
#     # Subset dataframe
#     #df.cl.XLim <- df.cl.XLim[,c("species","maxcl","indvalInMaxcl","freqInMaxcl",names(abu[,2:ncol(abu)]))]
#     df.cl.XLim <- df.cl.XLim[,c("species","maxcl","indvalInMaxcl","freqInMaxcl",paste0("freq_",seq(1:cl.cnt)))]# Need to build a new sequence here!
#     df.cl.XLim
#     # Add to list
#     topSp[[i]] <- df.cl.XLim
#   }
#   # Build dataframe from list
#   topSpdf <- do.call("rbind",topSp)
#   return(topSpdf)
# }
# indSpp <- vector("list", 4)
# indSpp <- indicatorSpecies(indval.pvalMin, indvalCutoff)
# 
# for (i in 1:length(indval.pvalMin)){
#   df.pvalMin <- indval.pvalMin[[i]]
#   cl <- unique(indval.pvalMin[[i]]$maxcl)
#   indSpp <- lapply(df.pvalMin, indicatorSpecies)
#   newList[[i]] <- indSpp
# }
# 
# # Order by maxcl
# # 
# 
# # # 5. Join with species look-up table
# # spLookup<-read.csv("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/LookupTbls/SpeciesLookUpTbl.csv")
# # #names(spLookup)[1]<-"Species_Code"
# # topSpdf<-merge(topSpdf, spLookup, by.x="species", by.y="Sp_cde")
# # # topSpdf<-topSpdf[,c(ncol(topSpdf), 2:(ncol(topSpdf)-1))]
# # # topSpdf<-topSpdf[order(topSpdf$maxcl, -topSpdf$indvalInMaxcl),]
# # head(topSpdf,3)
# # topSpdf <- topSpdf[c(1,11,12,2:10)]
# 
# # Save output
# saveRDS(indSpp, "IndicatorSpecies.RDS")
# 
