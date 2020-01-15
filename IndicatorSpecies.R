# Indicator Species Analysis

# start fresh
rm(list=ls())

library(reshape)
#library(data.table)
library(foreign)
library(vegan)
library(rgdal)
library(labdsv) # for indicator species analysis

#outdir<-"./forClusterAnalysis_/"

pValCutoff <- 0.05
indvalCutoff <- 0.15


# Read in RDS file
speciesFullCl <- readRDS("speciesFullCl.RDS")
dat <- speciesFullCl$HG
# 1. Read in species list
species <- readRDS("../../RDS/species.RDS")
# 2. Match species in species list as speciesNew (species list changes from region to region)
speciesNew <- species[species%in% names(dat)]

head(dat,3)
datNew <- dat[,speciesNew]
datNew[is.na(datNew)] <- 0
summary(datNew)
sumRows <- rowSums(datNew)
sumCols <- colSums(datNew)

# 3. #Run indval analysis
indvalList <- indval(datNew[,speciesNew], dat$cl) 

# 4. Organize indval output into a dataframe

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

indval.df <- indicatorValues(indvalList)
head(indval.df,3)


# Select only species that have a pval less than the predetermined cutoff?
if(!is.na(pValCutoff)){
  intabSub <- intab[intab$pval <= pValCutoff,]
} else { intabSub <- intab}

# New empty list
topSp <- list()


for (i in unique(intabSub$maxcl)){
  # Loop thru each cluster...
  intabSubx <- intabSub[intabSub$maxcl==i, ]
  # Select species
  intabSubx$species <- as.character( intabSubx$species )
  # Select matching indval values
  indCol <- grep(paste0("indval_",i,"$"), names(intabSubx))
  # Select matching frequency values
  freqCol <- grep(paste0("freq_",i,"$"), names(intabSubx)) 
  # Order by indval
  intabSubx <- intabSubx[order(-intabSubx[,indCol]), ]
  # Remove NA's, if any 
  intabSubx <- intabSubx[!is.na(intabSubx$species), ]
  # Select species where indval value is less than cutoff
  intabSubxLim <- intabSubx[intabSubx[,indCol] >= indvalCutoff, ]
  # If no species are less than the cutoff, then ...
  if(nrow(intabSubxLim)==0){
    intabSubxLim[1,c(1:2)] <- c( "no ind species",i )
    intabSubxLim$indvalInMaxcl <- NA
    intabSubxLim$freqInMaxcl <- NA
  # Else, select appropriate values 
  } else {
    intabSubxLim$indvalInMaxcl <- unlist( c(intabSubxLim[,indCol]) )
    intabSubxLim$freqInMaxcl <- unlist( c(intabSubxLim[,freqCol]) )
    intabSubxLim[,freqCol] <- NA
  }
  # Subset dataframe
  intabSubxLim <- intabSubxLim[,c("species","maxcl","indvalInMaxcl","freqInMaxcl",names(abu[,2:ncol(abu)]))]
  intabSubxLim
  # Add to list
  topSp[[i]]<-intabSubxLim
}
# Build dataframe from list
topSpdf <- do.call("rbind",topSp)
topSpdf

# 5. Join with species look-up table



#dat<-read.dbf(paste0(outdir,"3/SpatialPointsWithSpecies_cl__2019.04_HG_tryWards.dbf"))
#species<-readRDS(paste0(outdir,"species.rdata")) # This does not exist?


# dat_melt<-pivot_longer(dat, cols=1:103,names_to = "species",values_to = "present")
# 
# dat_melt<-melt(dat[,names(dat) %in% c("cl","TransDepth",species)], id.vars=c("TransDepth","cl"))
# dat_cast<-cast(dat_melt, ID+cl~variable, fun="mean")
# 
# #Remove species that do not occur in any site ## FLAG in input data why this is occuring
# specsums <- colSums(dat_cast[,names(dat_cast) %in% species])
# names(specsums) <- names(dat_cast[,names(dat_cast) %in% species])
# zeroSpec <- names(specsums[specsums==0])
# dat_cast <- dat_cast[,!names(dat_cast)%in% zeroSpec]
# speciesNew <- species[species%in% names(dat_cast)]

# There is an error here.  Need to bring in legendcluster df
legendcluster <- readRDS(paste0(outdir, "/legendcluster.rds"))
dat_cast<-dat_cast[dat_cast$cl %in% legendcluster$cl,]

#Run indval analysis
ind <- indval(dat_cast[,speciesNew], dat_cast$cl) 

indv <- data.frame(species=row.names(ind$indval), round(ind$indval,3))
names(indv) <- c("species",paste0("indval_",gsub("X","",names(indv[,2:ncol(indv)]))))
abu <- data.frame(species=row.names(ind$relabu), round((ind$relabu*100),1))
names(abu) <-c ("species",paste0("freq_",gsub("X","",names(abu[,2:ncol(abu)]))))
cls <- gsub("indval_","",names(indv)[2:length(names(indv))])
maxcls <- data.frame(species=names(ind$maxcls), maxcl=as.numeric(as.character(cls[ind$maxcls]))) #maxcls is the INDEX OF the class each species has the maximum indicator value for

pval <- data.frame(species=names(ind$pval),pval= ind$pval) #the class each species has the maximum indicator value for
intab <- cbind(maxcls,pval=pval[,-1],indv[,-1], abu[,-1])


if(!is.na(pValCutoff)){
  intabSub <- intab[intab$pval<=pValCutoff,]
} else { intabSub<-intab}

topSp<-list()
for (i in unique(intabSub$maxcl)){
  
  intabSubx<-intabSub[intabSub$maxcl==i,]
  intabSubx$species<-as.character(intabSubx$species)
  indCol<-grep(paste0("indval_",i,"$"),names(intabSubx))
  freqCol<-grep(paste0("freq_",i,"$"),names(intabSubx)) 
  intabSubx<-intabSubx[order(-intabSubx[,indCol]),]
  intabSubx<-intabSubx[!is.na(intabSubx$species),]
  intabSubxLim<-intabSubx[intabSubx[,indCol]>=indvalCutoff,]
  
  if(nrow(intabSubxLim)==0){
    intabSubxLim[1,c(1:2)]<-c("no ind species",i)
    intabSubxLim$indvalInMaxcl<-NA
    intabSubxLim$freqInMaxcl<-NA
    
  } else {
    intabSubxLim$indvalInMaxcl<-unlist(c(intabSubxLim[,indCol]))
    intabSubxLim$freqInMaxcl<-unlist(c(intabSubxLim[,freqCol]))
    intabSubxLim[,freqCol]<-NA
  }
  intabSubxLim<-intabSubxLim[,c("species","maxcl","indvalInMaxcl","freqInMaxcl",names(abu[,2:ncol(abu)]))]
  intabSubxLim
  topSp[[i]]<-intabSubxLim
}

topSpdf<-do.call("rbind",topSp)
topSpdf

spLookup<-read.csv("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/LookupTbls/SpeciesLookUpTbl.csv")
#names(spLookup)[1]<-"Species_Code"
topSpdf<-merge(topSpdf, spLookup, by.x="species", by.y="Sp_cde")
# topSpdf<-topSpdf[,c(ncol(topSpdf), 2:(ncol(topSpdf)-1))]
# topSpdf<-topSpdf[order(topSpdf$maxcl, -topSpdf$indvalInMaxcl),]
head(topSpdf,3)
topSpdf <- topSpdf[c(1,11,12,2:10)]


write.csv(topSpdf, paste0(outdir,"indicatorSpecies_tryWards.csv"),row.names=F)


#write.csv(matFullCl, paste0(outdir,"matFullCl_tryWards.csv"), row.names=F)
