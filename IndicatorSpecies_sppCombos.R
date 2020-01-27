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
library(vegan)


# Inputs
#-------
nSiteCl <- 0.20 # Threshold of sites a species is present in for each cluster 
max.order <- 3  # Maximum number of species combinations
At <- 0.6       # Threshold for positive predictive value
Bt <- 0.25      # Threshold for sensitivity


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
#outdir <- paste0(date,"_",region)
outdir <- paste0("nSites",nSiteCl,"_maxCombo",max.order,"_At",At,"_Bt",Bt)
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


# Read in RDS files --- #1 & #2
#------------------
# Species by region
speciesFullCl <- readRDS("../../RDS/speciesFullCl.RDS")
# Species list
species <- readRDS("../../RDS/species.RDS")
# Cluster frequency
cluster.freq <- readRDS("../../RDS/cluster.freq.RDS")


# Determine which species occur in atleast _ _ % of target clusters --- #3
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
thrshld <- data.frame( cl=as.integer(cl.freq$cl), thrshld = round(nSiteCl * cl.freq$Freq) )

# Determine number of species that meet frequency threshold
newSpp <- left_join(df, thrshld, by="cl")
newSpp$IN <- (newSpp$nSites - newSpp$thrshld)
newSpp <- filter(newSpp, IN>=0)
targetSpp <- unique(newSpp$Spp)
targetSpp


# Build cluster vector and new community data table  
#--------------------------------------------------
row.names(sppAll) <- sppAll$TransDepth
clusters <- sppAll$cl
sppObs <- sppAll[,targetSpp]
sppObs[is.na(sppObs)] <- 0
head(sppObs, 3)


# Multipatt() multi-level pattern analysis 
#-----------------------------------------
indval = multipatt(x = sppObs, cluster = clusters, max.order = 4, control = how(nperm=999))
summary(indval)


# Determine if the frequency of species in each site group was higher or lower than random
#-----------------------------------------------------------------------------------------
# Test the association btw species & each group of sites
# uses psidak correction 
# What am I comparing this output too?
prefsign <- signassoc( sppObs, cluster=clusters, alternative = "two.sided", control = how(nperm=199) )
head(prefsign)

 
# Determine the quantity coverage of the site group
#--------------------------------------------------
# The proportion of sites of a given site group where one or another indicator is found
indvalori <- multipatt(sppObs, clusters, duleg = TRUE, control = how(nperm=999))
summary(indvalori)

# Calculate the proportion of sites of the target site group where one or another indicator is found
#---------------------------------------------------------------------------------------------------
coverage(sppObs,indvalori)
coverage(sppObs, indvalori, At = 0.4, alpha = 0.05) # bound coverage by A & p-value
 

# Plot how coverage changes with 'A' threshold
#---------------------------------------------
# Add code to save as .png
par(mfrow = c(1,1))
png("Coverage_noSppCombos.png")
plotcoverage(x=sppObs, y=indvalori, group="1", lty=1)
plotcoverage(x=sppObs, y=indvalori, group="2", lty=2, col="blue", add=TRUE)
plotcoverage(x=sppObs, y=indvalori, group="3", lty=3, col="red", add=TRUE)
plotcoverage(x=sppObs, y=indvalori, group="4", lty=3, col="green", add=TRUE)
#plotcoverage(x=sppObs, y=indvalori, group="5", lty=3, col="purple", add=TRUE)
legend(x = 0.01, y=30,legend=c("group 1","group 2","group 3","group 4"),
       lty=c(1,2,3), col=c("black","blue","red","green"), bty="n")
dev.off()

# Species combinations as indicators of site groups --- #4
#--------------------------------------------------
# Build matrix with all possible species combinations
# My computer can't handle max.order=4
sppComb <- combinespecies(sppObs, max.order = max.order)$XC #...slow
dim(sppComb)
saveRDS(sppComb, "sppComb.RDS")

# Re-run mulitpatt with species combinations ADD CI
indvalspcomb = multipatt(sppComb, clusters, duleg = TRUE, control = how(nperm=999))#...slow
# List species with a significant association to one combination, including indval components
summary(indvalspcomb, indvalcomp = TRUE) 
saveRDS(indvalspcomb, "SpeciesComboMultipatt.RDS")

# Coverage - proportion of sites where one or another indicator is found
# To Do: look at options for this function
coverageSC <- coverage(sppComb,indvalspcomb)
coverageSD <- coverage(sppComb, indvalspcomb, At=0.6, alpha = 0.05)
# Create output with A, B, or indval stat > some threshold to reduce output
#...

# Plot how coverage changes with 'A' threshold
#---------------------------------------------
# Add code to save as .png
par(mfrow = c(1,1))
png("Coverage_SppCombos.png")
plotcoverage(x=sppComb, y=indvalspcomb, group="1", lty=1)
plotcoverage(x=sppComb, y=indvalspcomb, group="2", lty=2, col="blue", add=TRUE)
plotcoverage(x=sppComb, y=indvalspcomb, group="3", lty=3, col="red", add=TRUE)
plotcoverage(x=sppComb, y=indvalspcomb, group="4", lty=3, col="green", add=TRUE)
#plotcoverage(x=sppObs, y=indvalori, group="5", lty=3, col="purple", add=TRUE)
legend(x = 0.01, y=30,legend=c("group 1","group 2","group 3","group 4"),
       lty=c(1,2,3), col=c("black","blue","red","green"), bty="n")
dev.off()

# Determine indicators for each group 
#------------------------------------
# Determine strength of species site-group associations
# Square root of IndVal index is returned
B <- strassoc( sppObs, cluster=clusters, func="B" )

# Loop through clusters and select species with more than 20% of sensitivity for the first group
# Create empty lists for loop outputs
sc <- vector("list",length(unique(clusters)))
sc2 <- vector("list",length(unique(clusters)))

for (i in 1:length(unique(clusters))){
  label <- paste0("Group ",i)
  # Select species with more than 20% of sensitivity for group [[i]]
  sel <- which(B[,i]>0.2)
  names(sc)[[i]] <- i
  # Run indicator analysis with species combinations
  sc[[i]] <- indicators(X=sppObs[,sel], cluster=clusters, group=i, max.order=max.order, 
               verbose=TRUE, XC=TRUE, nboot=1000, At=At, Bt=Bt)
  print(sc[[i]]) 
  print(length(sel))
  # Determine if combinations improve coverage
  par(mfrow = c(1,1))
  filename <- paste0("Coverage.Group",i,".png")
  png(filename)
  plotcoverage(sc[[i]], main=label)
  plotcoverage(sc[[i]], max.order=1, add=TRUE, lty=2, col="red")
  legend(x=0.1, y=20, legend=c("Species combinations","Species singletons"),
         lty=c(1,2), col=c("black","red"), bty="n")
  dev.off()
  # Plot positive predictive power and sensitivity against the order of combinations
  filename <- paste0("PosPredPwr.Sensitivity.Grp",i,".png")
  png(filename)
  par(mfrow = c(1,2))
  plot(sc[[i]], type="A", main=label)
  plot(sc[[i]], type="B")
  dev.off()
  summary(sc[[i]])
  # Prune indicators to determine the best subset of indicators
  sc2[[i]] <- pruneindicators( sc[[i]], At=At, Bt=Bt, verbose=TRUE )
  names(sc2)[[i]] <- i
  print(sc2[[i]])
}

sc.4 <- sc[4]
new.cov <- coverage(sc.4, max.order=3)


# Clean up
par(mfrow = c(1,1))
saveRDS(sc, "Indicators4Clusters.RDS")
saveRDS(sc2, "PrunedIndicators4Clusters.RDS")
sc <- readRDS("Indicators4Clusters.RDS")
sc2 <- readRDS("PrunedIndicators4Clusters.RDS")
rm(cluster.freq) 

# 7.

# Need to build new dataframe for the output
#-------------------------------------------
i <- 1
i <- 2
# Table 1. Cl | Name | Sites | Spp | Ind | Valid | Final | Cover
# Cl = cluster number
# Name = descriptive name --- make it up
# Sites = number of sites (cluster.freq.RDS)
# Ind = number of candidate species --- what is candidate? from IndVal()?
# Valid = number of valid indicators
# Final = smallest set of valid indicators with the same coverage as the complete set 
# Cover = percentage coverage of the final set of valid indicators (output from pruneindicators())

nCl <- as.vector(integer())
nSites <- as.vector(integer())
nCandidate <- as.vector(integer())
nValid <- as.vector(integer())
nFinal <- as.vector(integer())
sppInd <- vector("list", nrow(cl.freq))

for (i in 1:nrow(cl.freq)){
  print(i)
  # 1 name - site group
  nCl[i] <- i
  # number of sites
  nSites[i] <- cl.freq$Freq[cl.freq$cl==(i)]
  # number of candidate species
  nCandidate[i] <- length(sc[[i]]$candidates)
  # Valid 
  df <- print(sc[[i]])
  sppInd[[i]] <- row.names(df)
  nValid[i] <- length(sppInd)
  # Final - output from prune
  df <- print(sc2[[i]])
  finalInd <- row.names(df)
  nFinal[i] <- length(finalInd)
}

site.Details <- as.data.frame( cbind(nCl, nSites, nCandidate, nValid, nFinal, coverageSC) )
colnames(site.Details) <- c("Cluster","Num.Sites","Num.Candidate.Spp","Num.Valid.Spp","Num.Final.Spp","Coverage")
write_csv(site.Details, path="SiteCharacteristics.csv")
saveRDS(site.Details, "SiteCharacteristics.RDS")

# Table 2. Cl | Name | Indicators | A (95CI) | B (95CI) | sqrt(IV)(95CI)
#-----------------------------------------------------------------------
valIndLst <- vector("list", nrow(cl.freq) )

for (i in 1:length(sc)){
  # 1 name - site group
  names(valIndLst)[i] <- i
  print(i)
  df <- print(sc[[i]])
  sppInd <- row.names(df)
  pos.predict <- sc[[i]]$A %>% 
    rename(A = stat,
           A.lci = lowerCI,
           A.uci = upperCI)
  sensitivity <- sc[[i]]$B %>% 
    rename(B = stat,
           B.lci = lowerCI,
           B.uci = upperCI)
  sqrtIV <- sc[[i]]$sqrtIV %>% 
    rename(sqrtIV = stat)
  valIndLst[[i]]$final <- cbind(sppInd, pos.predict, sensitivity, sqrtIV)
  write_csv(valIndLst[[i]]$final, path=paste0("ValidInd.Grp",i,".csv") )
  print(nrow(valIndLst[[i]]$final) )
}
saveRDS(valIndLst, "ValidIndicatorsTbl.RDS")

# Indicator power
#----------------
# IP values
ip <- indpower(sppObs)
diag(ip) <- NA
# And TIP values
(TIP <- rowMeans(ip, na.rm = T))
## p value calculation for a species from Halme et al. 2009
## i is ID for the species
i <- 1
fun <- function(x, i) indpower(x)[i,-i]
## 'c0' randomizes species occurrences
os <- oecosimu(sppObs, fun, "c0", i=i, nsimul=99)
## get z values from oecosimu output
z <- os$oecosimu$z
## p-value
(p <- sum(z) / sqrt(length(z)))
## 'heterogeneity' measure
(chi2 <- sum((z - mean(z))^2))
pchisq(chi2, df=length(z)-1)
## Halme et al.'s suggested output
out <- c(TIP=TIP[i], 
         significance=p,
         heterogeneity=chi2,
         minIP=min(fun(sppObs, i=i)),
         varIP=sd(fun(sppObs, i=i)^2))
out


# To Do *1:  Compare final with print(sc[[i]])
#            Look in examples
#       *2:  Build plots
#            Add an order column, Fig 4 (how many indicators are combos? etc.)

# Figure 4: Coverage of the cluster types obtained when considering indicators




# ## Run indicator analysis with species combinations for the first group,
# ## but forcing 'Orysp' to be in all combinations
# sc2= indicators(X=sppObs[,sel], cluster=clusters, group=1, verbose=TRUE, At=At, Bt=Bt, enableFixed=TRUE)
# 
# sc= indicators(X=sppComb[,sel], cluster=clusters, group=2, max.order = 2, verbose=TRUE, At=0.5, Bt=0.3)
# print(sc, sqrtIVt = 0.02) 
# summary(sc[[i]])
# 
# signassoc(X=sppComb[,sel], cluster=clusters, mode=1, control = how(nperm=999))
# ## Look for species whose abundance is significantly higher in sites belonging
# ## to one group as opposed to sites not belonging to it.
# signassoc(X=sppComb[,sel], cluster=clusters, mode=0, control = how(nperm=999))
# 
# 
# # Prune indicators to determine if coverage is changed with a smaller set of indicators
# # output does not match example in tutorial
# sc2 <- pruneindicators(sc, At=0.5, Bt=0.2, verbose=TRUE)
# print(sc2)
# 
# # predict indicators
# pcv <- predict(sc2, sppObs, cv=TRUE)
# pcv1 <- predict(sc, cv=TRUE)
# 
# # Compared predicted probabilities for each site --- ???
# data.frame(Group1 = as.numeric(speciesFullCl$ALL$cl==1), Prob = pcv, Prob_CV = pcv)
# 


