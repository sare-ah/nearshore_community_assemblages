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
library(indicspecies) # multipatt() multi-level pattern analysis
library(rstudioapi)

# Start timer
#------------
#start.time <- Sys.time()

# Inputs
#-------
nSiteCl <- 0.20 # Threshold of sites a species is present in for each cluster 
max.order <- 3  # Maximum number of species combinations
At <- 0.6       # Threshold for positive predictive value (specificity)
Bt <- 0.25      # Threshold for sensitivity (fidelity)
region <- "ALL"


# Get path for this script
#-------------------------
# Set working directory to one above script directory
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
getwd()

# Source functions
#-----------------
#source('functions.R')

# Set up new directory for all results
#-------------------------------------
date <- format(Sys.Date(), "%b_%d")
#outdir <- paste0(date,"_",region)
outdir <- paste0(region,".nSites",nSiteCl,"_maxCombo",max.order,"_At",At,"_Bt",Bt)
setwd( "../Results") 
if (dir.exists(outdir)){
  setwd(outdir)
} else{
  dir.create(path = outdir)
  setwd(outdir)
}
getwd()


# Read in RDS files --- #1 & #2
#------------------
# Species by region
speciesFullCl <- readRDS("../../RDS/speciesFullCl.RDS")

# Species list
#species <- readRDS("../../RDS/species.RDS")
# Cluster frequency
cluster.freq <- readRDS("../../RDS/cluster.freq.RDS")


# Determine which species occur in atleast _ _ % of target clusters --- #3
#----------------------------------------------------------------
# # Grab just one list element
spp <- speciesFullCl[region]
spp <- map_dfr(spp,identity)
d <- dim(spp)
maxCols <- d[2]-2 # don't count TransDepth & cluster

cl.freq <- cluster.freq[region]
cl.freq <- map_dfr(cl.freq,identity)

# Create a long pivot table of all species present in each trans
spp_long <- pivot_longer(spp, cols = 1:maxCols, names_to = "Spp", values_to = "Presence")
spp_long <- dplyr::filter( spp_long, Presence==1)

# Calculate number of sites for each cluster/species combination - reduce size of species combo matrix
df <- as.data.frame( spp_long %>%
  group_by(cl, Spp) %>%
  summarise(nSites = length(TransDepth)) )

# Threshold of sites for each cluster
thrshld <- data.frame( cl=as.integer(cl.freq$cl), thrshld = round(nSiteCl * cl.freq$Freq) )
thrshld

# Determine number of species that meet frequency threshold
# Check what happens to species with low in one cl and high in another
newSpp <- left_join(df, thrshld, by="cl")
newSpp$IN <- (newSpp$nSites - newSpp$thrshld)
newSpp <- filter(newSpp, IN>=0)
targetSpp <- unique(newSpp$Spp)
sort(targetSpp)


# Build cluster vector and new community data table  
#--------------------------------------------------
row.names(spp) <- spp$TransDepth
clusters <- spp$cl
sppObs <- spp[,targetSpp]
sppObs[is.na(sppObs)] <- 0
head(sppObs, 3)


# Determine the quantity coverage of the site group using multi-level pattern analysis 
#-------------------------------------------------------------------------------------
# The proportion of sites of a given site group where one or another indicator is found
indvalori <- multipatt(sppObs, clusters, duleg = TRUE, control = how(nperm=999))
summary(indvalori)
print(indvalori)

 
# Calculate the proportion of sites of the target site group where one or another indicator is found
#---------------------------------------------------------------------------------------------------
coverage(sppObs,indvalori)
coverage(sppObs, indvalori, At = At, Bt=Bt, alpha = 0.05) # bound coverage by At, Bt, & p-value
 

# Plot how coverage changes with 'A' threshold
#---------------------------------------------
par(mfrow = c(1,1))
png("Coverage_noSppCombos.png")
plotcoverage(x=sppObs, y=indvalori, group="1", lty=1)
plotcoverage(x=sppObs, y=indvalori, group="2", lty=2, col="blue", add=TRUE)
plotcoverage(x=sppObs, y=indvalori, group="3", lty=1, col="green", add=TRUE)
plotcoverage(x=sppObs, y=indvalori, group="4", lty=2, col="red", add=TRUE)
#plotcoverage(x=sppObs, y=indvalori, group="5", lty=3, col="purple", add=TRUE)
legend(x = 0.01, y=30,legend=c("group 1","group 2","group 3","group 4"),
       lty=c(1,2), col=c("black","blue","green","red"), bty="n")
dev.off()

# Species combinations as indicators of site groups --- #4
#--------------------------------------------------
# Build matrix with all possible species combinations
# My computer can't handle max.order=4
sppComb <- combinespecies(sppObs, max.order = max.order)$XC #...slow
dim(sppComb)
saveRDS(sppComb, "sppComb.RDS")
sppComb <- readRDS("sppComb.RDS")

# Re-run mulitpatt with species combinations ADD CI
indvalspcomb = multipatt(sppComb, clusters, duleg = TRUE, control = how(nperm=999))#...slow
# List species with a significant association to one combination, including indval components
summary(indvalspcomb, indvalcomp = TRUE) 
saveRDS(indvalspcomb, "SpeciesComboMultipatt.RDS")
indvalspcomb <- readRDS("SpeciesComboMultipatt.RDS")

# Coverage - proportion of sites where one or another indicator is found
coverageSC <- coverage(sppComb,indvalspcomb)
coverageSD <- coverage(sppComb, indvalspcomb, At=0.6, Bt=Bt, alpha = 0.05)


# Plot how coverage changes with 'A' threshold
#---------------------------------------------
# Use this figure to determine the appropriate At threshold value to use
par(mfrow = c(1,1))
png("Coverage_SppCombos.png")
plotcoverage(x=sppComb, y=indvalspcomb, group="1", lty=1)
plotcoverage(x=sppComb, y=indvalspcomb, group="2", lty=2, col="blue", add=TRUE)
plotcoverage(x=sppComb, y=indvalspcomb, group="3", lty=1, col="green", add=TRUE)
plotcoverage(x=sppComb, y=indvalspcomb, group="4", lty=2, col="red", add=TRUE)
#plotcoverage(x=sppObs, y=indvalori, group="5", lty=3, col="purple", add=TRUE)
legend(x = 0.01, y=30,legend=c("group 1","group 2","group 3","group 4"),
       lty=c(1,2), col=c("black","blue","green","red"), bty="n")
dev.off()

# palette="Dark2"

# Determine indicators for each group 
#------------------------------------
# Determine strength of species site-group associations
# Square root of IndVal index is returned
B <- strassoc( sppObs, cluster=clusters, func="B" )

# Loop through clusters and select species with more than 20% of sensitivity for the first group
# Create empty lists for loop outputs
sc <- vector("list",length(unique(clusters)))
sc2 <- vector("list",length(unique(clusters)))
indCombs <- vector("list",length(unique(clusters)))

for (i in 1:length(unique(clusters))){
  label <- paste0("Group ",i)
  # Select species with more than 20% of sensitivity for group [[i]]
  sel <- which(B[,i]>0.2)
  names(sc)[[i]] <- i
  X=sppObs[,sel]
  # Run indicator analysis with species combinations
  sc[[i]] <- indicators(X=X, cluster=clusters, group=i, max.order=max.order, 
               verbose=TRUE, XC=TRUE, nboot=1000, At=At, Bt=Bt)
  #indCombs[[i]] <- colnames(X) # Species used in species combos
  indCombs[[i]] <- data.frame(print(sc[[i]])) # Species and species combos
  # df <- print(sc[[i]])
  # Indicator <- as.character(row.names(df))
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

# Clean up
par(mfrow = c(1,1))
saveRDS(sc, "Indicators4Clusters.RDS")
saveRDS(sc2, "PrunedIndicators4Clusters.RDS")
saveRDS(indCombs, "IndicatorCombos.RDS")
sc <- readRDS("Indicators4Clusters.RDS")
sc2 <- readRDS("PrunedIndicators4Clusters.RDS")
indCombs <- readRDS("IndicatorCombos.RDS")
rm(cluster.freq) 


# # Change in coverage and number of species per indicator
# #-------------------------------------------------------
# # Select one element with a list of lists (nested list element selection)
# sens <- lapply(sc, '[',7)
# # indNmes <- sapply(indCombs, "[",1:4)
# 
# # Build empty dataframe
# i = 4
# #indCombs <- as.vector(c("LT", "AM", "AC+LT", "AC+NT"))
# 
# # indCov <- tibble("Cluster" = as.integer(),
# #                  "indSp" = as.character(),
# #                  "typeCnt" = as.numeric(),
# #                  "B" = as.numeric())
# # new.row <- indCov
# 
# new.row <- tibble("Cluster" = as.integer(),
#                  "indSp" = as.character(),
#                  "typeCnt" = as.numeric(),
#                  "B" = as.numeric())
# indCov <- vector("list",4)
# 
# 
# # Loop through each cluster and populate
# for (i in 1:length(sc)){
#   # Create a list of vectors for each new column for this cluster
#   col_vectors <- list(
#     indSp = row.names(indCombs[[i]]),
#     B = (100*sc[[i]]$B$stat) )
#   # Bind the list of vectors into new rows in a dataframe
#   new.row <- bind_rows(col_vectors)
#   # Add Cluster number
#   new.row$Cluster <- i
#   # Calculate Indicator type
#   new.row$typeCnt <- (1 + (str_count(new.row$indSp, pattern='\\+')) )
#   #new.row$typeCnt <- paste0(new.row$typeCnt, " Species" )
#   #indCov <- bind_rows(indCov, new.row)
#   indCov[[i]] <- new.row
# }
# indCov
# 
# 
# 
# 
# covAll <- vector("list",4)
# #coverageSD <- coverage(sppComb, indvalspcomb, At=0.6, Bt=Bt, alpha = 0.05)
# # 1         2         3         4 
# # 0.7374582 0.4233171 0.9925373 0.6380449 
# 
# for (i in 1:length(sc)){
#   # Build logic vectors
#   df <- indCov[[i]]
#   df <- df %>% 
#           mutate(ind.1 = ifelse(df$typeCnt==1, TRUE, FALSE)) %>%
#           mutate(ind.2 = ifelse(df$typeCnt==2, TRUE, FALSE)) %>%
#           mutate(ind.3 = ifelse(df$typeCnt==3, TRUE, FALSE))
#   covAll[[i]]$combo1 <- coverage(sc[[i]], indvalspcomb, selection=as.vector(df$ind.1), At=0.6, Bt=Bt, alpha = 0.05, type="lowerCI")
#   covAll[[i]]$combo2 <- coverage(sc[[i]], indvalspcomb, selection=as.vector(df$ind.2), At=0.6, Bt=Bt, alpha = 0.05, type="lowerCI")
#   covAll[[i]]$combo3 <- coverage(sc[[i]], indvalspcomb, selection=as.vector(df$ind.3), At=0.6, Bt=Bt, alpha = 0.05, type="lowerCI")
#   covAll[[i]]$Sum <- covAll[[i]]$combo3 + (covAll[[i]]$combo3 - covAll[[i]]$combo2) + (covAll[[i]]$combo2 - covAll[[i]]$combo3)
# }
# # Coverage - proportion of sites where one or another indicator is found
# coverageSC <- coverage(sppComb,indvalspcomb)
# coverageSD <- coverage(sppComb, indvalspcomb, At=0.6, Bt=Bt, alpha = 0.05)
# 
# cov.new <- coverage(sc, indvalspcomb, At=0.6, Bt=Bt, alpha = 0.05) 
# 
# 
# plotcoverage(sc[[i]], main=label)
# plotcoverage(sc[[i]], max.order=1, add=TRUE, lty=2, col="red", alpha = 0.05)
# legend(x=0.1, y=20, legend=c("Species combinations","Species singletons"),
#        lty=c(1,2), col=c("black","red"), bty="n")
# 
# 
# plot(sc[[i]], type="A", main=label)
# plot(sc[[i]], type="B")
# 
# 
# coverageSC <- coverage(sppComb,indvalspcomb)
# coverageSD <- coverage(sppComb, indvalspcomb, At=0.6, Bt=Bt, alpha = 0.05)
# covIndicators <- coverage(sc[[2]], indvalspcomb, selection=sc[[2]]$group.vec) #, At=0.6, Bt=Bt, alpha = 0.05)
# 
# 
# plotcoverage(indCombs)
# plotcoverage(indCombs, max.order=2, add=TRUE, lty=2)
# 
# 
# # Subset sppComb by species indicator vector
# for (i in 1:length(indCov)){
#   sel <- indCov[[i]]$indSp
#   newX <- sppComb[,sel]
#   coverageX[[i]] <- coverage(sppComb, indvalspcomb, selection=sel) #, At=0.6, Bt=Bt, alpha = 0.05)
# }

plot.indicators(sc[[2]])#x, type="sqrtIV", maxline=TRUE, ...)
# What is x$C - 

# 7.

# Need to build new dataframe for the output
#-------------------------------------------
# Table 1. Cl | Name | Sites | Spp | Ind | Valid | Final | Cover
nCl <- as.vector(integer()) #1 cluster 
nSites <- as.vector(integer()) #3 Number of sites (cluster.freq.RDS)
nCandidate <- as.vector(integer()) #4 Number of candidate species - what is candidate? from IndVal()?
nValid <- as.vector(integer()) #6 Number of valid indicators
nFinal <- as.vector(integer()) #7 Smallest set of valid indicators w/same coverage as the complete set (after pruning)
sppInd <- vector("list", nrow(cl.freq)) # List of valid indicators (single species + valid species combinations)
finalInd <- vector("list", nrow(cl.freq)) # List of final indicators (single species + valid species combinations)

for (i in 1:nrow(cl.freq)){
  print(i)
  # Cluster number
  nCl[i] <- i
  # Sites
  nSites[i] <- cl.freq$Freq[cl.freq$cl==(i)]
  # Candidate species
  nCandidate[i] <- length(sc[[i]]$candidates)
  # Valid indicators (species + species combinations)
  df <- print(sc[[i]])
  sppInd[[i]] <- row.names(df)
  nValid[i] <- length(sppInd[[i]])
  # Final indicators - output from prune()
  df <- print(sc2[[i]])
  finalInd[[i]] <- row.names(df)
  nFinal[i] <- length(finalInd[[i]])
}

site.Details <- as.data.frame( cbind(nCl, nSites, nCandidate, nValid, nFinal, coverageSC) )
colnames(site.Details) <- c("Cluster","Sites","Candidate.Spp","Valid.Spp+SppCombos","Pruned.Spp+SppCombos","Coverage")
site.Details
write_csv(site.Details, path="SiteCharacteristics.csv")
saveRDS(site.Details, "SiteCharacteristics.RDS")

# Table 2. Cl | Name | Indicators | A (95CI) | B (95CI) | sqrt(IV)(95CI)
#-----------------------------------------------------------------------
# Build empty dataframe to populate
valInd <- tibble("Cluster" = as.integer(),
                     "Indicator" = as.character(),
                     "A.lCI" = as.numeric(),
                     "A" = as.numeric(),
                     "A.uCI" = as.numeric(),
                     "B.lCI" = as.numeric(),
                     "B" = as.numeric(),
                     "B.uCI" = as.numeric(),
                     "sqrtIV.lCI" = as.numeric(),
                     "sqrtIV" = as.numeric(),
                     "sqrtIV.uCI" = as.numeric())
# Loop through each cluster and calculate columns for table
for (i in 1:length(sc)){
  Cluster <- i
  df <- print(sc[[i]])
  Indicator <- as.character(row.names(df))
  pos.predict <- sc[[i]]$A %>% 
    select(stat, lowerCI, upperCI) %>%
    rename(A.lCI = lowerCI,
           A = stat,
           A.uCI = upperCI)
  sensitivity <- sc[[i]]$B %>% 
    select(stat, lowerCI, upperCI) %>%
    rename(B.lCI = lowerCI,
           B = stat,
           B.uCI = upperCI)
  sqrtIV <- sc[[i]]$sqrtIV %>% 
    select(stat, lowerCI, upperCI) %>%
    rename(sqrtIV.lCI = lowerCI,
           sqrtIV = stat,
           sqrtIV.uCI = upperCI)
  new.row <- cbind(Cluster, Indicator, pos.predict, sensitivity, sqrtIV, stringsAsFactors=FALSE)
  valInd <- bind_rows(valInd, new.row)
}
valInd
df.1 <- dplyr::filter(valInd, Cluster==2)
write_csv(valInd, "ValidInd.Grp.csv")
write_csv(df.1, "ValidInd.Grp2.csv")
saveRDS(valInd, "ValidIndicatorsTbl.RDS")
valInd <- readRDS("ValidIndicatorsTbl.RDS")




new.df <- pivot_longer(valInd,col=3:11, names_to = "stat", values_to = "values") 
new.df$Cluster <- as.factor(new.df$Cluster)

# To add to plots:
# Do I want to plot sqrtIV or mutate to IndVal?

#png("PositivePredPwr_byCluster.png", res=72, height=800, width=1200)
A.plot <- new.df %>%
  filter(str_detect(stat,"A")) %>% 
  ggplot( aes(x=Indicator, y=values, fill=Cluster)) +
    geom_boxplot() +
    ylab("Positive Predictive Power (A)") +
    facet_wrap(~Cluster, scales = "free_x") +
    scale_fill_brewer(palette="Dark2") +
    theme_light() +
    theme(strip.text.x = element_text(color="black")) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) 
A.plot
ggsave("PositivePredPwr_byCluster.png", A.plot, dpi = 300, width = 12, height = 12, path = getwd())
  

B.plot <- new.df %>%
  filter(str_detect(stat,"B")) %>%
  ggplot( aes(x=Indicator, y=values, fill=Cluster)) +
    geom_boxplot() +
    ylab("Sensitivity (B)") +
    facet_wrap(~Cluster, scales = "free_x") +
    scale_fill_brewer(palette="Dark2")+
    theme_light() +
    theme(strip.text.x = element_text(color="black")) +
    theme(axis.text.x = element_text(angle=90, hjust=1))
B.plot
ggsave("Sensitivity_byCluster.png", A.plot, dpi = 300, width = 12, height = 12, path = getwd())


I.plot <- new.df %>%
  filter(str_detect(stat,"sqrt")) %>% 
  ggplot( aes(x=Indicator, y=values, fill=Cluster)) +
    geom_boxplot() +
    ylab("Square root Indicator Value") +
    facet_wrap(~Cluster, scales = "free_x") +
    scale_fill_brewer(palette="Dark2") +
    theme_light() +
    theme(strip.text.x = element_text(color="black")) +
    theme(axis.text.x = element_text(angle=90, hjust=1))
I.plot
ggsave("SqrtIndVal_byCluster.png", A.plot, dpi = 300, width = 12, height = 12, path = getwd())


plotcoverage(sppComb,indvalspcomb, by=0.5, type="stat")
# ggplot(coverageSC, aes(fill=typeCnt, y=B, x=Cluster)) + 
#   geom_bar(position="stack", stat="identity")

# Time taken to run code
#-----------------------
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


# # All metrics for each indicator in one cluster
# new.df <- new.df %>%
#   mutate(stat = replace(stat, str_detect(stat, "A"), "A")) %>% 
#   mutate(stat = replace(stat, str_detect(stat, "B"), "B")) %>%
#   mutate(stat = replace(stat, str_detect(stat, "sqrtIV"), "sqrtIV")) 
# 
# cl.plot <- dplyr::filter(new.df, Cluster==1)
# cl.plot$Indicator <- with(cl.plot, reorder(Indicator, values)) 
# cl.plot %>%
#   ggplot( aes(x=Indicator, y=values, fill=stat)) +
#   geom_boxplot() +
#   scale_fill_brewer(palette="Dark2") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle=90, hjust=1)) +
#   theme(panel.grid.major.y = element_line(colour = "grey"))



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



# # # predict indicators
# p <- predict(sc2[1], sppObs)
# pcv <- predict(sc2, sppObs, cv=TRUE)
# pcv1 <- predict(sc, cv=TRUE)
# 
# # Compared predicted probabilities for each site --- ???
# data.frame(Group1 = as.numeric(speciesFullCl$ALL$cl==1), Prob = pcv, Prob_CV = pcv)
# 

#**********************#
# Extras from tutorial
#**********************#

# # Determine if the frequency of species in each site group was higher or lower than random
# #-----------------------------------------------------------------------------------------
# # Test the association btw species & each group of sites
# # uses psidak correction 
# # What am I comparing this output too?
# prefsign <- signassoc( sppObs, cluster=clusters, alternative = "two.sided", control = how(nperm=199) )
# head(prefsign)







# # Indicator power
# #----------------
# # IP values
# ip <- indpower(sppObs)
# diag(ip) <- NA
# # And TIP values
# (TIP <- rowMeans(ip, na.rm = T))
# ## p value calculation for a species from Halme et al. 2009
# ## i is ID for the species
# i <- 1
# fun <- function(x, i) indpower(x)[i,-i]
# ## 'c0' randomizes species occurrences
# os <- oecosimu(sppObs, fun, "c0", i=i, nsimul=99)
# ## get z values from oecosimu output
# z <- os$oecosimu$z
# ## p-value
# (p <- sum(z) / sqrt(length(z)))
# ## 'heterogeneity' measure
# (chi2 <- sum((z - mean(z))^2))
# pchisq(chi2, df=length(z)-1)
# ## Halme et al.'s suggested output
# out <- c(TIP=TIP[i], 
#          significance=p,
#          heterogeneity=chi2,
#          minIP=min(fun(sppObs, i=i)),
#          varIP=sd(fun(sppObs, i=i)^2))
# out

plot.indicators<-function(x, type="sqrtIV", maxline=TRUE, ...) {
  A = x$A
  B = x$B
  sqrtIV=x$sqrtIV
  order = rowSums(x$C)
  if(is.data.frame(A)) {
    if(type=="IV") val = sqrtIV[,1]^2
    else if(type=="sqrtIV") val = sqrtIV[,1]
    else if(type=="A") val = A[,1]	
    else if(type=="B") val = B[,1]	
    else if(type=="LA") val = A[,2]	
    else if(type=="UA") val = A[,3]	
    else if(type=="LB") val = B[,2]	
    else if(type=="UB") val = B[,3]	
    else if(type=="LsqrtIV") val = sqrtIV[,2]	
    else if(type=="UsqrtIV") val = sqrtIV[,3]	
  } else {
    if(type=="IV") val = sqrtIV^2
    else if(type=="sqrtIV") val = sqrtIV
    else if(type=="A") val = A	
    else if(type=="B") val = B	   	
  }
  plot(order,val, type="n", axes=FALSE, xlab="Order", ylab=type,...)	
  points(order,val, pch=1, cex=0.5)	
  axis(1, at = order, labels=order)
  axis(2)
  if(maxline) lines(1:max(order),tapply(val,order,max), col="gray")
}
