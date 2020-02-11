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
library(indicspecies) # use indicspecies_1.7.6.tar.gz (older version)
library(labdsv)
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
outdir <- "Cluster8"
#outdir <- paste0(date,region,".nSites",nSiteCl,"_maxCombo",max.order,"_At",At,"_Bt",Bt)
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
speciesFullCl <- readRDS("speciesFullCl.RDS") # manually move correct matrix to working folder
#speciesFullCl <- readRDS("../../RDS/speciesFullCl.RDS")

# Species list
#species <- readRDS("../../RDS/species.RDS")
# Cluster frequency
cluster.freq <- readRDS("cluster.freq.RDS")


# Determine which species occur in atleast _ _ % of target clusters --- #3
#----------------------------------------------------------------
# # Grab just one list element
spp <- as.data.frame(speciesFullCl$ALL)
#spp <- map_dfr(spp,identity)
d <- dim(spp)
maxCols <- d[2]-2 # don't count TransDepth & cluster

cl.freq <- cluster.freq$ALL
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
row.names(spp) <- spp$TransDepth # To Do: row.names() on tibble is deprecated
clusters <- spp$cl
sppObs <- spp[,targetSpp]
sppObs[is.na(sppObs)] <- 0
head(sppObs, 3)

saveRDS(sppObs, "SpeciesMtrxUsed.RDS")
sppObs <- readRDS("SpeciesMtrxUsed.RDS")


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
# To do: Edit figure for number of clusters
par(mfrow = c(1,1))
png("Coverage_noSppCombos.png")
plotcoverage(x=sppObs, y=indvalori, group="1", lty=1)
plotcoverage(x=sppObs, y=indvalori, group="2", lty=2, col="blue", add=TRUE)
plotcoverage(x=sppObs, y=indvalori, group="3", lty=1, col="green", add=TRUE)
plotcoverage(x=sppObs, y=indvalori, group="4", lty=2, col="red", add=TRUE)
plotcoverage(x=sppObs, y=indvalori, group="5", lty=3, col="purple", add=TRUE)
#plotcoverage(x=sppObs, y=indvalori, group="6", lty=3, col="grey", add=TRUE)
legend(x = 0.01, y=30,legend=c("group 1","group 2","group 3","group 4","group 5"),
       lty=c(1,2,3), col=c("black","blue","green","red","purple"), bty="n")
dev.off()

# Species combinations as indicators of site groups --- #4
#--------------------------------------------------
# Build matrix with all possible species combinations
# My computer can't handle max.order=4
start.time <- Sys.time()
sppComb <- combinespecies(sppObs, max.order = max.order)$XC #...slow
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken # ~ 2 minutes

dim(sppComb)
saveRDS(sppComb, "sppComb.RDS")
sppComb <- readRDS("sppComb.RDS")

# Re-run mulitpatt with species combinations ADD CI
start.time <- Sys.time()
indvalspcomb = multipatt(sppComb, clusters, duleg = TRUE, control = how(nperm=999))#...slow
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken # ~ 45 minutes to 4 hours

# List species with a significant association to one combination, including indval components
summary(indvalspcomb, indvalcomp = TRUE) 
saveRDS(indvalspcomb, "SpeciesComboMultipatt.RDS")
indvalspcomb <- readRDS("SpeciesComboMultipatt.RDS")

# Coverage - proportion of sites where one or another indicator is found
coverage(sppComb,indvalspcomb)
covSppComb <- coverage(sppComb, indvalspcomb, At=At, Bt=Bt, alpha = 0.05)
covSppComb

# Plot how coverage changes with 'A' threshold
#---------------------------------------------
# Use this figure to determine the appropriate At threshold value to use
par(mfrow = c(1,1))
png("Coverage_SppCombos.png")
plotcoverage(x=sppComb, y=indvalspcomb, group="1", lty=1)
plotcoverage(x=sppComb, y=indvalspcomb, group="2", lty=2, col="blue", add=TRUE)
plotcoverage(x=sppComb, y=indvalspcomb, group="3", lty=1, col="green", add=TRUE)
plotcoverage(x=sppComb, y=indvalspcomb, group="4", lty=2, col="red", add=TRUE)
plotcoverage(x=sppComb, y=indvalspcomb, group="5", lty=3, col="purple", add=TRUE)
#plotcoverage(x=sppComb, y=indvalspcomb, group="6", lty=3, col="grey", add=TRUE)
legend(x = 0.01, y=30,legend=c("group 1","group 2","group 3","group 4","group 5"),
       lty=c(1,2,3), col=c("black","blue","green","red","purple"), bty="n")
dev.off()

# palette="Dark2"

# lapply(list, function)
# # Convert list element to a dataframe
# make_df <- function(x){
#   as.data.frame(table(x))
# }
# colorcount <- lapply(cl.list, make_df )

# Determine indicators for each group 
#------------------------------------
# Determine strength of species site-group associations
# Square root of IndVal index is returned
B <- strassoc( sppObs, cluster=clusters, func="B" )
B$Spp <- row.names(B)
B <- B %>%
  dplyr::select(Spp, everything())
write_csv(B, "SqrootOfIndval.csv")

# Loop through clusters and select species with more than 20% of sensitivity for the first group
# Create empty lists for loop outputs


# Loop through clusters and select candidate species with more than 20% of sensitivity 
candidates <- vector("list", length(unique(clusters)))
for (i in 1:length(unique(clusters))){
  # Select species with more than 20% of sensitivity for group [[i]]
  sel <- which(B[,i]>0.2)
  candidates[[i]]=sppObs[,sel]
}

# Determine which candidate species (and species combinations) are indicators
sc <- vector("list",length(unique(clusters)))
sc <- vector("list", 18)
names(sc) <- seq(1:length(sc))
for (i in 1:length(sc)){
 sc[[i]]$name <- paste0("Cluster ",i) 
}
for (i in 1:length(sc)){
  tryCatch(
    expr = {
      sc[[i]] <- indicators(X=candidates[[i]], cluster=clusters, group=i, max.order=max.order, 
                            verbose=TRUE, XC=TRUE, nboot=1000, At=At, Bt=Bt)
      sc[[i]]$name <- paste0("Cluster ",i) 
      # sc.10 <- indicators(X=candidates[[i]], cluster=clusters, group=i, max.order=max.order, 
   #                         verbose=TRUE, XC=TRUE, nboot=1000, At=At, Bt=Bt)
    },
    error = function(e){
      message(" No indicators for cluster ", i)
      print(e)
    },
    finally = {
      message(" Finished assessing cluster ", i)
      names(sc)[[i]] <- i
      summary(sc[[i]])
      print(names(sc))
    }
  )
}


# Create new list for the rest of the code
sc.new <- list(sc[[3]],sc[[7]],sc[[8]],sc[[9]])
names(sc.new) <- c("3","7","8","9")
sc <- sc.new

# Create empty list for species and species combinations 
indCombs <- vector("list",length(unique(clusters)))

# Extract species and species combinations from results
indCombs <- lapply( seq_along(1:length(unique(clusters))), function(i) {
  getNmes=try(row.names(print(sc[[i]])),TRUE)
  if(isTRUE(class(getNmes)=="try-error")) { 
    return(NULL) } 
  else { return(getNmes) } } )

# Add cluster names to list elements
names(indCombs) <- seq(1:length(unique(clusters)))

#new.seq <- c(3,8,9,10)

# Determine if combinations improve coverage
lapply( seq_along(1:length(sc)), function(i) {
    par(mfrow = c(1,1))
    filename <- paste0("Coverage.Group",i,".png")
    label <- i
    png(filename)
    plotcoverage(sc[[i]], main=label)
    plotcoverage(sc[[i]], max.order=1, add=TRUE, lty=2, col="red")
    legend(x=0.1, y=20, legend=c("Species combinations","Species singletons"),
       lty=c(1,2), col=c("black","red"), bty="n")
    dev.off()
  }
)

# # Plot positive predictive power and sensitivity against the order of combinations
# lapply( seq_along(1:length(unique(clusters))), function(i) {
#   if (!is.null(sc[[i]])){
#     filename <- paste0("PosPredPwr.Sensitivity.Grp",i,".png")
#     label <- i
#     png(filename)
#     par(mfrow = c(1,1))
#     plot(sc[[i]], type="A", main=label)
#     plot(sc[[i]], type="B")
#     dev.off() }
#   else {
#     print("No indicators for this cluster")
#     }
#   }
# )

# Create empty list for pruned species and species combinations
sc2 <- vector("list",length(sc))
names(sc2) <- seq(1:length(sc))

sc2 <- lapply( seq_along(1:length(unique(clusters))), function(i) {
  prune <- try(pruneindicators( sc[[i]], At=At, Bt=Bt, verbose=TRUE ),TRUE)
  if(isTRUE(class(prune)=="try-error")) { 
    return(NULL) } 
  else { return(prune) } } )

# Clean up
par(mfrow = c(1,1))
saveRDS(sc, "Indicators4Clusters.RDS")
saveRDS(sc2, "PrunedIndicators4Clusters.RDS")
saveRDS(indCombs, "IndicatorCombos.RDS")
sc <- readRDS("Indicators4Clusters.RDS")
sc2 <- readRDS("PrunedIndicators4Clusters.RDS")
indCombs <- readRDS("IndicatorCombos.RDS")
rm(cluster.freq) 


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

nCl <- unique(cl.freq$cl)

# Sites
nSites <- unlist(lapply(seq_along(nCl), function(i) {
                 value <- try(cl.freq$Freq[cl.freq$cl==(i)], TRUE)
                 if(isTRUE(class(value)=="try-error")) { 
                   return(NULL) } 
                 else { return(value) } } ))

# Candidate species
nCandidate <- unlist(lapply(seq_along(nCl), function(i) {
                    value <- try(length(sc[[i]]$candidates), TRUE)
                    if(isTRUE(class(value)=="try-error")) { 
                      return(NULL) } 
                    else { return(value) } } ))

# Valid indicators (species + species combinations)
nValid <- unlist(lapply(seq_along(nCl), function(i) {
                value <- try(length(print(sc[[i]])), TRUE)
                if(isTRUE(class(value)=="try-error")) { 
                  return(NULL) } 
                else { return(value) } } ))

# Final indicators - output from prune()
nFinal <- unlist(lapply(seq_along(nCl), function(i) {
                value <- try(length(print(sc2[[i]])), TRUE)
                if(isTRUE(class(value)=="try-error")) { 
                return(NULL) } 
                else { return(value) } } ))

site.Details <- as.data.frame( cbind(nCl, nSites, nCandidate, nValid, nFinal) )
#site.Details <- map_df(site.Details, identity)
colnames(site.Details) <- c("Cluster","Sites","Candidate.Spp","Valid.Spp+SppCombos","Pruned.Spp+SppCombos")
site.Details <- site.Details[order(site.Details$Cluster),]
site.Details
site.Details <- cbind(site.Details, covSppComb)
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
  tryCatch(
    expr = {
      Cluster <- i
      df <- print(sc[[i]])
      Indicator <- as.character(row.names(df))
      pos.predict <- sc[[i]]$A %>% 
        dplyr::select(stat, lowerCI, upperCI) %>%
        rename(A.lCI = lowerCI,
           A = stat,
           A.uCI = upperCI)
      sensitivity <- sc[[i]]$B %>% 
        dplyr::select(stat, lowerCI, upperCI) %>%
        rename(B.lCI = lowerCI,
           B = stat,
           B.uCI = upperCI)
      sqrtIV <- sc[[i]]$sqrtIV %>% 
        dplyr::select(stat, lowerCI, upperCI) %>%
        rename(sqrtIV.lCI = lowerCI,
           sqrtIV = stat,
           sqrtIV.uCI = upperCI)
      new.row <- cbind(Cluster, Indicator, pos.predict, sensitivity, sqrtIV, stringsAsFactors=FALSE)
      valInd <- bind_rows(valInd, new.row)
      },
      error = function(e){
        message(" No indicators for cluster ", i)
        print(e)
      },
      finally = {
        message(" Finished assessing cluster ", i)
        valInd
      }
  )
}
valInd
write_csv(valInd, "ValidInd.Grp.csv")
saveRDS(valInd, "ValidIndicatorsTbl.RDS")
valInd <- readRDS("ValidIndicatorsTbl.RDS")


# Create figures of confidence intervals for each indicator & indicator statistic
#--------------------------------------------------------------------------------
# To add to plots:
# Do I want to plot sqrtIV or mutate to IndVal?
# Melt to create a long table with indicators & stats
new.df <- pivot_longer(valInd,col=3:11, names_to = "stat", values_to = "values") 
new.df$Cluster <- as.factor(new.df$Cluster)


#png("PositivePredPwr_byCluster.png", res=72, height=800, width=1200)
A.plot <- new.df %>%
  dplyr::filter(str_detect(stat,"A")) %>% 
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


# To Do *1:  Compare final with print(sc[[i]])
#            Look in examples
#       *2:  Build plots
#            Add an order column, Fig 4 (how many indicators are combos? etc.)



# Save the entire workspace
#--------------------------
save.image(file="my_work_space_indicator.RData")
load("my_work_space_indicator.RData")


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


# Indicator power
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



