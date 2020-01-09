#######################################################################################################################
# Cluster analysis of Benthic Habitat Mapping Dive Survey Sites
#
# To do: change loops to lapply statements 


# Start fresh
rm(list=ls())

# Load packages
#--------------

library(simba)
library(dendroextras)
library(colorspace)
library(scales)
library(mapview)
library(rgdal)
library(sp)
library(dendextend)
library(rstudioapi)
library(tidyverse)

# Inputs
#-------

date <- format(Sys.Date(), "%b_%d")
region <- "All"
outdir <- paste0(date,"_",region)
dsn <- "SHP"

spThreshold <- 0.03 # Proportion of sites that a species must be found in
siteThreshold <- 3 # Minimum number of species / site
distance <- "simpson"
clusterMethod <- "ward.D2" 
seth <- 4 # number for cutoff height
nTopClusters <- 5 #integer for number of clusters


# Site by Species csv file
# myFile <- choose.files(default = "T:/Benthic_Habitat_Mapping/Data",
#                        caption = "Select site X species, with enviromental variables") 
# myFile <- "T:/Benthic_Habitat_Mapping/Data/Species by Site Matrices/qcsbyDepthCat_AllSpp.csv"

# Read in species by regions
#---------------------------
sppByRegion <- readRDS("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/RDS/sppByRegion.rds")
sppDF <- readRDS("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/RDS/sppByRegion_Dataframe.rds")

# Get path for this script
#-------------------------
# Set working directory to one above script directory
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
getwd()

# Source functions
#-----------------
source('CommunityAssemblages_functions.R')

# Set up new directory for all results
#-------------------------------------
setwd( "../Results") 
if (dir.exists(outdir)){
  setwd(outdir)
} else{
  dir.create(path = outdir)
  setwd(outdir)
  dir.create(dsn)
}
getwd()


# Import data
# -----------
# # Format is 1 row per sample (depth bin), with columns for species and environmental variables 
# matFull <- read.csv(myFile, header=T, sep=",", stringsAsFactors = F)
# head(matFull,3)

# Look-up table to match species codes to latin names
luTbl <- read.csv( "T:/Benthic_Habitat_Mapping/Data/Look-upTbls/SpeciesLookUpTbl.csv", header=T, sep=",", stringsAsFactors=F )

# Cluster Analysis
#-----------------

# --- Remove rare species and barren sites from each region --- #
# Create empty lists needed for outputs --- maybe we only need the listOfLists???
listOfLists <- vector("list", 4)
prep4Cluster <- vector("list", 4)
rareSpp <- vector("list", 4)
barrenSites <- vector("list", 4)

# Calculate summaries
for (i in 1:length(sppByRegion)){
  # Get region name
  names(listOfLists)[[i]] <- names(sppByRegion)[i]
  spp <- as.data.frame(sppByRegion[[i]])
  # Build list of all speciesXsites, rare species, barren sites
  listOfLists[[i]]$sppSiteList <- rare( speciesObs=spp, minSp=spThreshold, minSite=siteThreshold )
  # Build new list of species P/A for cluster analysis
  names(prep4Cluster)[[i]] <- names(sppByRegion)[i]
  prep4Cluster[[i]] <- listOfLists[[i]]$sppSiteList[1]
  # Build list of rare species 
  names(rareSpp)[[i]] <- names(sppByRegion)[i]
  rareSpp[[i]] <- listOfLists[[i]]$sppSiteList[2]
  # Build list of barren sites 
  names(barrenSites)[[i]] <- names(sppByRegion)[i]
  barrenSites[[i]] <- listOfLists[[i]]$sppSiteList[3]
}

### TO DO ###
# Outputs: Table of rare species by region (spp | NCC | HG | QCS | SoG)
# Outputs: Figure of barren sites
# Save as ? RDS, shp, csv

# # Keep all species
# forCl <- as.data.frame( sppSiteList[1] ) 
# 
# # Remove rare species
# remSp <- as.data.frame( sppSiteList[2] ) 
# as.character( remSp$species )
# write_csv( remSp, paste0(region,"_SpeciesRemoved.csv" ))
# 
# # Remove barren sites
# remSites <- as.data.frame( sppSiteList[3] ) 
# colnames(remSites) <- c( "TransDepth","SpCnt" )
# remSites$TransDepth <- as.character(remSites$TransDepth)
# remSites <- left_join( remSites, matFull, by="TransDepth" )
# remSites[1,1:6] # confirm that position columns are named x,y
# layer <- paste0(region,"_SitesRemoved" )
# remSites <- makeShp( data=remSites, dsn=dsn, layer=layer )
# mapview(remSites, zcol="SpCnt")
# 
# saveRDS( sppSiteList, "sppSiteList.rds" ) 

rm(barrenSites,rareSpp,sppByRegion,sppDF)

# --- Run cluster analysis and build dendrogram --- #
# Create empty list to store results
myCluster <- vector("list", 4)

# Run cluster analysis
for (i in 1:length(prep4Cluster)){
  # Get region name
  names(myCluster)[[i]] <- names(prep4Cluster)[i]
  print(names(myCluster)[[i]])
  spp <- as.data.frame(prep4Cluster[[i]])
  # Calculate B-sim (simpson) distance on the site x species matrix
  myCluster[[i]]$dist <- sim( spp,  method=distance )
  # Create dendrogram, using average grouping
  myCluster[[i]]$benthtree <- hclust( myCluster[[i]]$dist, method=clusterMethod ) 
  #plot( myCluster[[i]]$benthtree, hang=-1 )
  # Calculate cophenetic correlation value
  # how well does the dendrogram correlate with the original data?
  myCluster[[i]]$cophenetic <- cor( myCluster[[i]]$dist, cophenetic(myCluster[[i]]$benthtree) )
  print(myCluster[[i]]$cophenetic)
}

# Plot and save dendrograms as png and as separate objects to play with cutoffs
for (i in 1:length(myCluster)){
  print(names(myCluster)[[i]])
  region <- names(myCluster)[[i]]
  file_name = paste("dendrogram_", region, ".png", sep="")
  title_name = paste(region,"Cluster Dendrogram", sep=" ")
  png(file_name)
  plot(myCluster[[i]]$benthtree, hang=-1, main=region)
  dev.off()
}

# # --- Manually choose heights to cut trees --- #
hg <- myCluster$HG$benthtree # seth = 4
ncc <- myCluster$NCC$benthtree # seth = 5
qcs <- myCluster$QCS$benthtree # seth = 2
sog <- myCluster$SoG$benthtree # seth = 2

benthtree <- sog
plot(benthtree, hang=-1)

# Choose cutoff, play with h until the visual clusters are kept together
# Set h (height) for the dendrogram and cut the tree
#seth <- readline(prompt = "Set height for the dendrogram and cut the tree: ")
seth <- 2
rect.hclust(benthtree, h=seth, border="red") # Cutoff based on visual inspection of the tree

# Vector to entry heights from plots
# hts <- c(4,5,2,2)

# Add height to large list
myCluster$HG$height <- 4 # Height on y-axis
myCluster$NCC$height <- 5
myCluster$QCS$height <- 2
myCluster$SoG$height <- 2

# Final h choice to slice tree 
for (i in 1:length(myCluster)){
  # Get region name
  print(names(myCluster)[[i]])
  myCluster[[i]]$sliceTree <- dendroextras::slice(x=myCluster[[i]]$benthtree, h=myCluster[[i]]$height)
}

# Build list of sliced trees
cl.list=list(list(c(myCluster$HG$sliceTree)),list(c(myCluster$NCC$sliceTree)),list(c(myCluster$QCS$sliceTree)),list(c(myCluster$SoG$sliceTree)))
names(cl.list) <- names(myCluster)

# --- Determine number of clusters to capture 90% of samples --- #
# Get table of cluster memberships - number of sites in each cluster
# lapply(list, function)

# Convert list element to a dataframe
make_df <- function(x){
  as.data.frame(table(x))
}
colorcount <- lapply(cl.list, make_df )
#names(colorcount) <- names(myCluster)
c.names <- c("cl","Freq")
colorcount <- lapply(colorcount, setNames, c.names)
colorcount

# Order clusters by frequency order 
order.cl <- function(x){
  order <- data.frame( order=seq( 1:nrow(x) ),x[order(-x$Freq),] ) 
  order$cumsum <- cumsum( order$Freq )
  order$cumpercent <- round( order$cumsum/max(order$cumsum), 2 )
  order$percent <- round( order$Freq/sum(order$Freq),2 )
  #plot( order$order, order$Freq, ylab="n Samples" )
  # plot( order$order, order$cumpercent )
  # abline( h=0.9, col="red" ) 
  return(order)
}
par(mfrow = c(2, 2))
cluster.frq <- lapply(colorcount, order.cl)
cluster.frq

# Number of clusters that capture 90% of samples
nTopClusters <- c(5,4,5,4)

# Could not get %in% to work for elements within a list, so here goes some ugly code!
hg <- cluster.frq$HG
ncc <- cluster.frq$NCC
qcs <- cluster.frq$QCS
sog <- cluster.frq$SoG

#HG
hg.cluster <- hg$cl[hg$order %in% c( 1:5 )] 
colorscheme <- myColors(5)
colorHG <- colorcount$HG
for (i in c(1:5)){
  colorHG$assigned[colorHG$cl==hg.cluster[i]] <- colorscheme[i]
}
table(colorHG$assigned)
colorHG$assigned[is.na(colorHG$assigned)] <- "grey"
colorHG
#NCC
ncc.cluster <- ncc$cl[ncc$order %in% c( 1:4 )] 
colorscheme <- myColors(4)
colorNCC <- colorcount$NCC
for (i in c(1:4)){
  colorNCC$assigned[colorNCC$cl==ncc.cluster[i]] <- colorscheme[i]
}
table(colorNCC$assigned)
colorNCC$assigned[is.na(colorNCC$assigned)] <- "grey"
colorNCC
#QCS
qcs.cluster <- qcs$cl[qcs$order %in% c( 1:5 )] 
colorscheme <- myColors(5)
colorQCS <- colorcount$QCS
for (i in c(1:5)){
  colorQCS$assigned[colorQCS$cl==qcs.cluster[i]] <- colorscheme[i]
}
table(colorQCS$assigned)
colorQCS$assigned[is.na(colorQCS$assigned)] <- "grey"
colorQCS
#SoG
sog.cluster <- sog$cl[sog$order %in% c( 1:4 )] 
colorscheme <- myColors(4)
colorSoG <- colorcount$SoG
for (i in c(1:4)){
  colorSoG$assigned[colorSoG$cl==sog.cluster[i]] <- colorscheme[i]
}
table(colorSoG$assigned)
colorSoG$assigned[is.na(colorSoG$assigned)] <- "grey"
colorSoG

# Put is all together
colorcount <- list(colorHG,colorNCC,colorQCS,colorSoG)
names(colorcount) <- c("HG","NCC","QCS","SoG")

# Add colour assignments back to myCluster list
# There must be a better way than this
myCluster$HG$colours <- colorcount[1]
myCluster$NCC$colours <- colorcount[2]
myCluster$QCS$colours <- colorcount[3]
myCluster$SoG$colours <- colorcount[4]

# --- Create color-coded dendrogram....slow... --- #
# lapply(list, function)
# branches <- function(x){
#   trees <- colour_branches( myCluster[[i]]$benthtree, h=myCluster[[i]]$height,col=as.character(myCluster[[i]]$colors$assigned) )
#   #plot(tree) 
# }
# 
# colour_tree <- lapply(myCluster, branches)
# 
# colortree <- colour_branches( benthtree, h=seth,col=as.character(colorcount$assigned) )
# plot(colortree) 

for (i in 1:length(myCluster)){
  # Get region name
  print(names(myCluster)[[i]])
  myCluster[[i]]$colourTree <- colour_branches( myCluster[[i]]$benthtree, h=myCluster[[i]]$height,col=as.character(myCluster[[i]]$colors$assigned) )
}







#===========================================#

# Extract richness values
richness <- map(summariesByRegion, "sppRichness")
str(richness)

# Flatten list into a dataframe to assess differences between regions
flat1 <- map(map_if(richness,~class(.x)=="matrix",list),~map(.x,as.data.frame))
flat2 <- map_dfr(flat1,~map_dfr(.x,identity,.id="TransDepth"),identity,.id="region")
colnames(flat2)[3] <- "spCnt"
head(flat2)

# Plot summaries
par(mfrow = c(1, 1))  # Set up a 1 x 1 plotting space
#boxplot(richness, ylab="Count of Species")
p <- ggplot(flat2, aes(x=region, y=spCnt)) +
  geom_boxplot(notch = T,
               fill = "cornflowerblue", 
               alpha = 0.7) +
  ylab("Count of Species") +
  xlab("Region") +
  scale_y_continuous(breaks=seq(0,60,by=10))+
  theme_bw()
p
ggsave("CountOfSpeciesByRegion.png", p)

# Summary table
nUnitsSum <- flat2 %>%
  group_by(region) %>%
  mutate(units = length(spCnt))%>%
  mutate(mean = mean(spCnt)) %>%
  mutate(min = min(spCnt)) %>%
  mutate(max = max(spCnt)) %>%
  select(region, units, mean, min, max) %>%
  distinct(region, units, mean, min, max)
nUnitsSum 

# Are there significant differences between regions?
kw <- kruskal.test(spCnt ~ region, data=flat2)
kw # very large chi-squared value!

# Which regions are different?
pair.w <- pairwise.wilcox.test(flat2$spCnt, flat2$region,
                     p.adjust.method = "BH") # What are the benefits of the different adjustment methods?
pair.w


# --- Species accumulation curves --- #
par(mfrow = c(2, 2))  # Set up a 2 x 2 plotting space

# Calculate summaries
for (i in 1:length(sppByRegion)){
  # Get region name
  names(summariesByRegion)[[i]] <- names(sppByRegion)[i]
  spp <- as.data.frame(sppByRegion[[i]])
  label <- names(summariesByRegion)[[i]]
  # Number of sampling units 
  sac <- specaccum(spp)
  plot(sac, main=label, ylim=c(0,170), ylab="Number of Species", ci.col="light blue")
  #plot(sac, main=label, ylim=c(0,170), xlim=c(0,2100), ylab="Number of Species", ci.col="light blue")
}

# --- Frequency of occurence of each species --- #
# Can be found in summariesByRegion[[i]]$sppFreq 

# Build one large table of species X survey area
foo <- vector("list", 4)

for (i in 1:length(summariesByRegion)){
  # Get region name
  names(foo)[[i]] <- names(summariesByRegion)[i]
  foo[[i]] <- summariesByRegion[[i]]$sppFreq
}

# Create large dataframe with containing all regions
foo1 <- bind_rows(foo,.id="Region")
head(foo1)

# Frequency of Occurrence table (Species Code | HG | NCC | QCS | SoG)
foo2 <- dplyr::select(foo1, Region, Species_Code, Freq)
foo2 <- foo2 %>%
  pivot_wider(names_from = Region,
              values_from = Freq) %>%
  arrange(Species_Code) 
head(foo2)

# Add species names and save
foo2 <- right_join(luTbl, foo2, by="Species_Code")
foo2 <- dplyr::select(foo2, Name, Species_Code, HG, NCC, QCS, SoG)

outfile <- paste(getwd(),"FreqOfOccurence_byRegion.csv",sep="/")
write_csv(foo2, outfile)

# --- Species Rank Plots --- #
rank <- dplyr::select(foo1, Region, Freq, Rank)

# Plot
p <- ggplot(rank, aes(x=Rank, y=Freq)) +
  geom_bar(stat="identity", color = "steelblue") +
  ylab("Frequency of Occurrence") +
  xlab("Rank of Species") +
  facet_wrap( ~ Region)
ggsave("RankOfSpecies.png", p)

# ...by invert and by algae
invert <- dplyr::filter(foo1, grepl("I_", Species_Code) )
p <- ggplot(invert, aes(x=Rank, y=Freq)) +
  geom_bar(stat="identity", color = "steelblue") +
  ylab("Frequency of Occurrence") +
  xlab("Rank of Invertebrate Species") +
  facet_wrap( ~ Region)
ggsave("InvertRanks.png", p)

algae <- dplyr::filter(foo1, grepl("A_", Species_Code) )
p <- ggplot(algae, aes(x=Rank, y=Freq)) +
  geom_bar(stat="identity", color = "steelblue") +
  ylab("Frequency of Occurrence") +
  xlab("Rank of Algae Species") +
  facet_wrap( ~ Region)
ggsave("AlgaeRanks.png", p)

# Top and bottom 20 species for each region
# top20 | HG.spCde | NCC.spCde | QCS.spCde | SoG.spCde
head(foo1)

# Select rows and columns of interest
top20 <- dplyr::filter(foo1, Rank<21 )
top20 <- dplyr::select(top20, Region, Species_Code, Rank)
# Create pivot table
top20 <- top20 %>%
  pivot_wider(names_from = Region,
              values_from = Species_Code)
head(top20)
# Save output
outfile <- paste(getwd(),"Top20_byRegion.csv",sep="/")
write_csv(top20, outfile)

# bottom20 | HG.spCde | NCC.spCde | QCS.spCde | SoG.spCde
# Select rows and columns of interest
bot20 <- dplyr::filter(foo1, Rank>149 )
bot20 <- dplyr::select(bot20, Region, Species_Code, Rank)
# Create pivot table
bot20 <- bot20 %>%
  pivot_wider(names_from = Region,
              values_from = Species_Code)
head(bot20)
# Save output
#outfile <- paste(getwd(),"Bottom20_byRegion.csv",sep="/")
#write_csv(bot20, outfile)


# --- Depth range summaries --- #
# Start with sppDF RDS
head(sppDF,3)
# Recreate the depth category field and remove the transect HKey
df <- separate(sppDF, TransDepth, c(NA,"DepthCat"))
head(df,3)
# Melt dataframe to long format
dCat <- df %>%
  pivot_longer(
    cols = -c(DepthCat, Region),
               names_to = "species",
               values_to = "presence")
# Remove absences
dCat <- dplyr::filter(dCat, presence!=0)
# Separate inverts and algae
algae.dc <- dplyr::filter(dCat, grepl("A_", species) )
# Remove duplicate records
algae.dc <- unique(algae.dc)
# Get rid of A_prefix
algae.dc <- separate(algae.dc, species, c(NA,"algae"))
head(algae.dc)
# Plot and save
p <- ggplot(algae.dc, aes(x=as.numeric(DepthCat),y=reorder(algae,desc(algae)))) +
        #geom_bar(stat = "identity",colour="blue", width=0.5) +
        geom_line(linetype="solid", colour="blue") +
        #geom_point(colour="blue") +
        #theme_bw() +
        facet_grid( ~ Region) +
        theme(panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed")) +
        xlab("Depth Category") +
        ylab("Algae code")
ggsave("AlgaeDepthCat_byRegion.png", p)  

# Repeat for invertebrates
invert.dc <- dplyr::filter(dCat, grepl("I_", species) )
invert.dc <- unique(invert.dc)
invert.dc <- separate(invert.dc, species, c(NA,"invert"))
head(invert.dc)
p <- ggplot(invert.dc, aes(x=as.numeric(DepthCat),y=reorder(invert,desc(invert)))) +
  geom_line(linetype="solid", colour="blue") +
  facet_grid( ~ Region ) +
  xlab("Depth Category") +
  ylab("Invertebrate code") +
  theme(panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed")) +
  theme(axis.text.y = element_text(face = "bold", size = 5))
ggsave("InvertDepthCat_byRegion.png", p)  




####################################################
################ Extra Junk ########################
####################################################

# fix plot, add jittered points

# p <- ggplot(richlist) + 
#   geom_boxplot(y=y)
# p
# 
# series <- seq(1:10)
# 
# lst <- list()
# 
# set.seed(28100)
# for (t in series) {
#   lst[[t]] <- sample(c(1:20, NA), sample(1:20, 1))
# }
# 
# # Modify lst into data frames of varying dimension
# lst <- lapply(series, function(x) {
#   data.frame(series = factor(x, levels = series),
#              y = lst[[x]])
# })
# 
# 
# # Stack the data frames
# lst <- do.call(rbind, lst)
# 
# 
# 
# # Make the plot
# ggplot(lst, aes(x = richness, y = y)) +
#   geom_boxplot()

# S <- specnumber(spp) # observed number of species
# raremax <- min(rowSums(t(spp)))
# Srare <- rarefy(t(spp), raremax)
# plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# abline(0, 1)
# rarecurve(t(spp), step = 2, sample = raremax, col = "blue", cex = 0.6)
 
# specpool(spp)
# s <- sample(nrow(spp), 500)
# specpool(spp[s,])
# 
# rarecurve(spp)
# plot(specaccum(spp))

# sp1 <- specaccum(spp)
# sp2 <- specaccum(spp, "random")
# sp2
# summary(sp2)
# plot(sp1)
# boxplot(sp1, col="yellow", add=TRUE, pch="+")

# 
# diversity(spp, index="shannon")
# # To Do: Add map of sites

# data(BCI)
# H <- diversity(BCI)
# simp <- diversity(BCI, "simpson")
# 

# 
# 
# 
# # Plot on a map
# #-------------
# # coords <- dplyr::select(env, LonStart, LatStart)
# # ggplot(coords, aes(LonStart, LatStart)) + geom_point()
# # filename <- paste0(outdir,"SiteLocations.pdf")
# # ggsave(filename)
# 
# # Combine regional summaries 
# #---------------------------
