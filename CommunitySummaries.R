#######################################################################################################################
# Cluster analysis of Benthic Habitat Mapping Dive Survey Sites
# Summaries


# Start fresh
rm(list=ls())

# Load packages
#--------------
library(mapview)
library(tidyverse)
library(rstudioapi)
library(vegan)

# Inputs
#-------

date <- format(Sys.Date(), "%b_%d")
region <- "All"
outdir <- paste0(date,"_",region)
dsn <- "SHP"

# Site by Species csv file
# myFile <- choose.files(default = "T:/Benthic_Habitat_Mapping/Data",
#                        caption = "Select site X species, with enviromental variables") 
myFile <- "T:/Benthic_Habitat_Mapping/Data/Species by Site Matrices/qcsbyDepthCat_AllSpp.csv"

# Read in species by regions
#---------------------------
sppByRegion <- readRDS("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/RDS/sppByRegion.rds")

# Get path for this script
#-------------------------
# Set working directory to one above script directory
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# Source functions
#-----------------
source('CommunitySummaries_functions.R')

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
#luTbl <- read.csv( "T:/Benthic_Habitat_Mapping/Data/Look-upTbls/SpeciesLookUpTbl.csv", header=T, sep=",", stringsAsFactors=F )

# Summaries
#----------
# To Do: Biodiversity indices *** 

summariesByRegion <- vector("list", 4)
# Calculate summaries
for (i in 1:length(sppByRegion)){
  # Get region name
  names(summariesByRegion)[[i]] <- names(sppByRegion)[i]
  spp <- as.data.frame(sppByRegion[[i]])
  # Number of sampling units 
  summariesByRegion[[i]]$nUnits <- nrow(spp)
  # Number of species per sampling unit (Richness)
  summariesByRegion[[i]]$sppRichness <- rowSums(spp)
  # Frequency of occurence for each species
  summariesByRegion[[i]]$sppFreq <- spp.freq(speciesObs=spp)
  }

# Plot summaries
par(mfrow = c(1, 1))  # Set up a 1 x 1 plotting space

# Extract richness values
richness <- map(summariesByRegion, "sppRichness")

boxplot(richness, main="Number of Species per Sampling Unit")


# #hist(richness, prob=TRUE, main=" Number of Species per Sampling Unit",ylab="Density",xlab="Counts")
# 
# 
# 
# spp <- as.data.frame(sppByRegion[[1]])
# 
# S <- specnumber(spp) # observed number of species
# raremax <- min(rowSums(t(spp)))
# Srare <- rarefy(t(spp), raremax)
# plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# abline(0, 1)
# rarecurve(t(spp), step = 2, sample = raremax, col = "blue", cex = 0.6)
# 
# 
# 
# specpool(spp)
# s <- sample(nrow(spp), 500)
# specpool(spp[s,])
# 
# rarecurve(spp)
# plot(specaccum(spp))
# 
# # Species accumulation curve
# sac <- specaccum(spp)
# plot(sac, ci.type="polygon", ci.col="light blue")  
# 
# 
# sp1 <- specaccum(spp)
# sp2 <- specaccum(spp, "random")
# sp2
# summary(sp2)
# plot(sp1)
# boxplot(sp1, col="yellow", add=TRUE, pch="+")
# 
# 
# 
# diversity(spp, index="shannon")
# # To Do: Add map of sites
# 
# data(BCI)
# H <- diversity(BCI)
# simp <- diversity(BCI, "simpson")
# 
# # Species summaries
# #------------------
# #sppFinal <- spp.freq(speciesObs=spp, speciesLut=luTbl, samplingUnits=nUnits)
# 
# # Top and bottom 20 species
# top20 <- sppFinal %>%
#   arrange(desc(Freq)) %>%
#   slice(1:20)
# top20
# sppFinal <- dplyr::filter(sppFinal, Freq!=0)
# a <- nrow(sppFinal)
# b <- a -19
# bottom20 <- sppFinal %>%
#   arrange(desc(Freq)) %>%
#   slice(b:a)
# bottom20
# 
# 
# # Depth range summaries
# #----------------------
# # Subset 
# dCat <- dplyr::select(matFull, c(13,17:185))
# 
# # Melt dataframe
# dCat <- dRng %>%
#   pivot_longer(-DepthCat,
#                names_to = "species",
#                values_to = "presence")
# 
# dCat <- dplyr::filter(dCat, presence!=0)
# 
# # Plot and save
# # Plot is incorrect! Should be a histogram with species along x-axis & categories along the y-axis
# # May want to separate by inverts and algae for easier reading
# ggplot(dCat, aes(x=DepthCat,y=reorder(species,desc(species))))+
#   geom_line(colour="blue") +
#   theme_bw() +
#   scale_x_continuous(breaks=seq(-8,20,by=1))+
#   labs(x="Depth Category")+
#   theme(axis.title.y=element_blank(),
#         axis.text.y = element_text(face="italic"),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"))
# filename <- paste0(region,".DepthCategories")
# ggsave(file=paste0(outdir,filename,'.pdf',sep=''))
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
