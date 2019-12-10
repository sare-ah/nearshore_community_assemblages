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
luTbl <- read.csv( "T:/Benthic_Habitat_Mapping/Data/Look-upTbls/SpeciesLookUpTbl.csv", header=T, sep=",", stringsAsFactors=F )


# # Organize data into species & env
# #---------------------------------
# # Checks for complete recordset
# matFull <- unique(matFull)
# matFull <- dplyr::filter(matFull, x!=0)
# 
# # Set rownames
# rownames(matFull) <- matFull$TransDepth
# 
# # Divide columns into species and environmental variabls
# species <- names(matFull)[c(17:185)] # specify species columns, if different
# species
# idvars <- (matFull[,1]) # specify unique identifier for each row, if different
# 
# # Create matrix just of species columns
# spp <- matFull[,species]



# Regional summaries
#-------------------
# Biodiversity index *** To Do



output <- vector("list", 4)
for (i in 1:length(sppByRegion)){
  # Get region name
  names(output)[[i]] <- names(sppByRegion)[i]
  spp <- as.data.frame(sppByRegion[i])
  # Regional summaries
  #-------------------
  #
  # Site summaries
  #---------------
  # Number of sampling units 
  output[[i]]$nUnits <- nrow(spp)
  # Number of species per sampling unit
  output[[i]]$sppCnts <- as.data.frame(rowSums(spp))
  output[[i]]$freq <- spp.freq(speciesObs=spp)
  }


# colnames(sCnts) <- "sppCount"
# boxplot(sCnts$sppCount, horizontal = TRUE,main=c(region," Number of Species per Sampling Unit"))
# hist(sCnts$sppCount, prob=TRUE, main=c(region," Number of Species per Sampling Unit",ylab="Density",xlab="Counts"))




siteSummaries <- function(sites=sppByRegion[i]){
  nUnits <- nrow(sites)
  # Number of species per sampling unit
  sCnts <- as.data.frame(rowSums(sites))
  result <- list(nUnits,sCnts)
}
  
  


# To Do: Add map of sites


# Species summaries
#------------------
sppFinal <- spp.freq(speciesObs=spp, speciesLut=luTbl, samplingUnits=nUnits)

# Top and bottom 20 species
top20 <- sppFinal %>%
  arrange(desc(Freq)) %>%
  slice(1:20)
top20
sppFinal <- dplyr::filter(sppFinal, Freq!=0)
a <- nrow(sppFinal)
b <- a -19
bottom20 <- sppFinal %>%
  arrange(desc(Freq)) %>%
  slice(b:a)
bottom20


# Depth range summaries
#----------------------
# Subset 
dCat <- dplyr::select(matFull, c(13,17:185))

# Melt dataframe
dCat <- dRng %>%
  pivot_longer(-DepthCat,
               names_to = "species",
               values_to = "presence")

dCat <- dplyr::filter(dCat, presence!=0)

# Plot and save
# Plot is incorrect! Should be a histogram with species along x-axis & categories along the y-axis
# May want to separate by inverts and algae for easier reading
ggplot(dCat, aes(x=DepthCat,y=reorder(species,desc(species))))+
  geom_line(colour="blue") +
  theme_bw() +
  scale_x_continuous(breaks=seq(-8,20,by=1))+
  labs(x="Depth Category")+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(face="italic"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"))
filename <- paste0(region,".DepthCategories")
ggsave(file=paste0(outdir,filename,'.pdf',sep=''))


# Plot on a map
#-------------
# coords <- dplyr::select(env, LonStart, LatStart)
# ggplot(coords, aes(LonStart, LatStart)) + geom_point()
# filename <- paste0(outdir,"SiteLocations.pdf")
# ggsave(filename)

# Combine regional summaries 
#---------------------------
