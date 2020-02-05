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
# myFile <- "T:/Benthic_Habitat_Mapping/Data/Species by Site Matrices/qcsbyDepthCat_AllSpp.csv"

# Read in species by regions
#---------------------------
sppByRegion <- readRDS("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/RDS/sppByRegion.AllBC.rds")
sppDF <- readRDS("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/RDS/sppByRegion_Dataframe.rds")

nElements <- length(sppByRegion)

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

# Summaries
#----------
# To Do: Biodiversity indices *** 


# --- Build list of summary stats --- #
summariesByRegion <- vector("list", nElements)

# Calculate summaries
for (i in 1:length(sppByRegion)){
  # Get region name
  names(summariesByRegion)[[i]] <- names(sppByRegion)[i]
  spp <- as.data.frame(sppByRegion[[i]])
  # Number of sampling units 
  summariesByRegion[[i]]$nUnits <- nrow(spp)
  # Number of species per sampling unit (Richness)
  summariesByRegion[[i]]$sppRichness <- rowSums(spp)
  #print(min(rowSums(spp)))
  # Frequency of occurence for each species
  summariesByRegion[[i]]$sppFreq <- spp.freq(speciesObs=spp)
  }


# --- Number of species per sampling unit --- #

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

# Summary table - Sampling units and mean number of species observed per sampling area
nUnitsSum <- flat2 %>%
  group_by(region) %>%
  mutate(units = length(spCnt))%>%
  mutate(mean.SpCnt = mean(spCnt)) %>%
  mutate(min.SpCnt = min(spCnt)) %>%
  mutate(max.SpCnt = max(spCnt)) %>%
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
par(mfrow = c(3, 2))  # Set up a 2 x 2 plotting space

# Calculate summaries --- add save as png file
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
foo <- vector("list", nElements)

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
foo2 <- dplyr::select(foo2, Name, Species_Code, HG, NCC, QCS, SoG, ALL)

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
sppDF$TransDepth <- row.names(sppDF)
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
